#include <litho_invert/inversion/runner.h>
#include <litho_invert/inversion/gncg_optimizer.h>
#include <litho_invert/core/stats.h>
#include <filesystem>
#include <cmath>
#include <iostream>

namespace litho_invert {

InversionRunner::InversionRunner(const InversionConfig& config)
    : m_config(config)
{
    static const bool _banner = []() {
        std::cerr << "================================================\n"
                  << "LithoInv v1.0 -- Lithology-Constrained Joint Inversion\n"
                  << "Copyright (c) 2026 Thomas Bagley. MIT Licensed.\n"
                  << "Citation: Thomas Bagley (2026), preprint on EarthArXiv.\n"
                  << "================================================" << std::endl;
        return true;
    }();
    (void)_banner;
    setup();
}

void InversionRunner::setIterationCallback(IterationCallback cb) {
    m_callback = std::move(cb);
}

void InversionRunner::setup() {
    // setup() called from InversionRunner constructor
    // 1. Apply control-point downsampling if configured (must be done
    //    before building neighbors, which depends on the full mesh).
    std::cerr << "[DEBUG step 1]" << std::endl;
    if (m_config.controlPointStride > 0) {
        uint32_t coarseDofs = m_config.model->setControlPointStride(
            m_config.controlPointStride);
        std::cerr << "[DEBUG setControlPointStride returned " << coarseDofs << "]" << std::endl;
        if (coarseDofs == 0) {
            // Grid structure not suitable; fall back to full resolution
            m_config.controlPointStride = 0;
        }
    }
    // Debug: print dof count per surface
    for (int i = 0; i < m_config.model->groupMeshCount(); ++i) {
        auto* surf = m_config.model->groupMesh(i);
        std::cerr << "[DEBUG surface " << i << ": vertices=" << surf->vertexCount()
                  << " dofs=" << surf->dofCount()
                  << " interior=" << surf->interiorGridDim()
                  << " padding=" << surf->paddingRings()
                  << " stride=" << m_config.controlPointStride
                  << "]" << std::endl;
    }
    std::cerr << "[DEBUG totalDofCount=" << m_config.model->totalDofCount() << "]" << std::endl;

    // 2. Build neighbors on all surfaces (for regularization)
    std::cerr << "[DEBUG step 2: buildNeighbors, surfaces=" << m_config.model->groupMeshCount() << "]" << std::endl;
    for (int i = 0; i < m_config.model->groupMeshCount(); ++i) {
        m_config.model->groupMesh(i)->buildNeighbors();
    }

    // 2b. Set per-axis parameter scaling for well-conditioned optimization.
    //     Uniform axis scaling: 1 optimizer unit ≈ 10m in any direction.
    //     The physical gravity gradient is naturally larger in Z, so the
    //     optimizer will prioritize vertical movement without artificially
    //     amplifying XY.  This prevents lateral model expansion during
    //     early iterations when the smoothness gradient is near zero.
    for (int i = 0; i < m_config.model->groupMeshCount(); ++i) {
        auto* surf = m_config.model->groupMesh(i);
        surf->setAxisScale(0, 10.0);
        surf->setAxisScale(1, 10.0);
        surf->setAxisScale(2, 10.0);
    }

    // 3. Create GravityForward with model and observed data
    std::cerr << "[DEBUG step 3: GravityForward]" << std::endl;
    m_forward = std::make_shared<GravityForward>(m_config.model, m_config.observedData);

    // 3a. Scale FD step for the parameter scaling in use.
    //     Default 0.01 is tuned for unscaled Z_ONLY. With XYZ_FREE scaling,
    //     a unit step in parameter space corresponds to 50m XY / 10m Z.
    //     Use a larger FD step so the perturbation is physically meaningful.
    double fdStep = m_config.fdStep;
    m_forward->setFiniteDifferenceStep(fdStep);

    // 3b. ObjectiveFunction with forward and data
    std::cerr << "[DEBUG step 4: ObjectiveFunction]" << std::endl;

    // Apply default per-data uncertainties when data points have σ=0.
    // Without this, unweighted gravity (σ=0→w=1) and magnetics (σ=0→w=1)
    // produce raw misfit scales that differ by 10³–10⁶, causing one data
    // type to dominate the inversion entirely.
    if (m_config.gravityUncertainty > 0.0) {
        for (auto& pt : m_config.observedData) {
            if (pt.g_std <= 0.0) pt.g_std = m_config.gravityUncertainty;
        }
    }
    if (m_config.magneticUncertainty > 0.0) {
        for (auto& pt : m_config.magneticData) {
            if (pt.t_std <= 0.0) pt.t_std = m_config.magneticUncertainty;
        }
    }

    m_objective = std::make_shared<ObjectiveFunction>(m_forward, m_config.observedData);

    // 5. Create MagneticForward if magnetic data present, then JointObjective
    std::cerr << "[DEBUG step 5: magnetic=" << (!m_config.magneticData.empty() ? "yes" : "no") << "]" << std::endl;
    if (!m_config.magneticData.empty()) {
        m_magneticForward = std::make_shared<MagneticForward>(
            m_config.model, m_config.magneticData,
            m_config.magneticInclination, m_config.magneticDeclination,
            m_config.magneticField_nT);

        m_magneticForward->setRemanenceMode(m_config.remanenceMode);
        m_magneticForward->setFiniteDifferenceStep(fdStep);

        m_jointObjective = std::make_shared<JointObjective>(
            m_objective, m_magneticForward, m_config.magneticData,
            m_config.magneticWeight);
    }

    // 5b. Create EMActiveForward if active EM data present
    if (!m_config.activeEMData.empty()) {
        m_activeEMForward = std::make_shared<EMActiveForward>(
            m_config.model, m_config.emSources, m_config.emReceivers,
            m_config.activeEMData, m_config.emConfig);

        if (m_config.enablePaddingGroup) {
            m_activeEMForward->enablePadding(true, m_config.paddingDepth);
            m_activeEMForward->setPaddingConductivity(m_config.paddingConductivityInitial);
        }
    }

    // 5c. Create EMMTForward if MT data present
    if (!m_config.mtData.empty()) {
        m_mtForward = std::make_shared<EMMTForward>(
            m_config.model, m_config.mtStations,
            m_config.mtData, m_config.emConfig);

        if (m_config.enablePaddingGroup) {
            m_mtForward->enablePadding(true, m_config.paddingDepth);
            m_mtForward->setPaddingConductivity(m_config.paddingConductivityInitial);
        }
    }

    // 5d. If EM data present but no magnetic, still need JointObjective for EM
    if (m_jointObjective == nullptr &&
        (!m_config.activeEMData.empty() || !m_config.mtData.empty())) {
        // Create a minimal JointObjective with just gravity, then add EM terms
        m_jointObjective = std::make_shared<JointObjective>(
            m_objective, nullptr, MagneticData{}, 0.0);
    }

    // 5e. Add EM objectives to the joint objective if present
    if (m_jointObjective && m_activeEMForward) {
        m_jointObjective->addActiveEM(m_activeEMForward, m_config.activeEMData,
                                      m_config.activeEMWeight);
    }
    if (m_jointObjective && m_mtForward) {
        m_jointObjective->addMT(m_mtForward, m_config.mtData,
                                m_config.mtWeight);
    }

    // 6. Create SurfaceSmoothness with model, set weight = lambda
    std::cerr << "[DEBUG step 6: SurfaceSmoothness]" << std::endl;
    m_smoothness = std::make_shared<SurfaceSmoothness>(m_config.model);
    m_smoothness->setWeight(m_config.lambda);

    // 7. Add smoothness to objective
    std::cerr << "[DEBUG step 7: addRegularization]" << std::endl;
    m_objective->addRegularization(m_smoothness);

    // 7b. Reference model regularization (optional)
    if (m_config.enableReferenceModel) {
        std::cerr << "[DEBUG step 7b: ReferenceModelRegularization, lambda_ref="
                  << m_config.lambdaRef << "]" << std::endl;
        m_referenceReg = std::make_shared<ReferenceModelRegularization>(m_config.model);
        m_referenceReg->setWeight(m_config.lambdaRef);
        m_objective->addRegularization(m_referenceReg);
    }

    // 8. If constraints not empty, create ConstraintHandler
    std::cerr << "[DEBUG step 8: constraints=" << m_config.constraints.size() << "]" << std::endl;
    if (!m_config.constraints.empty()) {
        m_constraintHandler = std::make_shared<ConstraintHandler>(
            m_config.model, m_config.constraints);
        m_constraintHandler->setWeight(m_config.omega);
        m_objective->setConstraintHandler(m_constraintHandler);
    }

    // 9. Create geometry optimizer based on solver selection
    if (m_config.solver == "gncg") {
        std::cerr << "[DEBUG step 9: GNCGOptimizer]" << std::endl;
        auto gncg = std::make_shared<GNCGOptimizer>();
        gncg->setMaxIterations(m_config.maxIterations);
        gncg->setCGMaxIter(m_config.gncgCGMaxIter);
        gncg->setCGTolerance(m_config.gncgCGTolerance);
        gncg->setFDStep(m_config.fdStep);
        gncg->setArmijoC1(m_config.armijoC1);
        gncg->setLineSearchMaxIter(m_config.lineSearchMaxIter);
        gncg->setTolerance(m_config.tolerance);
        m_geometryOptimizer = gncg;
    } else {
        std::cerr << "[DEBUG step 9: LBFGSBOptimizer]" << std::endl;
        auto lbfgsb = std::make_shared<LBFGSBOptimizer>();
        lbfgsb->setHistorySize(m_config.lbfgsHistory);
        lbfgsb->setTolerance(m_config.tolerance);
        lbfgsb->setArmijoC1(m_config.armijoC1);
        lbfgsb->setLineSearchMaxIter(m_config.lineSearchMaxIter);
        lbfgsb->setDisableLineSearch(m_config.disableLineSearch);
        m_geometryOptimizer = lbfgsb;
    }

    std::cerr << "[DEBUG setup() done]" << std::endl;
}

InversionResult InversionRunner::run() {
    InversionResult result;
    result.finalModel = m_config.model;

    // 1. Assemble initial parameter vector from model surfaces
    VectorXd params = m_config.model->assembleParameterVector();
    size_t nData = m_config.observedData.size();

    // 2. Get geometry bounds
    VectorXd lower, upper;
    m_config.model->getBounds(lower, upper);

    // 2b. Capture reference for reference model regularization
    if (m_referenceReg) {
        if (m_config.referenceParams.size() > 0) {
            m_referenceReg->captureReference(m_config.referenceParams);
            std::cerr << "[DEBUG] Reference model set from config (" << m_config.referenceParams.size() << " DOFs)" << std::endl;
        } else {
            m_referenceReg->captureReference(params);
            std::cerr << "[DEBUG] Reference model captured (" << params.size() << " DOFs)" << std::endl;
        }
    }

    // 2c. Apply hard depth bounds if enabled
    if (m_config.enableDepthBounds) {
        double margin = m_config.depthBoundMargin;
        for (Index i = 0; i < params.size(); ++i) {
            lower[i] = std::max(lower[i], params[i] - margin);
            upper[i] = std::min(upper[i], params[i] + margin);
        }
        std::cerr << "[DEBUG] Depth bounds applied: +/- " << margin << " from starting position" << std::endl;
    }

    // 3. Record initial iteration (iteration 0)
    bool hasMagnetics = (m_jointObjective != nullptr && m_magneticForward != nullptr);
    bool hasActiveEM = (m_activeEMForward != nullptr);
    bool hasMT = (m_mtForward != nullptr);
    bool hasJoint = (m_jointObjective != nullptr);
    {
        InversionIteration initIter;
        initIter.iteration = 0;
        if (hasMagnetics) {
            auto comp = m_jointObjective->evaluateComponents(params);
            initIter.dataMisfit = comp.gravityMisfit;
            initIter.magneticMisfit = comp.magneticMisfit;
            initIter.regularization = comp.regularization;
            initIter.constraintPenalty = comp.constraintPenalty;
            initIter.totalObjective = comp.total;
            // Compute unweighted RMS from residuals
            {
                VectorXd gResid = m_jointObjective->gravityResiduals(params);
                initIter.rmsError = std::sqrt(gResid.squaredNorm() / static_cast<double>(nData));
                std::vector<Vector3d> gPositions;
                gPositions.reserve(nData);
                for (size_t i = 0; i < nData; ++i)
                    gPositions.push_back(m_config.observedData[i].position);
                auto dw = computeDurbinWatson(gPositions, gResid);
                initIter.dw_gravity_x = dw.dw_x;
                initIter.dw_gravity_y = dw.dw_y;
            }
        } else {
            auto components = m_objective->evaluateComponents(params);
            initIter.dataMisfit = components.dataMisfit;
            initIter.regularization = components.regularization;
            initIter.constraintPenalty = components.constraintPenalty;
            initIter.totalObjective = components.total;
            // Compute unweighted RMS and DW from residuals
            {
                VectorXd gResid = m_objective->residuals(params);
                initIter.rmsError = std::sqrt(gResid.squaredNorm() / static_cast<double>(nData));
                std::vector<Vector3d> gPositions;
                gPositions.reserve(nData);
                for (size_t i = 0; i < nData; ++i)
                    gPositions.push_back(m_config.observedData[i].position);
                auto dw = computeDurbinWatson(gPositions, gResid);
                initIter.dw_gravity_x = dw.dw_x;
                initIter.dw_gravity_y = dw.dw_y;
            }
        }
        result.history.push_back(initIter);
        if (m_callback) m_callback(initIter);

        // Pass initial objective value to L-BFGS-B for scale-aware convergence.
        if (m_config.solver != "gncg") {
            auto lbfgs = std::dynamic_pointer_cast<LBFGSBOptimizer>(m_geometryOptimizer);
            if (lbfgs) lbfgs->setReferenceObjective(initIter.totalObjective);
        }
    }

    // Eigenvalue-based Hessian scaling: after the first Jacobian evaluation
    // at the starting model, scale the magnetic weight so both data types
    // contribute comparable gradient magnitudes.
    if (hasMagnetics && m_config.enableEigenvalueScaling
        && m_config.magneticWeight > 0.0) {
        m_jointObjective->calibrateDataWeights(params);
    }

    // 4. Define geometry objective and gradient wrappers
    // Track last recorded parameters to deduplicate line-search evaluations
    // (L-BFGS line search calls objective several times per iteration; we only
    // want one history entry per actual geometry update).

    // 4a. Per-iteration export helper (exports all group meshes as .ts)
    auto exportIteration = [this]() {
        if (m_config.iterationExportDir.empty()) return;

        // Per-iteration subdirectory to keep line-search trial points
        // separate from validated geometry exports.
        char subdir[16];
        snprintf(subdir, sizeof(subdir), "iter_%03d", m_iterationCounter++);
        std::string iterDir = m_config.iterationExportDir + "/" + subdir;

        std::error_code ec;
        std::filesystem::create_directories(iterDir, ec);

        InversionExporter exporter(iterDir, "invert");
        if (!m_config.groupExportNames.empty()) {
            exporter.setGroupNaming(m_config.groupExportNames);
        }
        for (int g = 0; g < m_config.model->groupMeshCount(); ++g) {
            const auto* mesh = m_config.model->groupMesh(g);
            if (!mesh || mesh->vertexCount() == 0) continue;
            std::string name;
            if (!m_config.groupExportNames.empty()
                    && g < static_cast<int>(m_config.groupExportNames.size())) {
                name = m_config.groupExportNames[g];
            } else {
                name = "group_" + std::to_string(g);
            }
            exporter.exportClosedVolume(*mesh, name);
        }
    };

    VectorXd lastRecordedParams = params;

    auto geoObj = [this, &result, nData, hasMagnetics, hasActiveEM, hasMT,
                   hasJoint, &lastRecordedParams, &exportIteration](const VectorXd& p) -> double {
        // Line-search diagnostic export: snapshot geometry at every trial point
        if (!m_config.lineSearchExportDir.empty()) {
            char trialDir[32];
            snprintf(trialDir, sizeof(trialDir), "iter_%03d/trial_%03d",
                     m_iterationCounter, m_lineSearchTrialCounter++);
            std::string path = m_config.lineSearchExportDir + "/" + trialDir;
            std::error_code ec;
            std::filesystem::create_directories(path, ec);
            InversionExporter exporter(path, "ls");
            if (!m_config.groupExportNames.empty())
                exporter.setGroupNaming(m_config.groupExportNames);
            for (int g = 0; g < m_config.model->groupMeshCount(); ++g) {
                const auto* mesh = m_config.model->groupMesh(g);
                if (!mesh || mesh->vertexCount() == 0) continue;
                std::string name;
                if (!m_config.groupExportNames.empty()
                        && g < static_cast<int>(m_config.groupExportNames.size()))
                    name = m_config.groupExportNames[g];
                else
                    name = "group_" + std::to_string(g);
                exporter.exportClosedVolume(*mesh, name);
            }
        }

        InversionIteration iter;
        iter.iteration = static_cast<int>(result.history.size());

        if (hasJoint) {
            auto comp = m_jointObjective->evaluateComponents(p);
            iter.dataMisfit = comp.gravityMisfit;
            iter.magneticMisfit = comp.magneticMisfit;
            iter.regularization = comp.regularization;
            iter.constraintPenalty = comp.constraintPenalty;
            iter.totalObjective = comp.total;
            // Compute unweighted RMS and DW from gravity residuals
            {
                VectorXd gResid = m_jointObjective->gravityResiduals(p);
                iter.rmsError = std::sqrt(gResid.squaredNorm() / static_cast<double>(nData));
                std::vector<Vector3d> gPositions;
                gPositions.reserve(nData);
                for (size_t i = 0; i < nData; ++i)
                    gPositions.push_back(m_config.observedData[i].position);
                auto dw = computeDurbinWatson(gPositions, gResid);
                iter.dw_gravity_x = dw.dw_x;
                iter.dw_gravity_y = dw.dw_y;
            }

            // Compute Durbin-Watson on magnetic residuals
            {
                VectorXd mResid = m_jointObjective->magneticResiduals(p);
                size_t nM = m_config.magneticData.size();
                std::vector<Vector3d> mPositions;
                mPositions.reserve(nM);
                for (size_t i = 0; i < nM; ++i) {
                    mPositions.push_back(m_config.magneticData[i].position);
                }
                auto dw = computeDurbinWatson(mPositions, mResid);
                iter.dw_magnetic_x = dw.dw_x;
                iter.dw_magnetic_y = dw.dw_y;
            }

            // Active EM misfit + DW
            if (hasActiveEM) {
                iter.activeEMMisfit = comp.activeEMMisfit;
                VectorXd emResid = m_jointObjective->activeEMResiduals(p);
                size_t nEM = m_config.activeEMData.size();
                std::vector<Vector3d> emPositions;
                emPositions.reserve(nEM);
                for (size_t i = 0; i < nEM; ++i) {
                    emPositions.push_back(m_config.activeEMData[i].position);
                }
                auto dw = computeDurbinWatson(emPositions, emResid);
                iter.dw_activeEM_x = dw.dw_x;
                iter.dw_activeEM_y = dw.dw_y;
            }

            // MT misfit + DW
            if (hasMT) {
                iter.mtMisfit = comp.mtMisfit;
                VectorXd mtResid = m_jointObjective->mtResiduals(p);
                std::vector<Vector3d> mtPositions;
                for (size_t i = 0; i < m_config.mtStations.size(); ++i) {
                    mtPositions.push_back(m_config.mtStations[i].position);
                }
                auto dw = computeDurbinWatson(mtPositions, mtResid);
                iter.dw_mt_x = dw.dw_x;
                iter.dw_mt_y = dw.dw_y;
            }

            // Only record history when parameters actually change (dedup line-search calls)
            if (result.history.empty() || (p - lastRecordedParams).squaredNorm() > 1e-14) {
                lastRecordedParams = p;
                result.history.push_back(iter);
                if (m_callback) m_callback(iter);
            }
            return comp.total;
        } else {
            auto components = m_objective->evaluateComponents(p);
            iter.dataMisfit = components.dataMisfit;
            iter.regularization = components.regularization;
            iter.constraintPenalty = components.constraintPenalty;
            iter.totalObjective = components.total;
            // Compute unweighted RMS and DW from gravity residuals
            {
                VectorXd gResid = m_objective->residuals(p);
                iter.rmsError = std::sqrt(gResid.squaredNorm() / static_cast<double>(nData));
                std::vector<Vector3d> gPositions;
                gPositions.reserve(nData);
                for (size_t i = 0; i < nData; ++i)
                    gPositions.push_back(m_config.observedData[i].position);
                auto dw = computeDurbinWatson(gPositions, gResid);
                iter.dw_gravity_x = dw.dw_x;
                iter.dw_gravity_y = dw.dw_y;
            }

            if (result.history.empty() || (p - lastRecordedParams).squaredNorm() > 1e-14) {
                lastRecordedParams = p;
                result.history.push_back(iter);
                if (m_callback) m_callback(iter);
            }
            return components.total;
        }
    };

    auto geoGrad = [this, hasJoint](const VectorXd& p) -> VectorXd {
        if (hasJoint) {
            auto g = m_jointObjective->gradient(p);
            g = smoothGradient(g);
            return g;
        }
        auto g = m_objective->gradient(p);
        g = smoothGradient(g);
        return g;
    };

    // 5. Main phased loop: alternate geometry and property inversion
    int totalIters = 0;
    bool converged = false;

    // Compute overlap baselines for validation gates. Discretization of
    // shared contact surfaces produces small baseline counts that the
    // quality gate subtracts so it only flags regressions.
    size_t overlapBaseline = 0;
    size_t adjacentOverlapBaseline = 0;
    {
        int nM = m_config.model->groupMeshCount();
        for (int g = 0; g < nM; ++g) {
            const auto* mesh = m_config.model->groupMesh(g);
            if (!mesh) continue;
            for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
                int c = m_config.model->classifyPoint(mesh->vertex(vi).position);
                if (c < 0 || c == g) continue;
                if (std::abs(c - g) == 1)
                    ++adjacentOverlapBaseline;
                else
                    ++overlapBaseline;
            }
        }
        if (overlapBaseline > 0 || adjacentOverlapBaseline > 0) {
            std::cerr << "[DEBUG overlap baselines: non-adjacent=" << overlapBaseline
                      << " adjacent=" << adjacentOverlapBaseline << "]" << std::endl;
        }
    }

    while (totalIters < m_config.maxIterations && !converged
           && static_cast<int>(result.history.size()) < m_config.maxIterations * 3) {

        if (m_config.enableGeometryInversion) {
            // --- Phase A: Geometry optimization ---
            int remaining = m_config.maxIterations - totalIters;
            int phaseIters = std::min(m_config.propertyInversionInterval, remaining);

            if (phaseIters <= 0) break;

            m_geometryOptimizer->setMaxIterations(phaseIters);

            m_lineSearchTrialCounter = 0;
            std::cerr << "[DEBUG calling minimize, phaseIters=" << phaseIters << " params=" << params.size() << "]" << std::endl;
            auto geoResult = m_geometryOptimizer->minimize(
                geoObj, geoGrad, params, lower, upper);
            std::cerr << "[DEBUG minimize returned, iters=" << geoResult.iterations << " converged=" << geoResult.converged << "]" << std::endl;

            VectorXd prevParams = params;

            // Apply optimized geometry with validation gate.
            // When the step produces broken geometry (gaps, spikes,
            // edge tearing), dampen and retry.
            VectorXd bestParams = geoResult.params;
            bool accepted = false;
            const int maxRetries = 5;
            for (int retry = 0; retry <= maxRetries; ++retry) {
                VectorXd trialParams;
                if (retry == 0) {
                    trialParams = geoResult.params;
                } else {
                    // Dampen: reduce the step toward the previous params
                    double scale = 1.0 / (1 << retry);  // 1/2, 1/4, 1/8, ...
                    trialParams = prevParams + scale * (geoResult.params - prevParams);
                }
                m_config.model->applyParameterVector(trialParams);

                auto validation = m_config.model->validate(overlapBaseline,
                                                             adjacentOverlapBaseline);
                if (validation.passed) {
                    bestParams = trialParams;
                    accepted = true;
                    break;
                }
                std::cerr << "[WARNING] Model validation failed (retry " << retry
                          << "): " << validation.failureReason
                          << " (violation=" << validation.worstViolation << ")" << std::endl;
            }

            if (!accepted) {
                // All retries failed — revert to previous safe state
                m_config.model->applyParameterVector(prevParams);
                std::cerr << "[WARNING] All validation retries failed, reverting to "
                          << "previous model state" << std::endl;
            } else {
                params = bestParams;
                exportIteration();
            }

            if (geoResult.iterations == 0) {
                converged = true;
                break;
            }

            if (accepted) {
                converged = geoResult.converged;
                totalIters += geoResult.iterations;
            }
        } else {
            // Property-only mode: skip geometry, just do property inversion
            converged = true; // no geometry to converge
        }

        // --- Phase B: Property inversion (always run when enabled) ---
        if (m_config.enablePropertyInversion
                && totalIters < m_config.maxIterations) {
            runPropertyInversion(result);
            if (!m_config.enableGeometryInversion) {
                // In property-only mode, run one round and stop
                converged = true;
            } else {
                converged = false;
            }
        }
    }

    // 6. Compute final predicted data (gravity)
    result.predictedData = m_forward->compute(params);

    // 7. Compute final misfit and RMS
    result.converged = converged;
    result.totalIterations = totalIters;

    double rmsSum = 0.0;
    for (size_t i = 0; i < nData; ++i) {
        double res = m_config.observedData[i].g_obs - result.predictedData[static_cast<Index>(i)];
        rmsSum += res * res;
    }
    result.finalRMS = std::sqrt(rmsSum / static_cast<double>(nData));

    if (hasMagnetics) {
        auto comp = m_jointObjective->evaluateComponents(params);
        result.finalMisfit = comp.total;
    } else {
        auto finalComponents = m_objective->evaluateComponents(params);
        result.finalMisfit = finalComponents.total;
    }

    // 8. Record final densities, susceptibilities, conductivities, and remanence
    int nGroups = m_config.model->groupCount();
    result.finalDensities.resize(nGroups);
    result.finalSusceptibilities.resize(nGroups);
    result.finalConductivities.resize(nGroups);
    result.finalRemanences.resize(nGroups);
    result.finalRemanenceVectors.resize(nGroups);
    for (int i = 0; i < nGroups; ++i) {
        result.finalDensities[i] = m_config.model->group(i).density;
        result.finalSusceptibilities[i] = m_config.model->group(i).susceptibility;
        result.finalConductivities[i] = m_config.model->group(i).conductivity;
        result.finalRemanences[i] = m_config.model->group(i).remanence_magnitude;
        result.finalRemanenceVectors[i] = m_config.model->group(i).magnetization;
    }
    return result;
}

VectorXd InversionRunner::smoothGradient(const VectorXd& g) const {
    double w = m_config.gradientSmoothingWeight;
    if (w <= 0.0) return g;

    VectorXd smoothed = g;
    const auto* model = m_config.model.get();
    int nM = model->groupMeshCount();

    for (int mi = 0; mi < nM; ++mi) {
        const auto* mesh = model->groupMesh(mi);
        if (!mesh || mesh->vertexCount() == 0) continue;

        const auto& mappings = mesh->dofMappings();
        uint32_t nVerts = mesh->vertexCount();

        std::vector<double> vertGrad(nVerts * 3, 0.0);
        std::vector<int> vertHasAxis(nVerts * 3, 0);

        for (uint32_t li = 0; li < mappings.size(); ++li) {
            const auto& dm = mappings[li];
            int gi = model->globalDofIndex(mi, static_cast<int>(li));
            if (gi < 0 || gi >= g.size()) continue;
            uint32_t idx = dm.vertexIndex * 3 + dm.axis;
            vertGrad[idx] = g[gi];
            vertHasAxis[idx] = 1;
        }

        std::vector<double> smoothGrad(nVerts * 3, 0.0);
        for (uint32_t vi = 0; vi < nVerts; ++vi) {
            const auto& neighbors = mesh->neighborVertices(vi);
            for (int a = 0; a < 3; ++a) {
                uint32_t idx = vi * 3 + a;
                if (!vertHasAxis[idx]) continue;

                double neighborSum = 0.0;
                int nNbr = 0;
                for (uint32_t ni : neighbors) {
                    uint32_t nidx = ni * 3 + a;
                    if (vertHasAxis[nidx]) {
                        neighborSum += vertGrad[nidx];
                        ++nNbr;
                    }
                }

                if (nNbr > 0) {
                    double neighborMean = neighborSum / nNbr;
                    smoothGrad[idx] = (1.0 - w) * vertGrad[idx] + w * neighborMean;
                } else {
                    smoothGrad[idx] = vertGrad[idx];
                }
            }
        }

        for (uint32_t li = 0; li < mappings.size(); ++li) {
            const auto& dm = mappings[li];
            int gi = model->globalDofIndex(mi, static_cast<int>(li));
            if (gi < 0 || gi >= g.size()) continue;
            uint32_t idx = dm.vertexIndex * 3 + dm.axis;
            smoothed[gi] = smoothGrad[idx];
        }
    }

    return smoothed;
}

void InversionRunner::runPropertyInversion(InversionResult& /*result*/) {
    const int nGroups = m_config.model->groupCount();
    const int nPadding = 0;  // padding removed from gravity/magnetic; EM retains stub
    const size_t nDataG = m_config.observedData.size();

    // =====================================================================
    // Part A: Density optimization (gravity)
    // =====================================================================
    {
        const int nParams = nGroups;

        // Build unit response matrix for gravity (density)
        MatrixXd U(static_cast<Index>(nDataG), static_cast<Index>(nParams));
        for (int g = 0; g < nGroups; ++g) {
            U.col(g) = m_forward->computeGroupUnitResponse(g);
        }

        VectorXd propParams(nParams);
        for (int g = 0; g < nGroups; ++g) {
            propParams[g] = m_config.model->group(g).density;
        }

        // Capture reference for Tikhonov damping
        const VectorXd propRef = propParams;
        const double propLambda = m_config.propertyDamping;

        // Normalize by data RMS so the objective is ~O(N) and the
        // property_damping value is scale-invariant.
        double gravDataVar = 0.0;
        for (size_t i = 0; i < nDataG; ++i)
            gravDataVar += m_config.observedData[i].g_obs * m_config.observedData[i].g_obs;
        gravDataVar /= static_cast<double>(nDataG);
        const double gravDataRMS = std::sqrt(std::max(gravDataVar, 1.0));

        VectorXd propLower(nParams), propUpper(nParams);
        for (int g = 0; g < nGroups; ++g) {
            propLower[g] = m_config.propertyDensityMin;
            propUpper[g] = m_config.propertyDensityMax;
        }

        auto propObj = [&](const VectorXd& d) -> double {
            VectorXd gPred = U * d;
            double misfit = 0.0;
            for (size_t i = 0; i < nDataG; ++i) {
                double w = (m_config.observedData[i].g_std > 0)
                           ? (1.0 / m_config.observedData[i].g_std) : 1.0;
                double res = w * (m_config.observedData[i].g_obs - gPred[static_cast<Index>(i)]);
                misfit += 0.5 * res * res;
            }
            misfit /= gravDataRMS * gravDataRMS;
            if (propLambda > 0.0) {
                VectorXd diff = d - propRef;
                misfit += 0.5 * propLambda * diff.squaredNorm();
            }
            return misfit;
        };

        auto propGrad = [&](const VectorXd& d) -> VectorXd {
            VectorXd gPred = U * d;
            VectorXd Wr(static_cast<Index>(nDataG));
            double invScale2 = 1.0 / (gravDataRMS * gravDataRMS);
            for (size_t i = 0; i < nDataG; ++i) {
                double w = (m_config.observedData[i].g_std > 0)
                           ? (1.0 / m_config.observedData[i].g_std) : 1.0;
                Wr[static_cast<Index>(i)] = invScale2 * w * w
                    * (m_config.observedData[i].g_obs - gPred[static_cast<Index>(i)]);
            }
            VectorXd grad = -U.transpose() * Wr;
            if (propLambda > 0.0) {
                grad += propLambda * (d - propRef);
            }
            return grad;
        };

        LBFGSBOptimizer propOptimizer;
        propOptimizer.setMaxIterations(m_config.propertyInversionMaxIter);
        propOptimizer.setHistorySize(std::min(5, nParams));
        propOptimizer.setTolerance(m_config.tolerance);

        auto propResult = propOptimizer.minimize(
            propObj, propGrad, propParams, propLower, propUpper);

        for (int g = 0; g < nGroups; ++g) {
            m_config.model->group(g).density = propResult.params[g];
        }

        std::cerr << "[PROPINV density done: ";
        for (int g = 0; g < nGroups; ++g)
            std::cerr << m_config.model->group(g).density << " ";
        std::cerr << "]" << std::endl;
    }

    // =====================================================================
    // Part B: Magnetic property optimization (if magnetic data present)
    // =====================================================================
    if (m_magneticForward && !m_config.magneticData.empty()) {
        auto mode = m_config.remanenceMode;
        const size_t nDataM = m_config.magneticData.size();

        // Determine parameters-per-group from remanence mode
        int paramsPerGroup = 1; // χ only
        if (mode == RemanentMagnetizationMode::FixedVectorPerGroup)
            paramsPerGroup = 2; // χ + M_rem
        else if (mode == RemanentMagnetizationMode::VectorPerGroup)
            paramsPerGroup = 4; // χ + Mx + My + Mz

        const int nParams = paramsPerGroup * nGroups;

        MatrixXd U(static_cast<Index>(nDataM), static_cast<Index>(nParams));

        for (int g = 0; g < nGroups; ++g) {
            int col = paramsPerGroup * g;
            U.col(col) = m_magneticForward->computeGroupUnitResponse(g);

            if (mode == RemanentMagnetizationMode::FixedVectorPerGroup) {
                U.col(col + 1) = m_magneticForward->computeGroupUnitRemanenceResponse(g);
            } else if (mode == RemanentMagnetizationMode::VectorPerGroup) {
                U.col(col + 1) = m_magneticForward->computeGroupUnitResponseDirection(
                    g, Vector3d(1.0, 0.0, 0.0));
                U.col(col + 2) = m_magneticForward->computeGroupUnitResponseDirection(
                    g, Vector3d(0.0, 1.0, 0.0));
                U.col(col + 3) = m_magneticForward->computeGroupUnitResponseDirection(
                    g, Vector3d(0.0, 0.0, 1.0));
            }
        }

        VectorXd propParams(nParams);
        VectorXd propLower(nParams), propUpper(nParams);

        for (int g = 0; g < nGroups; ++g) {
            int col = paramsPerGroup * g;
            propParams[col] = m_config.model->group(g).susceptibility;
            propLower[col] = m_config.propertySusceptibilityMin;
            propUpper[col] = m_config.propertySusceptibilityMax;

            if (mode == RemanentMagnetizationMode::FixedVectorPerGroup) {
                propParams[col + 1] = m_config.model->group(g).remanence_magnitude;
                propLower[col + 1] = m_config.propertyRemanenceMin;
                propUpper[col + 1] = m_config.propertyRemanenceMax;
            } else if (mode == RemanentMagnetizationMode::VectorPerGroup) {
                propParams[col + 1] = m_config.model->group(g).magnetization.x();
                propParams[col + 2] = m_config.model->group(g).magnetization.y();
                propParams[col + 3] = m_config.model->group(g).magnetization.z();
                propLower[col + 1] = m_config.propertyRemanenceComponentMin;
                propUpper[col + 1] = m_config.propertyRemanenceComponentMax;
                propLower[col + 2] = m_config.propertyRemanenceComponentMin;
                propUpper[col + 2] = m_config.propertyRemanenceComponentMax;
                propLower[col + 3] = m_config.propertyRemanenceComponentMin;
                propUpper[col + 3] = m_config.propertyRemanenceComponentMax;
            }
        }

        // Capture reference for Tikhonov damping
        const VectorXd propMagRef = propParams;
        const double propMagLambda = m_config.propertyDamping;

        // Normalize by data RMS so the objective is ~O(N) and the
        // property_damping value is scale-invariant.  Without this,
        // magnetic data in nT produces misfit ~10^11 while the damping
        // term for χ ~0.01 with λ=0.01 contributes ~10^-8 — invisible.
        double magDataVar = 0.0;
        for (size_t i = 0; i < nDataM; ++i)
            magDataVar += m_config.magneticData[i].t_obs * m_config.magneticData[i].t_obs;
        magDataVar /= static_cast<double>(nDataM);
        const double magDataRMS = std::sqrt(std::max(magDataVar, 1.0));

        auto propObj = [&](const VectorXd& d) -> double {
            VectorXd tPred = U * d;
            double misfit = 0.0;
            for (size_t i = 0; i < nDataM; ++i) {
                double w = (m_config.magneticData[i].t_std > 0)
                           ? (1.0 / m_config.magneticData[i].t_std) : 1.0;
                double res = w * (m_config.magneticData[i].t_obs - tPred[static_cast<Index>(i)]);
                misfit += 0.5 * res * res;
            }
            misfit /= magDataRMS * magDataRMS;
            if (propMagLambda > 0.0) {
                VectorXd diff = d - propMagRef;
                misfit += 0.5 * propMagLambda * diff.squaredNorm();
            }
            return misfit;
        };

        auto propGrad = [&](const VectorXd& d) -> VectorXd {
            VectorXd tPred = U * d;
            VectorXd Wr(static_cast<Index>(nDataM));
            double invScale2 = 1.0 / (magDataRMS * magDataRMS);
            for (size_t i = 0; i < nDataM; ++i) {
                double w = (m_config.magneticData[i].t_std > 0)
                           ? (1.0 / m_config.magneticData[i].t_std) : 1.0;
                Wr[static_cast<Index>(i)] = invScale2 * w * w
                    * (m_config.magneticData[i].t_obs - tPred[static_cast<Index>(i)]);
            }
            VectorXd grad = -U.transpose() * Wr;
            if (propMagLambda > 0.0) {
                grad += propMagLambda * (d - propMagRef);
            }
            return grad;
        };

        LBFGSBOptimizer propOptimizer;
        propOptimizer.setMaxIterations(m_config.propertyInversionMaxIter);
        propOptimizer.setHistorySize(std::min(5, nParams));
        propOptimizer.setTolerance(m_config.tolerance);

        auto propResult = propOptimizer.minimize(
            propObj, propGrad, propParams, propLower, propUpper);

        // Write optimized values back to model
        for (int g = 0; g < nGroups; ++g) {
            int col = paramsPerGroup * g;
            m_config.model->group(g).susceptibility = propResult.params[col];

            if (mode == RemanentMagnetizationMode::FixedVectorPerGroup) {
                m_config.model->group(g).remanence_magnitude = propResult.params[col + 1];
            } else if (mode == RemanentMagnetizationMode::VectorPerGroup) {
                m_config.model->group(g).magnetization = Vector3d(
                    propResult.params[col + 1],
                    propResult.params[col + 2],
                    propResult.params[col + 3]);
            }
        }

        std::cerr << "[PROPINV magnetic done: ";
        for (int g = 0; g < nGroups; ++g) {
            std::cerr << "chi" << g << "=" << m_config.model->group(g).susceptibility << " ";
            if (mode == RemanentMagnetizationMode::FixedVectorPerGroup) {
                std::cerr << "M" << g << "=" << m_config.model->group(g).remanence_magnitude << " ";
            }
        }
        std::cerr << "]" << std::endl;
    }

    // =====================================================================
    // Part C: Conductivity optimization (active EM — if EM data present)
    // =====================================================================
    if (m_activeEMForward && !m_config.activeEMData.empty()) {
        const int nParams = nGroups + nPadding;
        const size_t nDataEM = m_config.activeEMData.size();

        // Build unit response matrix for active EM (conductivity)
        MatrixXd U(static_cast<Index>(nDataEM), static_cast<Index>(nParams));
        for (int g = 0; g < nGroups; ++g) {
            U.col(g) = m_activeEMForward->computeGroupUnitResponse(g);
        }
        if (m_config.enablePaddingGroup) {
            U.col(nGroups) = m_activeEMForward->computePaddingUnitResponse();
        }

        VectorXd propParams(nParams);
        for (int g = 0; g < nGroups; ++g) {
            propParams[g] = m_config.model->group(g).conductivity;
        }
        if (m_config.enablePaddingGroup) {
            propParams[nGroups] = m_activeEMForward->paddingConductivity();
        }

        VectorXd propLower(nParams), propUpper(nParams);
        for (int g = 0; g < nGroups; ++g) {
            propLower[g] = m_config.propertyConductivityMin;
            propUpper[g] = m_config.propertyConductivityMax;
        }
        if (m_config.enablePaddingGroup) {
            propLower[nGroups] = m_config.paddingConductivityLower;
            propUpper[nGroups] = m_config.paddingConductivityUpper;
        }

        auto propObj = [&](const VectorXd& d) -> double {
            VectorXd emPred = U * d;
            double misfit = 0.0;
            for (size_t i = 0; i < nDataEM; ++i) {
                double w = (m_config.activeEMData[i].d_std > 0)
                           ? (1.0 / m_config.activeEMData[i].d_std) : 1.0;
                double res = w * (m_config.activeEMData[i].d_obs
                                  - emPred[static_cast<Index>(i)]);
                misfit += 0.5 * res * res;
            }
            return misfit;
        };

        auto propGrad = [&](const VectorXd& d) -> VectorXd {
            VectorXd emPred = U * d;
            VectorXd Wr(static_cast<Index>(nDataEM));
            for (size_t i = 0; i < nDataEM; ++i) {
                double w = (m_config.activeEMData[i].d_std > 0)
                           ? (1.0 / m_config.activeEMData[i].d_std) : 1.0;
                Wr[static_cast<Index>(i)] = w * w
                    * (m_config.activeEMData[i].d_obs - emPred[static_cast<Index>(i)]);
            }
            return -U.transpose() * Wr;
        };

        LBFGSBOptimizer propOptimizer;
        propOptimizer.setMaxIterations(m_config.propertyInversionMaxIter);
        propOptimizer.setHistorySize(std::min(5, nParams));
        propOptimizer.setTolerance(m_config.tolerance);

        auto propResult = propOptimizer.minimize(
            propObj, propGrad, propParams, propLower, propUpper);

        for (int g = 0; g < nGroups; ++g) {
            m_config.model->group(g).conductivity = propResult.params[g];
        }
        if (m_config.enablePaddingGroup) {
            m_activeEMForward->setPaddingConductivity(propResult.params[nGroups]);
        }
    }

    // =====================================================================
    // Part D: Conductivity optimization (MT — if MT data present)
    // =====================================================================
    if (m_mtForward && !m_config.mtData.empty()) {
        const int nParams = nGroups + nPadding;
        const size_t nDataMT = m_config.mtData.size();

        MatrixXd U(static_cast<Index>(nDataMT), static_cast<Index>(nParams));
        for (int g = 0; g < nGroups; ++g) {
            U.col(g) = m_mtForward->computeGroupUnitResponse(g);
        }
        if (m_config.enablePaddingGroup) {
            U.col(nGroups) = m_mtForward->computePaddingUnitResponse();
        }

        VectorXd propParams(nParams);
        for (int g = 0; g < nGroups; ++g) {
            propParams[g] = m_config.model->group(g).conductivity;
        }
        if (m_config.enablePaddingGroup) {
            propParams[nGroups] = m_mtForward->paddingConductivity();
        }

        VectorXd propLower(nParams), propUpper(nParams);
        for (int g = 0; g < nGroups; ++g) {
            propLower[g] = m_config.propertyConductivityMin;
            propUpper[g] = m_config.propertyConductivityMax;
        }
        if (m_config.enablePaddingGroup) {
            propLower[nGroups] = m_config.paddingConductivityLower;
            propUpper[nGroups] = m_config.paddingConductivityUpper;
        }

        auto propObj = [&](const VectorXd& d) -> double {
            VectorXd mtPred = U * d;
            double misfit = 0.0;
            for (size_t i = 0; i < nDataMT; ++i) {
                double w = (m_config.mtData[i].z_std > 0)
                           ? (1.0 / m_config.mtData[i].z_std) : 1.0;
                double res = w * (m_config.mtData[i].zReal_obs
                                  - mtPred[static_cast<Index>(i)]);
                misfit += 0.5 * res * res;
            }
            return misfit;
        };

        auto propGrad = [&](const VectorXd& d) -> VectorXd {
            VectorXd mtPred = U * d;
            VectorXd Wr(static_cast<Index>(nDataMT));
            for (size_t i = 0; i < nDataMT; ++i) {
                double w = (m_config.mtData[i].z_std > 0)
                           ? (1.0 / m_config.mtData[i].z_std) : 1.0;
                Wr[static_cast<Index>(i)] = w * w
                    * (m_config.mtData[i].zReal_obs - mtPred[static_cast<Index>(i)]);
            }
            return -U.transpose() * Wr;
        };

        LBFGSBOptimizer propOptimizer;
        propOptimizer.setMaxIterations(m_config.propertyInversionMaxIter);
        propOptimizer.setHistorySize(std::min(5, nParams));
        propOptimizer.setTolerance(m_config.tolerance);

        auto propResult = propOptimizer.minimize(
            propObj, propGrad, propParams, propLower, propUpper);

        for (int g = 0; g < nGroups; ++g) {
            m_config.model->group(g).conductivity = propResult.params[g];
        }
        if (m_config.enablePaddingGroup) {
            m_mtForward->setPaddingConductivity(propResult.params[nGroups]);
        }
    }
}

} // namespace litho_invert
