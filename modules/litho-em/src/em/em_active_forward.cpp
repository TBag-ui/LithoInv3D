#include <litho_invert/em/em_active_forward.h>
#include <litho_invert/litho/lithology_model.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

namespace litho_invert {

EMActiveForward::EMActiveForward(std::shared_ptr<LithologyModel> model,
                                 const std::vector<EMSource>& sources,
                                 const std::vector<EMReceiver>& receivers,
                                 const ActiveEMData& data,
                                 const EMConfig& config)
    : EMForwardModel(model, config)
    , m_sources(sources)
    , m_receivers(receivers)
    , m_data(data)
{
    if (!m_solver) {
        m_solver = createEMSolver(config.solverMethod);
    }
}

size_t EMActiveForward::dataCount() const {
    return m_data.size();
}

size_t EMActiveForward::parameterCount() const {
    if (!m_model) return 0;
    return static_cast<size_t>(m_model->totalDofCount());
}

VectorXd EMActiveForward::compute(const VectorXd& params) {
    if (!m_model) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    // Apply the parameter vector to the model (updates surface vertex z-values)
    m_model->applyParameterVector(params);

    // Group data rows by (sourceIndex, receiverIndex) to minimize
    // solver calls.  Each unique pair is one solve.
    //
    // Map: (sourceIndex, receiverIndex) -> vector of data indices
    std::map<std::pair<int, int>, std::vector<size_t>> srGroups;
    for (size_t i = 0; i < m_data.size(); ++i) {
        srGroups[{m_data[i].sourceIndex, m_data[i].receiverIndex}].push_back(i);
    }

    VectorXd result(static_cast<Index>(m_data.size()));

    for (const auto& kv : srGroups) {
        int srcIdx = kv.first.first;
        int recvIdx = kv.first.second;
        const std::vector<size_t>& dataIndices = kv.second;

        // Build the active set for this source (trust-region subsetting)
        std::vector<double> activeVolume;
        std::unordered_set<int> activeGroups = buildActiveSet(srcIdx, activeVolume);

        // Collect unique gate indices for this source-receiver pair
        std::set<int> gateSet;
        for (size_t di : dataIndices) {
            gateSet.insert(m_data[di].gateIndex);
        }
        std::vector<int> gateIndices(gateSet.begin(), gateSet.end());

        // Compute response for this source-receiver pair
        VectorXd srResponse = computeSourceReceiver(
            srcIdx, recvIdx, gateIndices, &activeGroups);

        // Map the response back to individual data rows
        for (size_t di : dataIndices) {
            // Find which gate index in our sorted list corresponds to this datum
            int gIdx = m_data[di].gateIndex;
            auto it = std::find(gateIndices.begin(), gateIndices.end(), gIdx);
            if (it != gateIndices.end()) {
                int localIdx = static_cast<int>(std::distance(gateIndices.begin(), it));
                result[static_cast<Index>(di)] = srResponse[localIdx];
            }
        }
    }

    return result;
}

VectorXd EMActiveForward::computeSourceReceiver(
    int sourceIndex, int receiverIndex,
    const std::vector<int>& gateIndices,
    const std::unordered_set<int>* activeGroups) const
{
    // Determine which time gates / frequencies to compute
    int nResponses = static_cast<int>(gateIndices.size());
    VectorXd result(nResponses);

    // Build receiver position list (single receiver for this call)
    std::vector<Vector3d> recvPositions;
    recvPositions.push_back(m_receivers[receiverIndex].position);

    // Set up config with the correct gates/frequencies
    EMConfig gateConfig = m_config;
    gateConfig.timeGates_s.clear();
    for (int gi : gateIndices) {
        if (!m_config.timeGates_s.empty() && gi < static_cast<int>(m_config.timeGates_s.size())) {
            gateConfig.timeGates_s.push_back(m_config.timeGates_s[gi]);
        }
    }

    // If a solver is available, use it
    if (m_solver) {
        VectorXd solverResult = m_solver->solveActive(
            m_sources[sourceIndex].position,
            recvPositions,
            *m_model,
            gateConfig);

        // solverResult layout: all gates for recv 0, then all gates for recv 1, etc.
        if (solverResult.size() == nResponses) {
            return solverResult;
        }
    }

    // --- Fallback: approximate dipole response ---
    // When no solver is available (e.g. "ie" not yet implemented),
    // compute a simple conductive-sphere-in-halfspace approximation
    // for each gate.  This is dimensionally correct but not accurate
    // for real inversion — it exists so the framework compiles, links,
    // and runs end-to-end while the full IE solver is developed.
    //
    // For each gate/frequency, compute a single equivalent skin depth
    // and approximate the response as proportional to the conductivity-
    // weighted volume intersected by the decaying eddy currents.
    for (int i = 0; i < nResponses; ++i) {
        double response = 0.0;  // nT/s or pT depending on domain

        double sigma = 1e-4;  // default background
        if (activeGroups && !activeGroups->empty()) {
            // Use the maximum conductivity among active groups
            for (int gi : *activeGroups) {
                if (gi >= 0 && gi < m_model->groupCount()) {
                    sigma = std::max(sigma, m_model->group(gi).conductivity);
                }
            }
        }

        double skinD = EMForwardModel::skinDepthTimeDomain(sigma, gateConfig.timeGates_s[i]);

        // Simple geometric factor: receiver-source separation
        double R = (m_receivers[receiverIndex].position - m_sources[sourceIndex].position).norm();
        if (R < 1.0) R = 1.0;  // avoid singularity

        // Dimensionally: response ~ sigma * volume / R^3 (simplified dipole)
        // This is a placeholder; real IE solver replaces this.
        response = sigma * skinD * skinD * skinD / (R * R * R);

        result[i] = response;
    }

    return result;
}

std::unordered_set<int> EMActiveForward::buildActiveSet(
    int sourceIndex, std::vector<double>& activeVolume) const
{
    std::unordered_set<int> activeGroups;

    if (!m_model || m_model->groupCount() == 0) {
        activeVolume = {0, 0, 0, 0, 0, 0};
        return activeGroups;
    }

    // 1. Compute skin depth from the dominant conductivity in the model
    double sigma = dominantConductivity();
    double skinD = 1e9;  // default: infinite

    if (m_config.autoComputeSkinDepth) {
        if (!m_config.timeGates_s.empty()) {
            // Use the earliest gate for the shallowest (smallest) skin depth —
            // this is the most restrictive, so it's the safe choice for subsetting
            double tMin = *std::min_element(m_config.timeGates_s.begin(),
                                            m_config.timeGates_s.end());
            skinD = EMForwardModel::skinDepthTimeDomain(sigma, tMin);
        } else if (m_config.frequency_Hz > 0.0) {
            skinD = EMForwardModel::skinDepth(sigma, m_config.frequency_Hz);
        }
    }

    // 2. Determine margin multiplier based on iteration phase
    double marginMultiplier;
    double condThreshold = 0.0;
    if (m_currentIteration < m_config.subsettingIterationsWide) {
        // Wide phase: larger margin, all groups
        marginMultiplier = m_config.subsettingSkinDepthMultiplierWide;
        condThreshold = 0.0;  // include all groups
    } else {
        // Narrow phase: smaller margin, only conductive groups
        marginMultiplier = m_config.subsettingSkinDepthMultiplierNarrow;
        condThreshold = m_config.subsettingConductivityThreshold;
    }

    // 3. Build bounding box around source-receiver line
    const Vector3d& srcPos = m_sources[sourceIndex].position;
    double margin = marginMultiplier * skinD;

    double xmin = srcPos.x(), xmax = srcPos.x();
    double ymin = srcPos.y(), ymax = srcPos.y();
    double zmin = srcPos.z(), zmax = srcPos.z();

    // Expand to include all receivers for this source
    for (size_t r = 0; r < m_receivers.size(); ++r) {
        const Vector3d& rp = m_receivers[r].position;
        xmin = std::min(xmin, rp.x());
        xmax = std::max(xmax, rp.x());
        ymin = std::min(ymin, rp.y());
        ymax = std::max(ymax, rp.y());
        zmin = std::min(zmin, rp.z());
        zmax = std::max(zmax, rp.z());
    }

    // Apply margin
    xmin -= margin;  xmax += margin;
    ymin -= margin;  ymax += margin;
    zmin -= margin;  zmax += margin;  // note: positive-up z, so zmin < zmax

    activeVolume = {xmin, xmax, ymin, ymax, zmin, zmax};

    // 4. Filter litho groups
    for (int g = 0; g < m_model->groupCount(); ++g) {
        // Conductivity filter (narrow phase only)
        if (condThreshold > 0.0 && m_model->group(g).conductivity < condThreshold) {
            continue;
        }

        // Check if any vertex of the group's bounding surfaces falls within
        // the active volume.  Use the group's top and bottom surfaces.
        bool inVolume = false;

        // Check vertices of the group's bottom surface against the bounding box.
        const SurfaceMesh* botSurf = nullptr;

        if (g < m_model->groupMeshCount()) {
            botSurf = m_model->groupMesh(g);
        }
        if (g == 0) {
            // Group 0: top is the flat top (z=0 surface), bottom is surface[0]
            // Check surface[0] vertices only
            if (botSurf) {
                for (uint32_t v = 0; v < botSurf->vertexCount(); ++v) {
                    const Vector3d& p = botSurf->vertex(v).position;
                    if (p.x() >= xmin && p.x() <= xmax &&
                        p.y() >= ymin && p.y() <= ymax &&
                        p.z() >= zmin && p.z() <= zmax) {
                        inVolume = true;
                        break;
                    }
                }
            }
        } else if (g > 0 && g <= m_model->groupMeshCount()) {
            // Group g: top = surface[g-1], bottom = surface[g]
            if (botSurf) {
                for (uint32_t v = 0; v < botSurf->vertexCount(); ++v) {
                    const Vector3d& p = botSurf->vertex(v).position;
                    if (p.x() >= xmin && p.x() <= xmax &&
                        p.y() >= ymin && p.y() <= ymax &&
                        p.z() >= zmin && p.z() <= zmax) {
                        inVolume = true;
                        break;
                    }
                }
            }
        }

        if (inVolume) {
            activeGroups.insert(g);
        }
    }

    // Always include at least one group
    if (activeGroups.empty()) {
        for (int g = 0; g < m_model->groupCount(); ++g) {
            activeGroups.insert(g);
        }
    }

    return activeGroups;
}

double EMActiveForward::dominantConductivity() const {
    if (!m_model || m_model->groupCount() == 0) return 1e-4;

    // Use the maximum conductivity among all groups as the "dominant" value.
    // This ensures the skin depth used for subsetting is conservative
    // (shallowest penetration = largest margin needed).
    double maxSigma = 0.0;
    for (int g = 0; g < m_model->groupCount(); ++g) {
        maxSigma = std::max(maxSigma, m_model->group(g).conductivity);
    }
    return (maxSigma > 0.0) ? maxSigma : 1e-4;
}

VectorXd EMActiveForward::computeGroupUnitResponse(int groupIndex) {
    if (!m_model || groupIndex < 0 || groupIndex >= m_model->groupCount()) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    // Store original conductivities and set all to zero
    std::vector<double> originalConductivities(m_model->groupCount());
    for (int g = 0; g < m_model->groupCount(); ++g) {
        originalConductivities[g] = m_model->group(g).conductivity;
        m_model->group(g).conductivity = 0.0;
    }
    // Set only the target group to unit conductivity
    m_model->group(groupIndex).conductivity = 1.0;

    // Compute response with unit conductivity in the target group
    VectorXd unitResponse = compute(m_model->assembleParameterVector());

    // Restore original conductivities
    for (int g = 0; g < m_model->groupCount(); ++g) {
        m_model->group(g).conductivity = originalConductivities[g];
    }

    return unitResponse;
}

VectorXd EMActiveForward::computePaddingUnitResponse() const {
    // Placeholder: padding contribution for active EM.
    // The deep half-space background response — typically small if the
    // padding is deep enough that eddy currents don't reach it.
    return VectorXd::Zero(static_cast<Index>(m_data.size()));
}

void EMActiveForward::enablePadding(bool enabled, double paddingDepth) {
    m_paddingEnabled = enabled;
    m_paddingDepth = paddingDepth;
    if (enabled) {
        buildFlatDeepBottom();
    }
}

void EMActiveForward::buildFlatDeepBottom() {
    // Same pattern as GravityForward::buildFlatDeepBottom():
    // create a deep flat surface at m_paddingDepth using the topology
    // of the deepest model surface.
    if (!m_model || m_model->groupMeshCount() == 0) return;

    const SurfaceMesh* deepestSurf = m_model->groupMesh(m_model->groupMeshCount() - 1);
    m_flatDeepBottom = std::make_shared<SurfaceMesh>();
    m_flatDeepBottom->setName("em_padding_bottom");
    for (uint32_t v = 0; v < deepestSurf->vertexCount(); ++v) {
        Vector3d pos = deepestSurf->vertex(v).position;
        pos.z() = m_paddingDepth;
        m_flatDeepBottom->addVertex(pos);
    }
    for (uint32_t t = 0; t < deepestSurf->triangleCount(); ++t) {
        const Triangle& tri = deepestSurf->triangle(t);
        m_flatDeepBottom->addTriangle(tri.v0, tri.v1, tri.v2);
    }
    m_flatDeepBottom->buildNeighbors();
}

} // namespace litho_invert
