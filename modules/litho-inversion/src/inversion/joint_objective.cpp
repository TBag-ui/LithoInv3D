#include <litho_invert/inversion/joint_objective.h>
#include <cmath>
#include <iostream>

namespace litho_invert {

JointObjective::JointObjective(std::shared_ptr<ObjectiveFunction> gravityObj,
                                std::shared_ptr<MagneticForward> magForward,
                                const MagneticData& magData,
                                double alpha_mag)
    : m_gravityObj(std::move(gravityObj))
    , m_magForward(std::move(magForward))
    , m_gravityData(m_gravityObj->data())
    , m_magData(magData)
    , m_alpha_mag(alpha_mag)
{
}

void JointObjective::addActiveEM(std::shared_ptr<EMActiveForward> activeEMForward,
                                  const ActiveEMData& data,
                                  double alpha)
{
    m_activeEMForward = std::move(activeEMForward);
    m_activeEMData = data;
    m_alpha_activeEM = alpha;
}

void JointObjective::addMT(std::shared_ptr<EMMTForward> mtForward,
                            const MTData& data,
                            double alpha)
{
    m_mtForward = std::move(mtForward);
    m_mtData = data;
    m_alpha_mt = alpha;
}

double JointObjective::evaluate(const VectorXd& params) {
    auto c = evaluateComponents(params);
    return c.total;
}

VectorXd JointObjective::gradient(const VectorXd& params) {
    // Gravity data misfit + regularization + constraint gradients (single Jacobian)
    VectorXd totalGrad = m_gravityObj->gradient(params);

    // Magnetic gradient (skip if weight is zero — Jacobian is very expensive)
    if (m_magForward && !m_magData.empty() && m_alpha_mag > 0.0) {
        MatrixXd Jmag = m_magForward->computeJacobian(params);
        VectorXd predMag = m_magForward->compute(params);
        Index nMag = static_cast<Index>(m_magData.size());
        VectorXd magWr(nMag);
        for (size_t i = 0; i < m_magData.size(); ++i) {
            double w = (m_magData[i].t_std > 0) ? (1.0 / m_magData[i].t_std) : 1.0;
            double r = m_magData[i].t_obs - predMag[static_cast<Index>(i)];
            magWr[static_cast<Index>(i)] = w * w * r;
        }
        totalGrad += -m_alpha_mag * Jmag.transpose() * magWr;
    }

    // Active EM gradient (finite-difference Jacobian)
    if (m_activeEMForward && !m_activeEMData.empty()) {
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [grad] activeEM Jacobian start, nd=" << m_activeEMData.size()
                  << " np=" << params.size() << std::endl;
#endif
        MatrixXd Jem = m_activeEMForward->computeJacobian(params);
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [grad] activeEM Jacobian done, J=" << Jem.rows() << "x" << Jem.cols() << std::endl;
#endif
        VectorXd predEM = m_activeEMForward->compute(params);
        Index nEM = static_cast<Index>(m_activeEMData.size());
        VectorXd emWr(nEM);
        for (size_t i = 0; i < m_activeEMData.size(); ++i) {
            double w = (m_activeEMData[i].d_std > 0)
                       ? (1.0 / m_activeEMData[i].d_std) : 1.0;
            double r = m_activeEMData[i].d_obs - predEM[static_cast<Index>(i)];
            emWr[static_cast<Index>(i)] = w * w * r;
        }
        totalGrad += -m_alpha_activeEM * Jem.transpose() * emWr;
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [grad] activeEM gradient done" << std::endl;
#endif
    }

    // MT gradient (finite-difference Jacobian)
    if (m_mtForward && !m_mtData.empty()) {
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [grad] MT Jacobian start, nd=" << m_mtData.size()
                  << " np=" << params.size() << std::endl;
#endif
        MatrixXd Jmt = m_mtForward->computeJacobian(params);
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [grad] MT Jacobian done, J=" << Jmt.rows() << "x" << Jmt.cols() << std::endl;
#endif
        VectorXd predMT = m_mtForward->compute(params);
        Index nMT = static_cast<Index>(m_mtData.size());
        VectorXd mtWr(nMT);
        for (size_t i = 0; i < m_mtData.size(); ++i) {
            double w = (m_mtData[i].z_std > 0)
                       ? (1.0 / m_mtData[i].z_std) : 1.0;
            double r = m_mtData[i].zReal_obs - predMT[static_cast<Index>(i)];
            mtWr[static_cast<Index>(i)] = w * w * r;
        }
        totalGrad += -m_alpha_mt * Jmt.transpose() * mtWr;
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [grad] MT gradient done" << std::endl;
#endif
    }

    return totalGrad;
}

JointObjective::JointComponents JointObjective::evaluateComponents(const VectorXd& params) {
    JointComponents c;

    // Gravity component + regularization + constraints (single evaluation)
    auto gravComp = m_gravityObj->evaluateComponents(params);
    c.gravityMisfit = gravComp.dataMisfit;
    c.regularization = gravComp.regularization;
    c.constraintPenalty = gravComp.constraintPenalty;

    // Magnetic component (skip if weight is zero — expensive compute)
    if (m_magForward && !m_magData.empty() && m_alpha_mag > 0.0) {
        VectorXd predMag = m_magForward->compute(params);
        double magMisfit = 0.0;
        for (size_t i = 0; i < m_magData.size(); ++i) {
            double w = (m_magData[i].t_std > 0) ? (1.0 / m_magData[i].t_std) : 1.0;
            double res = w * (m_magData[i].t_obs - predMag[static_cast<Index>(i)]);
            magMisfit += 0.5 * res * res;
        }
        c.magneticMisfit = magMisfit;
    }

    // Active EM component
    if (m_activeEMForward && !m_activeEMData.empty()) {
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [eval] activeEM compute start, dataCount=" << m_activeEMData.size() << std::endl;
#endif
        VectorXd predEM = m_activeEMForward->compute(params);
        double emMisfit = 0.0;
        for (size_t i = 0; i < m_activeEMData.size(); ++i) {
            double w = (m_activeEMData[i].d_std > 0)
                       ? (1.0 / m_activeEMData[i].d_std) : 1.0;
            double res = w * (m_activeEMData[i].d_obs
                              - predEM[static_cast<Index>(i)]);
            emMisfit += 0.5 * res * res;
        }
        c.activeEMMisfit = emMisfit;
    }

    // MT component
    if (m_mtForward && !m_mtData.empty()) {
#ifdef LITHO_INVERT_DEBUG
        std::cerr << "[JOINT] " "  [eval] MT compute start, dataCount=" << m_mtData.size() << std::endl;
#endif
        VectorXd predMT = m_mtForward->compute(params);
        double mtMisfit = 0.0;
        for (size_t i = 0; i < m_mtData.size(); ++i) {
            double w = (m_mtData[i].z_std > 0)
                       ? (1.0 / m_mtData[i].z_std) : 1.0;
            double r = w * (m_mtData[i].zReal_obs - predMT[static_cast<Index>(i)]);
            mtMisfit += 0.5 * r * r;
        }
        c.mtMisfit = mtMisfit;
    }

    c.total = c.gravityMisfit
              + m_alpha_mag * c.magneticMisfit
              + m_alpha_activeEM * c.activeEMMisfit
              + m_alpha_mt * c.mtMisfit
              + c.regularization + c.constraintPenalty;
    return c;
}

double JointObjective::gravityMisfit(const VectorXd& params) {
    return m_gravityObj->dataMisfit(params);
}

double JointObjective::magneticMisfit(const VectorXd& params) {
    if (!m_magForward || m_magData.empty()) return 0.0;
    VectorXd predMag = m_magForward->compute(params);
    double misfit = 0.0;
    for (size_t i = 0; i < m_magData.size(); ++i) {
        double w = (m_magData[i].t_std > 0) ? (1.0 / m_magData[i].t_std) : 1.0;
        double res = w * (m_magData[i].t_obs - predMag[static_cast<Index>(i)]);
        misfit += 0.5 * res * res;
    }
    return misfit;
}

double JointObjective::activeEMMisfit(const VectorXd& params) {
    if (!m_activeEMForward || m_activeEMData.empty()) return 0.0;
    VectorXd predEM = m_activeEMForward->compute(params);
    double misfit = 0.0;
    for (size_t i = 0; i < m_activeEMData.size(); ++i) {
        double w = (m_activeEMData[i].d_std > 0)
                   ? (1.0 / m_activeEMData[i].d_std) : 1.0;
        double res = w * (m_activeEMData[i].d_obs
                          - predEM[static_cast<Index>(i)]);
        misfit += 0.5 * res * res;
    }
    return misfit;
}

double JointObjective::mtMisfit(const VectorXd& params) {
    if (!m_mtForward || m_mtData.empty()) return 0.0;
    VectorXd predMT = m_mtForward->compute(params);
    double misfit = 0.0;
    for (size_t i = 0; i < m_mtData.size(); ++i) {
        double w = (m_mtData[i].z_std > 0)
                   ? (1.0 / m_mtData[i].z_std) : 1.0;
        double r = w * (m_mtData[i].zReal_obs - predMT[static_cast<Index>(i)]);
        misfit += 0.5 * r * r;
    }
    return misfit;
}

VectorXd JointObjective::gravityResiduals(const VectorXd& params) {
    return m_gravityObj->residuals(params);
}

VectorXd JointObjective::magneticResiduals(const VectorXd& params) {
    if (!m_magForward || m_magData.empty()) return VectorXd();
    Index n = static_cast<Index>(m_magData.size());
    VectorXd predicted = m_magForward->compute(params);
    VectorXd resid(n);
    for (Index i = 0; i < n; ++i) {
        resid[i] = m_magData[static_cast<size_t>(i)].t_obs - predicted[i];
    }
    return resid;
}

VectorXd JointObjective::activeEMResiduals(const VectorXd& params) {
    if (!m_activeEMForward || m_activeEMData.empty()) return VectorXd();
    Index n = static_cast<Index>(m_activeEMData.size());
    VectorXd predicted = m_activeEMForward->compute(params);
    VectorXd resid(n);
    for (Index i = 0; i < n; ++i) {
        resid[i] = m_activeEMData[static_cast<size_t>(i)].d_obs - predicted[i];
    }
    return resid;
}

VectorXd JointObjective::mtResiduals(const VectorXd& params) {
    if (!m_mtForward || m_mtData.empty()) return VectorXd();
    Index n = static_cast<Index>(m_mtData.size());
    VectorXd predicted = m_mtForward->compute(params);
    VectorXd resid(n);
    for (Index i = 0; i < n; ++i) {
        resid[i] = m_mtData[static_cast<size_t>(i)].zReal_obs - predicted[i];
    }
    return resid;
}

size_t JointObjective::parameterCount() const {
    return m_gravityObj->parameterCount();
}

size_t JointObjective::gravityDataCount() const {
    return m_gravityObj->dataCount();
}

size_t JointObjective::magneticDataCount() const {
    return m_magData.size();
}

size_t JointObjective::activeEMDataCount() const {
    return m_activeEMData.size();
}

size_t JointObjective::mtDataCount() const {
    return m_mtData.size();
}

// =========================================================================
// Data-weight calibration: balance gravity and magnetic gradient magnitudes.
//
// Primary: power iteration on JᵀWd²J (eigenvalue-based, à la SimPEG).
// Fallback: gradient-norm ratio if power iteration degenerates.
// =========================================================================
void JointObjective::calibrateDataWeights(const VectorXd& params) {
    if (!m_magForward || m_magData.empty() || m_alpha_mag <= 0.0) return;

    constexpr int kPowerIters = 50;
    constexpr double kPowerTol = 1e-8;

    auto& gravFwd = m_gravityObj->forwardModel();
    MatrixXd Jgrav = gravFwd.computeJacobian(params);
    MatrixXd Jmag = m_magForward->computeJacobian(params);
    Index nParam = Jgrav.cols();
    Index nGrav = Jgrav.rows();
    Index nMag = Jmag.rows();

    std::cout << "[Data scaling] J_grav: " << nGrav << "x" << nParam
              << "  |J|_max=" << std::scientific << Jgrav.cwiseAbs().maxCoeff()
              << "  |J|_mean=" << Jgrav.cwiseAbs().mean() << std::defaultfloat << std::endl;
    std::cout << "[Data scaling] J_mag:  " << nMag << "x" << nParam
              << "  |J|_max=" << std::scientific << Jmag.cwiseAbs().maxCoeff()
              << "  |J|_mean=" << Jmag.cwiseAbs().mean() << std::defaultfloat << std::endl;

    // --- Helper: power iteration for largest eigenvalue of JᵀW²J ---
    auto powerIteration = [&](const MatrixXd& J, Index nObs,
                              auto weightFn) -> double {
        VectorXd v = VectorXd::Random(nParam);
        v.normalize();
        double lambdaOld = 0.0;
        for (int iter = 0; iter < kPowerIters; ++iter) {
            VectorXd w = J * v;
            for (Index i = 0; i < nObs; ++i)
                w[i] *= weightFn(i);
            v = J.transpose() * w;
            double lambda = v.norm();
            if (lambda < 1e-60) return 0.0;
            v /= lambda;
            if (std::abs(lambda - lambdaOld) / lambda < kPowerTol) break;
            lambdaOld = lambda;
        }
        VectorXd w2 = J * v;
        for (Index i = 0; i < nObs; ++i)
            w2[i] *= weightFn(i);
        VectorXd Hv = J.transpose() * w2;
        return std::max(v.dot(Hv), 0.0);
    };

    double eigGrav = powerIteration(Jgrav, nGrav,
        [&](Index i) {
            double sig = m_gravityData[static_cast<size_t>(i)].g_std;
            return (sig > 0) ? (1.0 / (sig * sig)) : 1.0;
        });

    double eigMag = powerIteration(Jmag, nMag,
        [&](Index i) {
            double sig = m_magData[static_cast<size_t>(i)].t_std;
            return (sig > 0) ? (1.0 / (sig * sig)) : 1.0;
        });

    std::cout << "[Data scaling] λmax_grav=" << std::scientific << eigGrav
              << "  λmax_mag=" << eigMag << std::defaultfloat << std::endl;

    double oldAlpha = m_alpha_mag;
    double ratio = (eigMag > 1e-60) ? (eigGrav / eigMag) : 0.0;

    if (std::isfinite(ratio) && ratio > 0.0) {
        // Apply eigenvalue ratio but clamp the adjustment to prevent one
        // data type from being zeroed out by extreme unit-scale differences
        // (e.g. gravity in mGal vs magnetics in nT produce ratios ~1e-10).
        constexpr double kMinAdjust = 1e-6;
        constexpr double kMaxAdjust = 1e6;
        double clipped = std::max(kMinAdjust, std::min(ratio, kMaxAdjust));
        m_alpha_mag *= clipped;
        std::cout << "[Data scaling] Eigenvalue ratio=" << ratio
                  << "  clipped=" << clipped
                  << "  α_mag: " << oldAlpha << " → " << m_alpha_mag << std::endl;
        return;
    }

    // --- Fallback: gradient-norm scaling ---
    VectorXd predGrav = gravFwd.compute(params);
    VectorXd predMag = m_magForward->compute(params);

    VectorXd wResGrav(nGrav);
    for (Index i = 0; i < nGrav; ++i) {
        double sig = m_gravityData[static_cast<size_t>(i)].g_std;
        double w = (sig > 0) ? (1.0 / sig) : 1.0;
        wResGrav[i] = w * w * (m_gravityData[static_cast<size_t>(i)].g_obs - predGrav[i]);
    }
    VectorXd wResMag(nMag);
    for (Index i = 0; i < nMag; ++i) {
        double sig = m_magData[static_cast<size_t>(i)].t_std;
        double w = (sig > 0) ? (1.0 / sig) : 1.0;
        wResMag[i] = w * w * (m_magData[static_cast<size_t>(i)].t_obs - predMag[i]);
    }

    double normGrav = (Jgrav.transpose() * wResGrav).norm();
    double normMag  = (Jmag.transpose()  * wResMag).norm();
    std::cout << "[Data scaling] Eigenvalue degenerate — falling back to gradient-norm"
              << std::endl;
    std::cout << "[Data scaling] ||∇grav||=" << std::scientific << normGrav
              << "  ||∇mag||=" << normMag << std::defaultfloat << std::endl;

    if (normGrav < 1e-60) {
        std::cout << "[Data scaling] Gravity gradient near-zero, keeping α_mag="
                  << m_alpha_mag << std::endl;
    } else if (normMag > 1e-60) {
        double gradRatio = normGrav / normMag;
        m_alpha_mag *= gradRatio;
        std::cout << "[Data scaling] Gradient ratio=" << gradRatio
                  << "  α_mag: " << oldAlpha << " → " << m_alpha_mag << std::endl;
    } else {
        std::cout << "[Data scaling] Both gradients near-zero, keeping α_mag="
                  << m_alpha_mag << std::endl;
    }
}

} // namespace litho_invert

