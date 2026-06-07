#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/inversion/objective.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/em/em_data.h>
#include <litho_invert/em/em_active_forward.h>
#include <litho_invert/em/em_mt_forward.h>
#include <memory>

namespace litho_invert {

// =========================================================================
// JointObjective — weighted sum of multiple geophysical data misfits.
//
// Pattern: start with a gravity ObjectiveFunction (always present as the
// reference, weight=1), then add magnetic, active EM, and/or MT terms with
// user-set alpha weights:
//
//   total = misfit_grav + α_mag * misfit_mag + α_aem * misfit_activeEM
//         + α_mt * misfit_mt + regularization + constraints
//
// Any combination works — gravity-only, gravity+mag, gravity+EM, all four, etc.
// Empty data vectors simply skip their term.
//
// EXTENSION POINT: To add a new data type (e.g. seismic traveltime):
//   1. Add a misfit field to JointComponents
//   2. Add an add<Type>() method following the addActiveEM() pattern
//   3. Update evaluateComponents() and gradient() to include the new term
//   4. Add a residuals() accessor for DW computation
// =========================================================================
class JointObjective {
public:
    // Construct with gravity (always required) and optionally magnetic.
    // Pass magForward=nullptr and empty magData for gravity-only.
    JointObjective(std::shared_ptr<ObjectiveFunction> gravityObj,
                   std::shared_ptr<MagneticForward> magForward,
                   const MagneticData& magData,
                   double alpha_mag = 1.0);

    // --- Add optional EM objectives ---
    // Call these after construction to add active EM and/or MT terms.
    void addActiveEM(std::shared_ptr<EMActiveForward> activeEMForward,
                     const ActiveEMData& data,
                     double alpha = 1.0);

    void addMT(std::shared_ptr<EMMTForward> mtForward,
               const MTData& data,
               double alpha = 1.0);

    // --- Evaluation ---
    double evaluate(const VectorXd& params);
    VectorXd gradient(const VectorXd& params);

    // Component-wise breakdown for iteration logging.
    // Fields that are zero indicate the corresponding data type is absent.
    struct JointComponents {
        double gravityMisfit = 0.0;
        double magneticMisfit = 0.0;
        double activeEMMisfit = 0.0;
        double mtMisfit = 0.0;
        double regularization = 0.0;
        double constraintPenalty = 0.0;
        double total = 0.0;
    };
    JointComponents evaluateComponents(const VectorXd& params);

    // --- Individual misfits ---
    double gravityMisfit(const VectorXd& params);
    double magneticMisfit(const VectorXd& params);
    double activeEMMisfit(const VectorXd& params);
    double mtMisfit(const VectorXd& params);

    // --- Residual vectors (for DW spatial autocorrelation computation) ---
    VectorXd gravityResiduals(const VectorXd& params);
    VectorXd magneticResiduals(const VectorXd& params);
    VectorXd activeEMResiduals(const VectorXd& params);
    VectorXd mtResiduals(const VectorXd& params);

    // --- Parameter / data counts ---
    size_t parameterCount() const;
    size_t gravityDataCount() const;
    size_t magneticDataCount() const;
    size_t activeEMDataCount() const;
    size_t mtDataCount() const;

    // --- Weight setters (can be tuned during inversion callbacks) ---
    void setMagneticWeight(double alpha) { m_alpha_mag = alpha; }
    double magneticWeight() const { return m_alpha_mag; }
    void setActiveEMWeight(double alpha) { m_alpha_activeEM = alpha; }
    double activeEMWeight() const { return m_alpha_activeEM; }
    void setMTWeight(double alpha) { m_alpha_mt = alpha; }
    double mtWeight() const { return m_alpha_mt; }

    // --- Eigenvalue-based Hessian scaling (SimPEG ScalingMultipleDataMisfits_ByEig) ---
    // Computes eig_max(JᵀWdJ) for gravity and magnetics at the current params,
    // then scales magnetic misfit by eig_max_grav / eig_max_mag so both data
    // types contribute comparable gradient magnitudes at the starting model.
    // Call once at iteration 0 before the first geometry optimisation step.
    void calibrateDataWeights(const VectorXd& params);

    // --- Accessors ---
    std::shared_ptr<ObjectiveFunction> gravityObjective() const { return m_gravityObj; }
    std::shared_ptr<MagneticForward> magneticForward() const { return m_magForward; }
    std::shared_ptr<EMActiveForward> activeEMForward() const { return m_activeEMForward; }
    std::shared_ptr<EMMTForward> mtForward() const { return m_mtForward; }

    const GravityData& gravityData() const { return m_gravityData; }
    const MagneticData& magneticData() const { return m_magData; }
    const ActiveEMData& activeEMData() const { return m_activeEMData; }
    const MTData& mtData() const { return m_mtData; }

private:
    std::shared_ptr<ObjectiveFunction> m_gravityObj;
    std::shared_ptr<MagneticForward> m_magForward;
    std::shared_ptr<EMActiveForward> m_activeEMForward;
    std::shared_ptr<EMMTForward> m_mtForward;

    // Gravity data reference is valid for the lifetime of m_gravityObj
    const GravityData& m_gravityData;

    // Stored by value: magnetic and EM data are copied in so references
    // remain valid regardless of how the caller manages the originals
    MagneticData m_magData;
    ActiveEMData m_activeEMData;
    MTData m_mtData;

    double m_alpha_mag = 1.0;       // magnetic misfit weight
    double m_alpha_activeEM = 1.0;  // active EM misfit weight
    double m_alpha_mt = 1.0;        // MT misfit weight
};

} // namespace litho_invert
