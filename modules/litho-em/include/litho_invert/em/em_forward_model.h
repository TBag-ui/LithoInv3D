#pragma once
#include <litho_invert/forward/forward_model.h>
#include <litho_invert/em/em_data.h>
#include <litho_invert/em/em_solver.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>
#include <vector>

namespace litho_invert {

// =========================================================================
// EMForwardModel — abstract base for all EM forward models.
//
// Extends ForwardModel (which provides compute(), computeJacobian(),
// computeBoth()) with conductivity-aware methods shared by active-source
// and passive-source EM.
//
// Subclasses:
//   EMActiveForward — airborne + large-loop TEM/FEM
//   EMMTForward     — magnetotelluric (plane-wave)
//
// EXTENSION POINT: To add a new EM method (e.g. CSEM, borehole EM):
//   1. Extend this class
//   2. Implement compute(), dataCount(), parameterCount()
//   3. Override applyParameters() if your method uses different
//      parameterization of the model
//   4. Add registration in JointObjective
// =========================================================================
class EMForwardModel : public ForwardModel {
public:
    EMForwardModel(std::shared_ptr<LithologyModel> model,
                   const EMConfig& config);

    virtual ~EMForwardModel() = default;

    // --- Subsetting interface ---
    // Set the current inversion iteration.  The forward model uses this
    // to decide between wide/narrow trust-region margins.
    virtual void setCurrentIteration(int iteration) { m_currentIteration = iteration; }
    int currentIteration() const { return m_currentIteration; }

    // --- Subsetting configuration ---
    // Read/write access so the user can tune during inversion callbacks.
    EMConfig& config() { return m_config; }
    const EMConfig& config() const { return m_config; }

    // --- Solver ---
    // Set a custom solver (default is created from config.solverMethod).
    void setSolver(std::shared_ptr<EMSolver> solver) { m_solver = std::move(solver); }
    std::shared_ptr<EMSolver> solver() const { return m_solver; }

    // --- Conductivity property inversion ---
    // Compute unit-conductivity response for a single lithology group.
    // Used during property inversion phases to build the response matrix.
    virtual VectorXd computeGroupUnitResponse(int groupIndex) = 0;

    // --- Padding conductivity ---
    virtual void setPaddingConductivity(double sigma) { m_paddingConductivity = sigma; }
    double paddingConductivity() const { return m_paddingConductivity; }

    // --- Skin depth ---
    // δ = 503 / sqrt(σ · f)  [m, S/m, Hz]
    // For time-domain: δ = 503 · sqrt(t_gate / σ)  [approximate]
    static double skinDepth(double conductivity_Sm, double frequency_Hz);
    static double skinDepthTimeDomain(double conductivity_Sm, double timeGate_s);

protected:
    std::shared_ptr<LithologyModel> m_model;
    std::shared_ptr<EMSolver> m_solver;
    EMConfig m_config;
    int m_currentIteration = 0;
    double m_paddingConductivity = 1e-4;  // S/m
};

} // namespace litho_invert

