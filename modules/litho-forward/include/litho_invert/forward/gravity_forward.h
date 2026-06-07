#pragma once
#include <litho_invert/forward/forward_model.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>

namespace litho_invert {

class GravityForward : public ForwardModel {
public:
    GravityForward(std::shared_ptr<LithologyModel> model,
                   const GravityData& data);

    VectorXd compute(const VectorXd& params) override;
    MatrixXd computeJacobian(const VectorXd& params) override;
    size_t dataCount() const override;
    size_t parameterCount() const override;

    // Compute vertical gravity of a closed triangulated mesh at one observation
    // point. All triangles must have CCW winding for outward-pointing normals.
    // Uses the Nagy (2000) face integral formula. Result in mGal.
    static double gravityClosedMesh(const Vector3d& obs,
                                     const SurfaceMesh& closedMesh,
                                     double density_g_per_cm3);

    // Compute vertical gravity contribution of a single triangle facet.
    // Returns the dimensionless geometric factor (without G*rho scaling).
    static double gravityFacet(const Vector3d& obs,
                                const Vector3d& a,
                                const Vector3d& b,
                                const Vector3d& c);

    // Compute unit-density gravity response of a single lithology group
    // at all observation points (uses current geometry, does NOT apply params).
    VectorXd computeGroupUnitResponse(int groupIndex) const;

    // Gravitational constant in SI units
    static constexpr double G_SI = 6.67430e-11;
    // Density conversion: g/cm^3 -> kg/m^3
    static constexpr double DENSITY_SCALE = 1000.0;
    // Conversion: m/s^2 -> mGal
    static constexpr double M_S2_TO_MGAL = 1e5;

private:
    std::shared_ptr<LithologyModel> m_model;
    const GravityData& m_data;

    double computeUnitGravity(int groupIndex, const Vector3d& obs) const;

    // Analytic Jacobian: face-level vertex gradients accumulated to global DOFs.
    // ~10 face evaluations per face instead of 2·N_dof full forward evals.
    MatrixXd computeGravityAnalyticJacobian(const VectorXd& params) const;
};

} // namespace litho_invert

