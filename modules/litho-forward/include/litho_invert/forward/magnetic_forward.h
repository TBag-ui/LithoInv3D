#pragma once
#include <litho_invert/forward/forward_model.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>

namespace litho_invert {

class MagneticForward : public ForwardModel {
public:
    MagneticForward(std::shared_ptr<LithologyModel> model,
                    const MagneticData& data,
                    double inc_deg, double dec_deg, double field_nT);

    VectorXd compute(const VectorXd& params) override;
    MatrixXd computeJacobian(const VectorXd& params) override;
    size_t dataCount() const override;
    size_t parameterCount() const override;

    Vector3d earthFieldUnitVector() const;

    // Compute total-field magnetic anomaly of a closed triangulated mesh.
    // All triangles must have CCW winding for outward-pointing normals.
    // Uses the Okabe/Plouff analytical polyhedron formula. Result in nT.
    static double magneticClosedMesh(const Vector3d& obs,
                                      const SurfaceMesh& closedMesh,
                                      double susceptibility_SI,
                                      const Vector3d& fieldUnitVector,
                                      double fieldStrength_nT);

    // Same as above but with arbitrary magnetization vector M (A/m).
    static double magneticClosedMeshVector(const Vector3d& obs,
                                            const SurfaceMesh& closedMesh,
                                            double magnetization_amp_Apm,
                                            const Vector3d& magDirection,
                                            const Vector3d& fieldUnitVector);

    // Magnetic field contribution of a single triangle facet.
    // Returns the dimensionless geometric factor.
    static double magneticFacet(const Vector3d& obs,
                                 const Vector3d& a, const Vector3d& b,
                                 const Vector3d& c,
                                 const Vector3d& magDirection,
                                 const Vector3d& fieldUnitVector);

    // Build effective magnetization vector for a group (A/m, direction).
    void buildGroupMagnetization(int groupIndex,
                                  double& magAmp_Apm,
                                  Vector3d& magDirection) const;

    // Unit-susceptibility response of a single group at all observation points.
    VectorXd computeGroupUnitResponse(int groupIndex) const;

    // Unit response for a specified magnetization direction.
    VectorXd computeGroupUnitResponseDirection(int groupIndex,
                                                const Vector3d& magDirection) const;

    // Unit response for the fixed remnant direction.
    VectorXd computeGroupUnitRemanenceResponse(int groupIndex) const;

    // --- Remanent magnetization mode ---
    void setRemanenceMode(RemanentMagnetizationMode mode) { m_remanenceMode = mode; }
    RemanentMagnetizationMode remanenceMode() const { return m_remanenceMode; }

private:
    std::shared_ptr<LithologyModel> m_model;
    const MagneticData& m_data;

    Vector3d m_fieldUnitVector;
    double m_fieldStrength_nT;

    RemanentMagnetizationMode m_remanenceMode = RemanentMagnetizationMode::EffectiveSusceptibility;

    double computeUnitMagnetic(int groupIndex, const Vector3d& obs) const;

    // Analytic Jacobian: face-level vertex gradients accumulated to global DOFs.
    MatrixXd computeMagneticAnalyticJacobian(const VectorXd& params) const;
};

} // namespace litho_invert

