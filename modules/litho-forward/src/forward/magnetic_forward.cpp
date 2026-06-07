#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/core/geometry.h>
#include <cmath>
#include <cassert>

namespace litho_invert {

namespace {
    constexpr double kEps = 1e-15;
    constexpr double kPi = 3.14159265358979323846;

    // Compute magneticFacet AND its gradients w.r.t. the three vertices.
    // Uses local forward-difference FD on just this face (10 face evaluations).
    void magneticFacetWithGradient(const Vector3d& obs,
                                   const Vector3d& a, const Vector3d& b, const Vector3d& c,
                                   const Vector3d& magDirection,
                                   const Vector3d& fieldUnitVector,
                                   double& value,
                                   Vector3d& grad_a, Vector3d& grad_b, Vector3d& grad_c) {
        constexpr double h = 1e-4;

        auto eval = [&](const Vector3d& va, const Vector3d& vb, const Vector3d& vc) {
            return MagneticForward::magneticFacet(obs, va, vb, vc,
                                                   magDirection, fieldUnitVector);
        };

        value = eval(a, b, c);

        // Vertex a
        {
            Vector3d ap = a;
            ap.x() += h; grad_a.x() = (eval(ap, b, c) - value) / h;
            ap = a; ap.y() += h; grad_a.y() = (eval(ap, b, c) - value) / h;
            ap = a; ap.z() += h; grad_a.z() = (eval(ap, b, c) - value) / h;
        }
        // Vertex b
        {
            Vector3d bp = b;
            bp.x() += h; grad_b.x() = (eval(a, bp, c) - value) / h;
            bp = b; bp.y() += h; grad_b.y() = (eval(a, bp, c) - value) / h;
            bp = b; bp.z() += h; grad_b.z() = (eval(a, bp, c) - value) / h;
        }
        // Vertex c
        {
            Vector3d cp = c;
            cp.x() += h; grad_c.x() = (eval(a, b, cp) - value) / h;
            cp = c; cp.y() += h; grad_c.y() = (eval(a, b, cp) - value) / h;
            cp = c; cp.z() += h; grad_c.z() = (eval(a, b, cp) - value) / h;
        }
    }
}

// =========================================================================
// MagneticForward implementation
// =========================================================================

MagneticForward::MagneticForward(std::shared_ptr<LithologyModel> model,
                                 const MagneticData& data,
                                 double inc_deg, double dec_deg, double field_nT)
    : m_model(std::move(model)), m_data(data), m_fieldStrength_nT(field_nT) {
    double I = inc_deg * kPi / 180.0;
    double D = dec_deg * kPi / 180.0;
    double cosI = std::cos(I);
    m_fieldUnitVector = Vector3d(cosI * std::cos(D),
                                  cosI * std::sin(D),
                                  -std::sin(I));
}

size_t MagneticForward::dataCount() const {
    return m_data.size();
}

size_t MagneticForward::parameterCount() const {
    return static_cast<size_t>(m_model->totalDofCount());
}

Vector3d MagneticForward::earthFieldUnitVector() const {
    return m_fieldUnitVector;
}

VectorXd MagneticForward::compute(const VectorXd& params) {
    m_model->applyParameterVector(params);

    VectorXd result(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        const Vector3d& obs = m_data[i].position;
        double dt = 0.0;
        for (int g = 0; g < m_model->groupCount(); ++g) {
            dt += computeUnitMagnetic(g, obs);
        }
        result[static_cast<Index>(i)] = dt;
    }
    return result;
}

MatrixXd MagneticForward::computeJacobian(const VectorXd& params) {
    return computeMagneticAnalyticJacobian(params);
}

void MagneticForward::buildGroupMagnetization(int groupIndex,
                                              double& magAmp_Apm,
                                              Vector3d& magDirection) const {
    const auto& group = m_model->group(groupIndex);

    double inducedAmp = group.susceptibility * m_fieldStrength_nT / 100.0;
    Vector3d totalM = inducedAmp * m_fieldUnitVector;

    switch (m_remanenceMode) {
    case RemanentMagnetizationMode::EffectiveSusceptibility:
        break;

    case RemanentMagnetizationMode::FixedVectorPerGroup: {
        double I = group.remanence_inclination * kPi / 180.0;
        double D = group.remanence_declination * kPi / 180.0;
        double cosI = std::cos(I);
        Vector3d R_hat(cosI * std::cos(D),
                        cosI * std::sin(D),
                        -std::sin(I));
        totalM += group.remanence_magnitude * R_hat;
        break;
    }

    case RemanentMagnetizationMode::VectorPerGroup:
        totalM += group.magnetization;
        break;
    }

    magAmp_Apm = totalM.norm();
    if (magAmp_Apm < kEps) {
        magDirection = m_fieldUnitVector;
    } else {
        magDirection = totalM / magAmp_Apm;
    }
}

double MagneticForward::computeUnitMagnetic(int groupIndex, const Vector3d& obs) const {
    assert(groupIndex >= 0 && groupIndex < m_model->groupCount());
    const SurfaceMesh* mesh = m_model->groupMesh(groupIndex);
    if (!mesh || mesh->triangleCount() == 0) return 0.0;

    double magAmp;
    Vector3d magDir;
    buildGroupMagnetization(groupIndex, magAmp, magDir);

    return magneticClosedMeshVector(obs, *mesh, magAmp, magDir, m_fieldUnitVector);
}

VectorXd MagneticForward::computeGroupUnitResponse(int groupIndex) const {
    assert(groupIndex >= 0 && groupIndex < m_model->groupCount());
    const SurfaceMesh* mesh = m_model->groupMesh(groupIndex);
    if (!mesh || mesh->triangleCount() == 0) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    VectorXd result(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        result[static_cast<Index>(i)] = magneticClosedMesh(
            m_data[i].position, *mesh, 1.0,
            m_fieldUnitVector, m_fieldStrength_nT);
    }
    return result;
}

VectorXd MagneticForward::computeGroupUnitResponseDirection(int groupIndex,
                                                          const Vector3d& magDirection) const {
    assert(groupIndex >= 0 && groupIndex < m_model->groupCount());
    const SurfaceMesh* mesh = m_model->groupMesh(groupIndex);
    if (!mesh || mesh->triangleCount() == 0) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    VectorXd result(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        result[static_cast<Index>(i)] = magneticClosedMeshVector(
            m_data[i].position, *mesh, 1.0,
            magDirection, m_fieldUnitVector);
    }
    return result;
}

VectorXd MagneticForward::computeGroupUnitRemanenceResponse(int groupIndex) const {
    const auto& group = m_model->group(groupIndex);
    double I = group.remanence_inclination * kPi / 180.0;
    double D = group.remanence_declination * kPi / 180.0;
    double cosI = std::cos(I);
    Vector3d R_hat(cosI * std::cos(D),
                    cosI * std::sin(D),
                    -std::sin(I));
    return computeGroupUnitResponseDirection(groupIndex, R_hat);
}

// =========================================================================
// Analytic Jacobian
// =========================================================================

MatrixXd MagneticForward::computeMagneticAnalyticJacobian(const VectorXd& params) const {
    const size_t nd = dataCount();
    const size_t np = parameterCount();
    const int nGroups = m_model->groupCount();

    m_model->applyParameterVector(params);

    MatrixXd J(static_cast<Index>(nd), static_cast<Index>(np));
    J.setZero();

    // Pre-build per-mesh vertex-active mask
    struct MeshInfo {
        const SurfaceMesh* mesh = nullptr;
        int meshIndex = 0;
    };
    std::vector<MeshInfo> meshInfos(nGroups);
    for (int g = 0; g < nGroups; ++g) {
        meshInfos[g].mesh = m_model->groupMesh(g);
        meshInfos[g].meshIndex = g;
    }

    // Pre-build magnetization info per group
    std::vector<double> magAmps(nGroups);
    std::vector<Vector3d> magDirs(nGroups);
    for (int g = 0; g < nGroups; ++g) {
        buildGroupMagnetization(g, magAmps[g], magDirs[g]);
    }

    // Pre-build localDof → globalDof lookup (avoid linear search in inner loop)
    std::vector<std::vector<int>> localToGlobal(nGroups);
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* mesh = meshInfos[g].mesh;
        if (!mesh) continue;
        size_t nLocal = mesh->dofMappings().size();
        localToGlobal[g].resize(nLocal, -1);
        for (size_t ld = 0; ld < nLocal; ++ld) {
            localToGlobal[g][ld] = m_model->globalDofIndex(g, static_cast<int>(ld));
        }
    }

    constexpr double kInv4Pi = 1.0 / (4.0 * kPi);
    // Scale factor: magnetization(A/m) → nT: M * 100 / (4*pi), then negated
    // magneticClosedMeshVector returns: -scale * geometricFactor
    // where scale = magAmp * 100.0 / (4.0 * kPi)
    // So the physical-to-geometric mapping is: dt = -scale * geom
    // and ∂(dt)/∂v = -scale * ∂(geom)/∂v

    // Per-observation loop
    for (Index iObs = 0; iObs < static_cast<Index>(nd); ++iObs) {
        const Vector3d& obs = m_data[static_cast<size_t>(iObs)].position;

        // Per-mesh vertex gradient accumulators
        std::vector<std::vector<Vector3d>> vertexGrads(nGroups);
        for (int g = 0; g < nGroups; ++g) {
            const SurfaceMesh* mesh = meshInfos[g].mesh;
            if (!mesh || mesh->triangleCount() == 0) continue;
            vertexGrads[g].assign(mesh->vertexCount(), Vector3d::Zero());
        }

        // Accumulate face gradients
        for (int g = 0; g < nGroups; ++g) {
            const SurfaceMesh* mesh = meshInfos[g].mesh;
            if (!mesh || mesh->triangleCount() == 0) continue;
            if (magAmps[g] < kEps) continue;

            // Scale: ∂(dt)/∂(geom) = -scale
            double faceScale = -magAmps[g] * 100.0 * kInv4Pi;

            for (uint32_t t = 0; t < mesh->triangleCount(); ++t) {
                const Triangle& tri = mesh->triangle(t);
                const Vector3d& a = mesh->vertex(tri.v0).position;
                const Vector3d& b = mesh->vertex(tri.v1).position;
                const Vector3d& c = mesh->vertex(tri.v2).position;

                double fVal;
                Vector3d ga, gb, gc;
                magneticFacetWithGradient(obs, a, b, c,
                                          magDirs[g], m_fieldUnitVector,
                                          fVal, ga, gb, gc);

                auto& vg = vertexGrads[g];
                vg[tri.v0] += ga * faceScale;
                vg[tri.v1] += gb * faceScale;
                vg[tri.v2] += gc * faceScale;
            }
        }

        // Map vertex gradients → local DOF gradients → global Jacobian row
        for (int g = 0; g < nGroups; ++g) {
            const SurfaceMesh* mesh = meshInfos[g].mesh;
            if (!mesh || mesh->dofMappings().empty()) continue;

            const auto& vg = vertexGrads[g];
            const auto& mappings = mesh->dofMappings();

            for (uint32_t localDof = 0; localDof < static_cast<uint32_t>(mappings.size()); ++localDof) {
                const auto& dm = mappings[localDof];
                double physGrad = vg[dm.vertexIndex][dm.axis];
                if (std::abs(physGrad) < 1e-30) continue;

                int globalDof = localToGlobal[g][localDof];
                if (globalDof < 0) continue;

                double axisScale = mesh->axisScale(dm.axis);
                J(iObs, globalDof) += physGrad * axisScale;
            }
        }
    }

    return J;
}

// =========================================================================
// static: magneticFacet
// =========================================================================

double MagneticForward::magneticFacet(const Vector3d& obs,
                                       const Vector3d& a, const Vector3d& b,
                                       const Vector3d& c,
                                       const Vector3d& magDirection,
                                       const Vector3d& fieldUnitVector) {
    Vector3d edge1 = b - a;
    Vector3d edge2 = c - a;
    Vector3d n = edge1.cross(edge2);
    double area_x2 = n.norm();
    if (area_x2 < kEps) return 0.0;
    Vector3d n_hat = n / area_x2;

    double omega = solidAngle(a, b, c, obs);

    double mn = magDirection.dot(n_hat);
    double fn = fieldUnitVector.dot(n_hat);
    double contribution = mn * fn * omega;

    const Vector3d* verts[3] = {&a, &b, &c};
    for (int i = 0; i < 3; ++i) {
        const Vector3d& p1 = *verts[i];
        const Vector3d& p2 = *verts[(i + 1) % 3];
        Vector3d e = p2 - p1;
        double L = e.norm();
        if (L < kEps) continue;
        Vector3d e_hat = e / L;
        Vector3d t = n_hat.cross(e_hat);

        double mt = magDirection.dot(t);
        double ft = fieldUnitVector.dot(t);
        double Lambda = lineIntegralTerm(p1, p2, obs);
        contribution += mt * ft * Lambda;
    }

    return contribution;
}

// =========================================================================
// static: magneticClosedMesh (induced-only)
// =========================================================================

double MagneticForward::magneticClosedMesh(const Vector3d& obs,
                                            const SurfaceMesh& closedMesh,
                                            double susceptibility_SI,
                                            const Vector3d& fieldUnitVector,
                                            double fieldStrength_nT) {
    double magAmp = susceptibility_SI * fieldStrength_nT / 100.0;
    return magneticClosedMeshVector(obs, closedMesh, magAmp,
                                     fieldUnitVector, fieldUnitVector);
}

// =========================================================================
// static: magneticClosedMeshVector (arbitrary magnetization)
// =========================================================================

double MagneticForward::magneticClosedMeshVector(const Vector3d& obs,
                                                  const SurfaceMesh& closedMesh,
                                                  double magnetization_amp_Apm,
                                                  const Vector3d& magDirection,
                                                  const Vector3d& fieldUnitVector) {
    if (closedMesh.triangleCount() == 0) return 0.0;

    double geometricFactor = 0.0;
    for (uint32_t t = 0; t < closedMesh.triangleCount(); ++t) {
        const Triangle& tri = closedMesh.triangle(t);
        const Vector3d& a = closedMesh.vertex(tri.v0).position;
        const Vector3d& b = closedMesh.vertex(tri.v1).position;
        const Vector3d& c = closedMesh.vertex(tri.v2).position;
        geometricFactor += magneticFacet(obs, a, b, c,
                                          magDirection, fieldUnitVector);
    }

    double scale = magnetization_amp_Apm * 100.0 / (4.0 * kPi);
    return -scale * geometricFactor;
}

} // namespace litho_invert

