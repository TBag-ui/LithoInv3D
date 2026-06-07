#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/core/geometry.h>
#include <cmath>
#include <cassert>

namespace litho_invert {

namespace {
    constexpr double kEps = 1e-15;

    // Contribution of a single triangular face to the vertical gravity component.
    // Returns the dimensionless geometric factor n_z * (d_sum - h * omega).
    double faceNagyContribution(const Vector3d& obs,
                                 const Vector3d& a,
                                 const Vector3d& b,
                                 const Vector3d& c) {
        Vector3d n = (b - a).cross(c - a);
        double area_x2 = n.norm();
        if (area_x2 < kEps) return 0.0;
        Vector3d unit_n = n / area_x2;

        double nz = unit_n.z();
        double omega = solidAngle(a, b, c, obs);
        double h = (a - obs).dot(unit_n);

        double d_sum = 0.0;
        const Vector3d* verts[3] = {&a, &b, &c};
        for (int i = 0; i < 3; ++i) {
            const Vector3d& p1 = *verts[i];
            const Vector3d& p2 = *verts[(i + 1) % 3];
            Vector3d e = p2 - p1;
            double L = e.norm();
            if (L < kEps) continue;
            Vector3d e_hat = e / L;
            Vector3d perp = unit_n.cross(e_hat);
            double d_perp = (obs - p1).dot(perp);
            double Lambda = lineIntegralTerm(p1, p2, obs);
            d_sum += d_perp * Lambda;
        }

        return nz * (d_sum - h * omega);
    }

    // Compute faceNagyContribution AND its gradients w.r.t. the three vertices.
    // Uses local forward-difference FD on just this face (10 face evaluations).
    // h=1e-4 (~0.1 mm in model coordinates) gives stable gradients for
    // Nagy polyhedron faces at typical inversion scales (100-5000 m).
    void faceNagyWithGradient(const Vector3d& obs,
                              const Vector3d& a, const Vector3d& b, const Vector3d& c,
                              double& value,
                              Vector3d& grad_a, Vector3d& grad_b, Vector3d& grad_c) {
        constexpr double h = 1e-4;
        value = faceNagyContribution(obs, a, b, c);

        // Vertex a
        {
            Vector3d ap = a;
            ap.x() += h;
            grad_a.x() = (faceNagyContribution(obs, ap, b, c) - value) / h;
            ap = a; ap.y() += h;
            grad_a.y() = (faceNagyContribution(obs, ap, b, c) - value) / h;
            ap = a; ap.z() += h;
            grad_a.z() = (faceNagyContribution(obs, ap, b, c) - value) / h;
        }
        // Vertex b
        {
            Vector3d bp = b;
            bp.x() += h;
            grad_b.x() = (faceNagyContribution(obs, a, bp, c) - value) / h;
            bp = b; bp.y() += h;
            grad_b.y() = (faceNagyContribution(obs, a, bp, c) - value) / h;
            bp = b; bp.z() += h;
            grad_b.z() = (faceNagyContribution(obs, a, bp, c) - value) / h;
        }
        // Vertex c
        {
            Vector3d cp = c;
            cp.x() += h;
            grad_c.x() = (faceNagyContribution(obs, a, b, cp) - value) / h;
            cp = c; cp.y() += h;
            grad_c.y() = (faceNagyContribution(obs, a, b, cp) - value) / h;
            cp = c; cp.z() += h;
            grad_c.z() = (faceNagyContribution(obs, a, b, cp) - value) / h;
        }
    }
}

// =========================================================================
// GravityForward implementation
// =========================================================================

GravityForward::GravityForward(std::shared_ptr<LithologyModel> model,
                               const GravityData& data)
    : m_model(std::move(model)), m_data(data) {
}

size_t GravityForward::dataCount() const {
    return m_data.size();
}

size_t GravityForward::parameterCount() const {
    return static_cast<size_t>(m_model->totalDofCount());
}

VectorXd GravityForward::compute(const VectorXd& params) {
    m_model->applyParameterVector(params);

    VectorXd result(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        const Vector3d& obs = m_data[i].position;
        double gz = 0.0;
        for (int g = 0; g < m_model->groupCount(); ++g) {
            gz += computeUnitGravity(g, obs);
        }
        result(static_cast<Index>(i)) = gz;
    }
    return result;
}

MatrixXd GravityForward::computeJacobian(const VectorXd& params) {
    return computeGravityAnalyticJacobian(params);
}

double GravityForward::computeUnitGravity(int groupIndex, const Vector3d& obs) const {
    assert(groupIndex >= 0 && groupIndex < m_model->groupCount());
    const SurfaceMesh* mesh = m_model->groupMesh(groupIndex);
    if (!mesh || mesh->triangleCount() == 0) return 0.0;
    double rho = m_model->group(groupIndex).density;
    return gravityClosedMesh(obs, *mesh, rho);
}

VectorXd GravityForward::computeGroupUnitResponse(int groupIndex) const {
    assert(groupIndex >= 0 && groupIndex < m_model->groupCount());
    const SurfaceMesh* mesh = m_model->groupMesh(groupIndex);
    if (!mesh || mesh->triangleCount() == 0) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    VectorXd result(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        result[static_cast<Index>(i)] = gravityClosedMesh(
            m_data[i].position, *mesh, 1.0);
    }
    return result;
}

// -----------------------------------------------------------------------
// gravityFacet
// -----------------------------------------------------------------------

double GravityForward::gravityFacet(const Vector3d& obs,
                                     const Vector3d& a,
                                     const Vector3d& b,
                                     const Vector3d& c) {
    return faceNagyContribution(obs, a, b, c);
}

// -----------------------------------------------------------------------
// gravityClosedMesh — iterate all triangles of a closed mesh
// -----------------------------------------------------------------------

double GravityForward::gravityClosedMesh(const Vector3d& obs,
                                          const SurfaceMesh& closedMesh,
                                          double density_g_per_cm3) {
    if (closedMesh.triangleCount() == 0) return 0.0;

    double total = 0.0;
    for (uint32_t t = 0; t < closedMesh.triangleCount(); ++t) {
        const Triangle& tri = closedMesh.triangle(t);
        const Vector3d& a = closedMesh.vertex(tri.v0).position;
        const Vector3d& b = closedMesh.vertex(tri.v1).position;
        const Vector3d& c = closedMesh.vertex(tri.v2).position;
        total += faceNagyContribution(obs, a, b, c);
    }

    double scale = G_SI * density_g_per_cm3 * DENSITY_SCALE * M_S2_TO_MGAL;
    return total * scale;
}

// =========================================================================
// Analytic Jacobian
// =========================================================================

MatrixXd GravityForward::computeGravityAnalyticJacobian(const VectorXd& params) const {
    const size_t nd = dataCount();
    const size_t np = parameterCount();
    const int nGroups = m_model->groupCount();

    // Apply current parameters so mesh vertices are at the right positions
    m_model->applyParameterVector(params);

    MatrixXd J(static_cast<Index>(nd), static_cast<Index>(np));
    J.setZero();

    // Pre-build per-mesh vertex-active mask: which vertices have any DOF
    struct MeshInfo {
        std::vector<uint8_t> vertexHasAnyDof; // bool per vertex
        std::vector<std::vector<int>> vertexDofAxes; // which axes have DOFs (0=x,1=y,2=z)
        const SurfaceMesh* mesh;
        int meshIndex;
    };
    std::vector<MeshInfo> meshInfos(nGroups);
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* mesh = m_model->groupMesh(g);
        if (!mesh) continue;
        MeshInfo& info = meshInfos[g];
        info.mesh = mesh;
        info.meshIndex = g;
        uint32_t nv = mesh->vertexCount();
        info.vertexHasAnyDof.assign(nv, 0);
        info.vertexDofAxes.resize(nv);
        for (const auto& dm : mesh->dofMappings()) {
            info.vertexHasAnyDof[dm.vertexIndex] = 1;
            info.vertexDofAxes[dm.vertexIndex].push_back(static_cast<int>(dm.axis));
        }
    }

    // Pre-build localDof → globalDof lookup (avoid linear search in inner loop)
    std::vector<std::vector<int>> localToGlobal(nGroups);
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* mesh = m_model->groupMesh(g);
        if (!mesh) continue;
        size_t nLocal = mesh->dofMappings().size();
        localToGlobal[g].resize(nLocal, -1);
        for (size_t ld = 0; ld < nLocal; ++ld) {
            localToGlobal[g][ld] = m_model->globalDofIndex(g, static_cast<int>(ld));
        }
    }

    const double scale = G_SI * DENSITY_SCALE * M_S2_TO_MGAL;

    // Per-observation loop
    for (Index iObs = 0; iObs < static_cast<Index>(nd); ++iObs) {
        const Vector3d& obs = m_data[static_cast<size_t>(iObs)].position;

        // Per-mesh vertex gradient accumulators (allocated once per thread)
        std::vector<std::vector<Vector3d>> vertexGrads(nGroups);
        for (int g = 0; g < nGroups; ++g) {
            const MeshInfo& info = meshInfos[g];
            if (!info.mesh || info.mesh->triangleCount() == 0) continue;
            vertexGrads[g].assign(info.mesh->vertexCount(), Vector3d::Zero());
        }

        // Accumulate face gradients into per-vertex gradient vectors
        for (int g = 0; g < nGroups; ++g) {
            const MeshInfo& info = meshInfos[g];
            if (!info.mesh || info.mesh->triangleCount() == 0) continue;

            double rho = m_model->group(g).density;
            double faceScale = scale * rho;

            for (uint32_t t = 0; t < info.mesh->triangleCount(); ++t) {
                const Triangle& tri = info.mesh->triangle(t);
                const Vector3d& a = info.mesh->vertex(tri.v0).position;
                const Vector3d& b = info.mesh->vertex(tri.v1).position;
                const Vector3d& c = info.mesh->vertex(tri.v2).position;

                double fVal;
                Vector3d ga, gb, gc;
                faceNagyWithGradient(obs, a, b, c, fVal, ga, gb, gc);

                auto& vg = vertexGrads[g];
                vg[tri.v0] += ga * faceScale;
                vg[tri.v1] += gb * faceScale;
                vg[tri.v2] += gc * faceScale;
            }
        }

        // Map vertex gradients → local DOF gradients → global Jacobian row
        for (int g = 0; g < nGroups; ++g) {
            const MeshInfo& info = meshInfos[g];
            if (!info.mesh || info.mesh->dofMappings().empty()) continue;

            const auto& vg = vertexGrads[g];
            const auto& mappings = info.mesh->dofMappings();

            for (uint32_t localDof = 0; localDof < static_cast<uint32_t>(mappings.size()); ++localDof) {
                const auto& dm = mappings[localDof];
                double physGrad = vg[dm.vertexIndex][dm.axis];
                if (std::abs(physGrad) < 1e-30) continue;

                int globalDof = localToGlobal[g][localDof];
                if (globalDof < 0) continue;

                double axisScale = info.mesh->axisScale(dm.axis);
                J(iObs, globalDof) += physGrad * axisScale;
            }
        }
    }

    return J;
}

// =========================================================================
// ForwardModel base class default implementations
// =========================================================================

MatrixXd ForwardModel::computeJacobian(const VectorXd& params) {
    const size_t nd = dataCount();
    const size_t np = parameterCount();

    MatrixXd J(static_cast<Index>(nd), static_cast<Index>(np));
    const double h = m_fdStep;

    for (size_t j = 0; j < np; ++j) {
        VectorXd params_plus = params;
        VectorXd params_minus = params;
        params_plus(static_cast<Index>(j)) += h;
        params_minus(static_cast<Index>(j)) -= h;

        VectorXd f_plus = compute(params_plus);
        VectorXd f_minus = compute(params_minus);

        for (size_t i = 0; i < nd; ++i) {
            J(static_cast<Index>(i), static_cast<Index>(j)) =
                (f_plus(static_cast<Index>(i)) - f_minus(static_cast<Index>(i))) / (2.0 * h);
        }
    }

    return J;
}

void ForwardModel::computeBoth(const VectorXd& params,
                                VectorXd& predicted,
                                MatrixXd& jacobian) {
    predicted = compute(params);
    jacobian = computeJacobian(params);
}

} // namespace litho_invert

