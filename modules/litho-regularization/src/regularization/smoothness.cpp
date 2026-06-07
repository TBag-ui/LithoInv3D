#include <litho_invert/regularization/smoothness.h>
#include <litho_invert/surface/surface_mesh.h>
#include <cmath>
#include <unordered_map>

namespace litho_invert {

SurfaceSmoothness::SurfaceSmoothness(std::shared_ptr<LithologyModel> model)
    : m_model(std::move(model))
{
    // Build neighbor relationships for all group meshes
    for (int s = 0; s < m_model->groupMeshCount(); ++s) {
        SurfaceMesh* surf = m_model->groupMesh(s);
        if (surf->vertexCount() > 0) {
            surf->buildNeighbors();
        }
    }
}

double SurfaceSmoothness::evaluate(const VectorXd& params) {
    m_model->applyParameterVector(params);

    double penalty = 0.0;

    for (int s = 0; s < m_model->groupMeshCount(); ++s) {
        const SurfaceMesh* surf = m_model->groupMesh(s);
        const uint32_t nVerts = surf->vertexCount();
        if (nVerts == 0) continue;

        for (uint32_t vi = 0; vi < nVerts; ++vi) {
            const Vertex& vert = surf->vertex(vi);
            if (vert.freedom == VertexFreedom::FIXED) continue;

            const auto& nbrs = surf->neighborVertices(vi);
            if (nbrs.empty()) continue;

            // Only consider neighbors with the same freedom type.
            // This filters out side-wall connections in closed meshes,
            // where a contact-surface vertex is adjacent to a closure-
            // surface vertex at a very different depth.
            Vector3d mean = Vector3d::Zero();
            int nSame = 0;
            for (uint32_t nbr : nbrs) {
                if (surf->vertex(nbr).freedom == vert.freedom) {
                    mean += surf->vertex(nbr).position;
                    ++nSame;
                }
            }
            if (nSame == 0) continue;
            mean /= static_cast<double>(nSame);

            const Vector3d& pos = vert.position;

            switch (vert.freedom) {
                case VertexFreedom::FIXED:
                    break;
                case VertexFreedom::Z_ONLY: {
                    double dz = pos.z() - mean.z();
                    penalty += 0.5 * dz * dz;
                    break;
                }
                case VertexFreedom::X_ONLY: {
                    double dx = pos.x() - mean.x();
                    penalty += 0.5 * dx * dx;
                    break;
                }
                case VertexFreedom::Y_ONLY: {
                    double dy = pos.y() - mean.y();
                    penalty += 0.5 * dy * dy;
                    break;
                }
                case VertexFreedom::ALONG_VECTOR: {
                    Vector3d mv = vert.moveVector.normalized();
                    double diff = (pos - mean).dot(mv);
                    penalty += 0.5 * diff * diff;
                    break;
                }
                case VertexFreedom::XY_FREE: {
                    double dx = pos.x() - mean.x();
                    double dy = pos.y() - mean.y();
                    penalty += 0.5 * (dx * dx + dy * dy);
                    break;
                }
                case VertexFreedom::XZ_FREE: {
                    double dx = pos.x() - mean.x();
                    double dz = pos.z() - mean.z();
                    penalty += 0.5 * (dx * dx + dz * dz);
                    break;
                }
                case VertexFreedom::YZ_FREE: {
                    double dy = pos.y() - mean.y();
                    double dz = pos.z() - mean.z();
                    penalty += 0.5 * (dy * dy + dz * dz);
                    break;
                }
                case VertexFreedom::XYZ_FREE: {
                    penalty += 0.5 * (pos - mean).squaredNorm();
                    break;
                }
            }
        }
    }

    return m_weight * penalty;
}

VectorXd SurfaceSmoothness::gradient(const VectorXd& params) {
    m_model->applyParameterVector(params);

    const uint32_t nDof = m_model->totalDofCount();
    VectorXd grad = VectorXd::Zero(static_cast<Index>(nDof));

    for (int s = 0; s < m_model->groupMeshCount(); ++s) {
        const SurfaceMesh* surf = m_model->groupMesh(s);
        const auto& dofMappings = surf->dofMappings();
        const uint32_t nVerts = surf->vertexCount();

        if (nVerts == 0) continue;

        // ---------------------------------------------------------------
        // Step 1: Precompute means and diffs for all vertices,
        //         using only same-freedom neighbors to exclude side walls.
        // ---------------------------------------------------------------
        std::vector<Vector3d> diffs(nVerts, Vector3d::Zero());
        std::vector<uint32_t> nNbrs(nVerts, 0);

        for (uint32_t vi = 0; vi < nVerts; ++vi) {
            const Vertex& vert = surf->vertex(vi);
            if (vert.freedom == VertexFreedom::FIXED) continue;

            const auto& nbrs = surf->neighborVertices(vi);

            Vector3d mean = Vector3d::Zero();
            uint32_t nSame = 0;
            for (uint32_t nbr : nbrs) {
                if (surf->vertex(nbr).freedom == vert.freedom) {
                    mean += surf->vertex(nbr).position;
                    ++nSame;
                }
            }
            if (nSame == 0) continue;
            mean /= static_cast<double>(nSame);
            nNbrs[vi] = nSame;

            diffs[vi] = vert.position - mean;
        }

        // When control-point stride is active, compute full-resolution
        // per-vertex gradients for all free axes and downsample through
        // the chain rule to control-point DOFs.
        if (surf->controlPointStride() > 0) {
            // Per-vertex smoothness gradient (3D vector per vertex)
            std::vector<Vector3d> fullVertexGrad(nVerts, Vector3d::Zero());

            // Own term: ∂(0.5·||v - mean||²)/∂v = v - mean = diffs[vi],
            // filtered by the vertex freedom type.
            for (uint32_t vi = 0; vi < nVerts; ++vi) {
                if (nNbrs[vi] == 0) continue;
                const Vertex& vert = surf->vertex(vi);
                switch (vert.freedom) {
                    case VertexFreedom::Z_ONLY:
                        fullVertexGrad[vi].z() = diffs[vi].z();
                        break;
                    case VertexFreedom::X_ONLY:
                        fullVertexGrad[vi].x() = diffs[vi].x();
                        break;
                    case VertexFreedom::Y_ONLY:
                        fullVertexGrad[vi].y() = diffs[vi].y();
                        break;
                    case VertexFreedom::XY_FREE:
                        fullVertexGrad[vi].x() = diffs[vi].x();
                        fullVertexGrad[vi].y() = diffs[vi].y();
                        break;
                    case VertexFreedom::XZ_FREE:
                        fullVertexGrad[vi].x() = diffs[vi].x();
                        fullVertexGrad[vi].z() = diffs[vi].z();
                        break;
                    case VertexFreedom::YZ_FREE:
                        fullVertexGrad[vi].y() = diffs[vi].y();
                        fullVertexGrad[vi].z() = diffs[vi].z();
                        break;
                    case VertexFreedom::XYZ_FREE:
                        fullVertexGrad[vi] = diffs[vi];
                        break;
                    case VertexFreedom::ALONG_VECTOR: {
                        Vector3d mv = vert.moveVector.normalized();
                        fullVertexGrad[vi] = diffs[vi].dot(mv) * mv;
                        break;
                    }
                    default: break;
                }
            }

            // Neighbor term: vertex ui contributes -diffs[ui]/nNbrs
            // to each same-freedom neighbor's gradient, per active axis.
            for (uint32_t ui = 0; ui < nVerts; ++ui) {
                if (nNbrs[ui] == 0) continue;
                const Vertex& uVert = surf->vertex(ui);
                Vector3d contrib = Vector3d::Zero();
                switch (uVert.freedom) {
                    case VertexFreedom::Z_ONLY:
                        contrib.z() = -diffs[ui].z() / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::X_ONLY:
                        contrib.x() = -diffs[ui].x() / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::Y_ONLY:
                        contrib.y() = -diffs[ui].y() / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::XY_FREE:
                        contrib.x() = -diffs[ui].x() / static_cast<double>(nNbrs[ui]);
                        contrib.y() = -diffs[ui].y() / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::XZ_FREE:
                        contrib.x() = -diffs[ui].x() / static_cast<double>(nNbrs[ui]);
                        contrib.z() = -diffs[ui].z() / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::YZ_FREE:
                        contrib.y() = -diffs[ui].y() / static_cast<double>(nNbrs[ui]);
                        contrib.z() = -diffs[ui].z() / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::XYZ_FREE:
                        contrib = -diffs[ui] / static_cast<double>(nNbrs[ui]);
                        break;
                    case VertexFreedom::ALONG_VECTOR: {
                        Vector3d mv_u = uVert.moveVector.normalized();
                        double diffProj = diffs[ui].dot(mv_u);
                        contrib = -diffProj * mv_u / static_cast<double>(nNbrs[ui]);
                        break;
                    }
                    default: continue;
                }

                const auto& nbrs = surf->neighborVertices(ui);
                for (uint32_t nbr : nbrs) {
                    if (surf->vertex(nbr).freedom != uVert.freedom) continue;
                    const Vertex& nbrVert = surf->vertex(nbr);
                    if (nbrVert.freedom == VertexFreedom::ALONG_VECTOR) {
                        Vector3d mv_nbr = nbrVert.moveVector.normalized();
                        fullVertexGrad[nbr] += contrib.dot(mv_nbr) * mv_nbr;
                    } else {
                        fullVertexGrad[nbr] += contrib;
                    }
                }
            }

            // Map per-vertex gradients to local DOFs: one entry per
            // control-point vertex per free axis.
            std::vector<double> localGrad(dofMappings.size(), 0.0);
            for (uint32_t d = 0; d < static_cast<uint32_t>(dofMappings.size()); ++d) {
                uint32_t vi = dofMappings[d].vertexIndex;
                uint8_t ax = dofMappings[d].axis;
                if (vi < nVerts && ax < 3) {
                    localGrad[d] = fullVertexGrad[vi][ax];
                }
            }

            // Route local DOF contributions to global gradient.
            // Chain rule: dϕ/d(param) = dϕ/d(pos) * d(pos)/d(param)
            // where d(pos)/d(param) = axisScale for that axis.
            for (uint32_t d = 0; d < static_cast<uint32_t>(dofMappings.size()); ++d) {
                int gIdx = m_model->globalDofIndex(s, static_cast<int>(d));
                if (gIdx >= 0) {
                    grad[gIdx] += localGrad[d] * surf->axisScale(dofMappings[d].axis);
                }
            }
        } else {
            // ---------------------------------------------------------------
            // Step 2: Build (vertex -> local DOF entries) mapping
            // ---------------------------------------------------------------
            struct DofEntry {
                Index dofIndex;  // local DOF index within this mesh
                uint8_t axis;
            };
            std::vector<std::vector<DofEntry>> vertexDofs(nVerts);

            for (uint32_t d = 0; d < static_cast<uint32_t>(dofMappings.size()); ++d) {
                uint32_t vi = dofMappings[d].vertexIndex;
                vertexDofs[vi].push_back({
                    static_cast<Index>(d),
                    dofMappings[d].axis
                });
            }

            // Accumulate into a local gradient first, then route to global
            std::vector<double> localGrad(dofMappings.size(), 0.0);

            // ---------------------------------------------------------------
            // Step 3: Own term
            // ---------------------------------------------------------------
            for (uint32_t vi = 0; vi < nVerts; ++vi) {
                if (nNbrs[vi] == 0) continue;
                const Vertex& vert = surf->vertex(vi);

                for (const auto& entry : vertexDofs[vi]) {
                    if (vert.freedom == VertexFreedom::ALONG_VECTOR) {
                        Vector3d mv = vert.moveVector.normalized();
                        localGrad[entry.dofIndex] += diffs[vi].dot(mv);
                    } else {
                        localGrad[entry.dofIndex] += diffs[vi][entry.axis];
                    }
                }
            }

            // ---------------------------------------------------------------
            // Step 4: Neighbor term (only same-freedom neighbors)
            // ---------------------------------------------------------------
            for (uint32_t ui = 0; ui < nVerts; ++ui) {
                if (nNbrs[ui] == 0) continue;
                const Vertex& uVert = surf->vertex(ui);
                const auto& nbrs = surf->neighborVertices(ui);

                Vector3d contrib = Vector3d::Zero();
                bool hasContrib = false;

                switch (uVert.freedom) {
                    case VertexFreedom::FIXED:
                        continue;
                    case VertexFreedom::Z_ONLY:
                        contrib.z() = -diffs[ui].z() / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    case VertexFreedom::X_ONLY:
                        contrib.x() = -diffs[ui].x() / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    case VertexFreedom::Y_ONLY:
                        contrib.y() = -diffs[ui].y() / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    case VertexFreedom::ALONG_VECTOR: {
                        Vector3d mv_u = uVert.moveVector.normalized();
                        double diffProj = diffs[ui].dot(mv_u);
                        contrib = -diffProj * mv_u / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    }
                    case VertexFreedom::XY_FREE:
                        contrib.x() = -diffs[ui].x() / static_cast<double>(nNbrs[ui]);
                        contrib.y() = -diffs[ui].y() / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    case VertexFreedom::XZ_FREE:
                        contrib.x() = -diffs[ui].x() / static_cast<double>(nNbrs[ui]);
                        contrib.z() = -diffs[ui].z() / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    case VertexFreedom::YZ_FREE:
                        contrib.y() = -diffs[ui].y() / static_cast<double>(nNbrs[ui]);
                        contrib.z() = -diffs[ui].z() / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                    case VertexFreedom::XYZ_FREE:
                        contrib = -diffs[ui] / static_cast<double>(nNbrs[ui]);
                        hasContrib = true;
                        break;
                }

                if (!hasContrib) continue;

                for (uint32_t nbr : nbrs) {
                    const Vertex& nbrVert = surf->vertex(nbr);
                    if (nbrVert.freedom != uVert.freedom) continue;
                    for (const auto& entry : vertexDofs[nbr]) {
                        if (nbrVert.freedom == VertexFreedom::ALONG_VECTOR) {
                            Vector3d mv_nbr = nbrVert.moveVector.normalized();
                            localGrad[entry.dofIndex] += contrib.dot(mv_nbr);
                        } else {
                            localGrad[entry.dofIndex] += contrib[entry.axis];
                        }
                    }
                }
            }

            // Route local DOF contributions to global gradient.
            // Chain rule: dϕ/d(param) = dϕ/d(pos) * d(pos)/d(param)
            // where d(pos)/d(param) = axisScale for that axis.
            for (uint32_t d = 0; d < static_cast<uint32_t>(dofMappings.size()); ++d) {
                int gIdx = m_model->globalDofIndex(s, static_cast<int>(d));
                if (gIdx >= 0) {
                    grad[gIdx] += localGrad[d] * surf->axisScale(dofMappings[d].axis);
                }
            }
        }
    }

    return m_weight * grad;
}

} // namespace litho_invert

