#include <litho_invert/litho/lithology_model.h>
#include <limits>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <set>
#include <iostream>

namespace litho_invert {

// ---------------------------------------------------------------------------
// Group management
// ---------------------------------------------------------------------------

int LithologyModel::addGroup(const LithoGroup& group) {
    if (!group.isValid()) {
        return -1;
    }
    auto dup = std::find_if(m_groups.begin(), m_groups.end(),
                            [&group](const LithoGroup& g) { return g.id == group.id; });
    if (dup != m_groups.end()) {
        return -1;
    }
    m_groups.push_back(group);
    return static_cast<int>(m_groups.size()) - 1;
}

void LithologyModel::removeGroup(int groupId) {
    auto it = std::find_if(m_groups.begin(), m_groups.end(),
                           [groupId](const LithoGroup& g) { return g.id == groupId; });
    if (it == m_groups.end()) return;
    int idx = static_cast<int>(std::distance(m_groups.begin(), it));
    m_groups.erase(it);
    if (idx < static_cast<int>(m_groupMeshes.size())) {
        m_groupMeshes.erase(m_groupMeshes.begin() + idx);
    }
    m_globalDofsDirty = true;
}

const std::vector<LithoGroup>& LithologyModel::groups() const {
    return m_groups;
}

const LithoGroup& LithologyModel::group(int index) const {
    return m_groups.at(index);
}

int LithologyModel::groupCount() const {
    return static_cast<int>(m_groups.size());
}

LithoGroup& LithologyModel::group(int index) {
    return m_groups.at(index);
}

// ---------------------------------------------------------------------------
// Group mesh management
// ---------------------------------------------------------------------------

void LithologyModel::setGroupMesh(int groupIndex, std::shared_ptr<SurfaceMesh> mesh) {
    if (groupIndex < 0) return;
    if (static_cast<size_t>(groupIndex) >= m_groupMeshes.size()) {
        m_groupMeshes.resize(groupIndex + 1);
    }
    m_groupMeshes[groupIndex] = std::move(mesh);
    m_globalDofsDirty = true;
}

SurfaceMesh* LithologyModel::groupMesh(int index) {
    if (index < 0 || static_cast<size_t>(index) >= m_groupMeshes.size()) return nullptr;
    return m_groupMeshes[index].get();
}

const SurfaceMesh* LithologyModel::groupMesh(int index) const {
    if (index < 0 || static_cast<size_t>(index) >= m_groupMeshes.size()) return nullptr;
    return m_groupMeshes[index].get();
}

int LithologyModel::groupMeshCount() const {
    return static_cast<int>(m_groupMeshes.size());
}

// ---------------------------------------------------------------------------
// Label grid
// ---------------------------------------------------------------------------

void LithologyModel::setLabelGrid(const std::vector<int>& labels,
                                   int nx, int ny, int nz,
                                   double x0, double y0, double z0,
                                   double dx, double dy, double dz) {
    m_labelGrid = labels;
    m_labelNx = nx;
    m_labelNy = ny;
    m_labelNz = nz;
    m_labelX0 = x0;
    m_labelY0 = y0;
    m_labelZ0 = z0;
    m_labelDx = dx;
    m_labelDy = dy;
    m_labelDz = dz;
}

bool LithologyModel::hasLabelGrid() const {
    return !m_labelGrid.empty() && m_labelNx > 0 && m_labelNy > 0 && m_labelNz > 0;
}

// ---------------------------------------------------------------------------
// Global DOF deduplication
// ---------------------------------------------------------------------------

void LithologyModel::rebuildGlobalDofs() const {
    m_globalDofRefs.clear();

    // Collect all local DOF mappings from all group meshes
    struct Candidate {
        uint32_t meshIndex;
        uint32_t localDofIndex;
        int64_t qx, qy, qz;  // quantized position (µm)
        uint8_t axis;
    };
    std::vector<Candidate> candidates;

    for (size_t mi = 0; mi < m_groupMeshes.size(); ++mi) {
        if (!m_groupMeshes[mi]) continue;
        const auto& mappings = m_groupMeshes[mi]->dofMappings();
        for (uint32_t li = 0; li < static_cast<uint32_t>(mappings.size()); ++li) {
            const auto& dm = mappings[li];
            const Vector3d& pos = m_groupMeshes[mi]->vertex(dm.vertexIndex).position;
            Candidate c;
            c.meshIndex = static_cast<uint32_t>(mi);
            c.localDofIndex = li;
            // Quantize to ~1 µm for spatial hash
            c.qx = static_cast<int64_t>(std::round(pos.x() * 1e6));
            c.qy = static_cast<int64_t>(std::round(pos.y() * 1e6));
            c.qz = static_cast<int64_t>(std::round(pos.z() * 1e6));
            c.axis = dm.axis;
            candidates.push_back(c);
        }
    }

    // Hash key combines position and axis
    auto makeKey = [](int64_t x, int64_t y, int64_t z, uint8_t axis) -> uint64_t {
        uint64_t k = static_cast<uint64_t>(x);
        k ^= (static_cast<uint64_t>(y) << 21);
        k ^= (static_cast<uint64_t>(z) << 42);
        k ^= (static_cast<uint64_t>(axis) << 60);
        return k;
    };

    std::unordered_map<uint64_t, size_t> keyToGlobalIdx;

    for (const auto& c : candidates) {
        uint64_t key = makeKey(c.qx, c.qy, c.qz, c.axis);
        auto it = keyToGlobalIdx.find(key);
        if (it == keyToGlobalIdx.end()) {
            // First occurrence: create new global DOF
            keyToGlobalIdx[key] = m_globalDofRefs.size();
            m_globalDofRefs.push_back({{c.meshIndex, c.localDofIndex}});
        } else {
            // Duplicate: append to existing global DOF ref list
            m_globalDofRefs[it->second].push_back({c.meshIndex, c.localDofIndex});
        }
    }

    m_globalDofsDirty = false;
}

int LithologyModel::globalDofIndex(int meshIndex, int localDofIndex) const {
    if (m_globalDofsDirty) rebuildGlobalDofs();
    if (meshIndex < 0 || localDofIndex < 0) return -1;
    auto mi = static_cast<uint32_t>(meshIndex);
    auto li = static_cast<uint32_t>(localDofIndex);
    for (size_t g = 0; g < m_globalDofRefs.size(); ++g) {
        for (const auto& ref : m_globalDofRefs[g]) {
            if (ref.meshIndex == mi && ref.localDofIndex == li) {
                return static_cast<int>(g);
            }
        }
    }
    return -1;
}

// ---------------------------------------------------------------------------
// DOF / Parameter vector
// ---------------------------------------------------------------------------

uint32_t LithologyModel::totalDofCount() const {
    if (m_globalDofsDirty) rebuildGlobalDofs();
    return static_cast<uint32_t>(m_globalDofRefs.size());
}

VectorXd LithologyModel::assembleParameterVector() const {
    if (m_globalDofsDirty) rebuildGlobalDofs();
    VectorXd params(m_globalDofRefs.size());

    for (size_t g = 0; g < m_globalDofRefs.size(); ++g) {
        // Read from the first ref — all refs are kept in sync by applyParameterVector
        const auto& ref = m_globalDofRefs[g][0];
        const SurfaceMesh* mesh = m_groupMeshes[ref.meshIndex].get();
        const auto& dm = mesh->dofMappings()[ref.localDofIndex];
        const Vector3d& pos = mesh->vertex(dm.vertexIndex).position;
        double s = mesh->axisScale(dm.axis);
        params[static_cast<Eigen::Index>(g)] = pos[dm.axis] / s;
    }
    return params;
}

void LithologyModel::applyParameterVector(const VectorXd& params) {
    if (m_globalDofsDirty) rebuildGlobalDofs();

    for (size_t g = 0; g < m_globalDofRefs.size(); ++g) {
        double val = params[static_cast<Eigen::Index>(g)];
        for (const auto& ref : m_globalDofRefs[g]) {
            SurfaceMesh* mesh = m_groupMeshes[ref.meshIndex].get();
            const auto& dm = mesh->dofMappings()[ref.localDofIndex];
            Vertex& v = mesh->vertex(dm.vertexIndex);
            double s = mesh->axisScale(dm.axis);
            v.position[dm.axis] = val * s;
        }
    }

    // Interpolate non-control-point vertices for meshes with stride active
    for (auto& mesh : m_groupMeshes) {
        if (mesh && mesh->controlPointStride() > 0) {
            mesh->interpolateFromControlPoints();
        }
    }

    // Synchronize shared vertices across meshes so adjacent group boundaries
    // stay gap-free.  For each shared spatial position we compute the average
    // across all mesh instances and write it back, eliminating any divergence
    // introduced by mesh-dependent CP interpolation.
    for (const auto& group : m_sharedVertexSync) {
        if (group.size() < 2) continue;

        // Compute average position across all instances
        Vector3d avg = Vector3d::Zero();
        for (const auto& ref : group)
            avg += m_groupMeshes[ref.meshIndex]->vertex(ref.vertexIndex).position;
        avg /= static_cast<double>(group.size());

        // Write back
        for (const auto& ref : group)
            m_groupMeshes[ref.meshIndex]->vertex(ref.vertexIndex).position = avg;
    }
}

void LithologyModel::getBounds(VectorXd& lower, VectorXd& upper) const {
    if (m_globalDofsDirty) rebuildGlobalDofs();
    lower.resize(m_globalDofRefs.size());
    upper.resize(m_globalDofRefs.size());

    for (size_t g = 0; g < m_globalDofRefs.size(); ++g) {
        double lo = -1e12;
        double hi = 1e12;

        for (const auto& ref : m_globalDofRefs[g]) {
            const SurfaceMesh* mesh = m_groupMeshes[ref.meshIndex].get();
            const auto& dm = mesh->dofMappings()[ref.localDofIndex];
            double s = mesh->axisScale(dm.axis);

            if (dm.axis == 2) {
                // Z axis: use per-mesh depth bounds, take tightest intersection
                lo = std::max(lo, mesh->minDepth() / s);
                hi = std::min(hi, mesh->maxDepth() / s);
            } else if (dm.axis == 0) {
                // X axis: use model AABB if applyModelBounds was called
                lo = std::max(lo, m_modelBoundsMin.x() / s);
                hi = std::min(hi, m_modelBoundsMax.x() / s);
            } else if (dm.axis == 1) {
                // Y axis: use model AABB if applyModelBounds was called
                lo = std::max(lo, m_modelBoundsMin.y() / s);
                hi = std::min(hi, m_modelBoundsMax.y() / s);
            }
        }

        lower[static_cast<Eigen::Index>(g)] = lo;
        upper[static_cast<Eigen::Index>(g)] = hi;
    }
}

// ---------------------------------------------------------------------------
// Point classification
// ---------------------------------------------------------------------------

int LithologyModel::classifyPoint(const Vector3d& point) const {
    if (m_groupMeshes.empty() || m_groups.empty()) {
        return -1;
    }
    if (point.z() < m_bottomDepth) {
        return -1;
    }

    // Primary path: label grid O(1) lookup
    if (hasLabelGrid()) {
        int ix = static_cast<int>((point.x() - m_labelX0) / m_labelDx);
        int iy = static_cast<int>((point.y() - m_labelY0) / m_labelDy);
        int iz = static_cast<int>((point.z() - m_labelZ0) / m_labelDz);

        if (ix >= 0 && ix < m_labelNx && iy >= 0 && iy < m_labelNy && iz >= 0 && iz < m_labelNz) {
            int idx = iz * (m_labelNy * m_labelNx) + iy * m_labelNx + ix;
            int gid = m_labelGrid[idx];
            if (gid >= 0) {
                return gid;
            }
        }
    }

    // Fallback: ray casting (Möller-Trumbore) along +x axis.
    //
    // A tiny sub-micron offset is applied to y and z so the ray cannot
    // pass exactly through mesh edges/vertices.  This breaks the symmetry
    // that causes double-counting when a ray is aligned with a shared
    // edge (e.g. the boundary between two side-wall quads).
    const Vector3d rayDir(1.0, 0.0, 0.0);
    const Vector3d rayOrigin(point.x(),
                              point.y() + 1.234567e-6,
                              point.z() + 2.345678e-6);

    for (size_t gi = 0; gi < m_groupMeshes.size(); ++gi) {
        const SurfaceMesh* mesh = m_groupMeshes[gi].get();
        if (!mesh || mesh->triangleCount() == 0) continue;

        int hits = 0;
        for (uint32_t ti = 0; ti < mesh->triangleCount(); ++ti) {
            const Triangle& tri = mesh->triangle(ti);
            const Vector3d& v0 = mesh->vertex(tri.v0).position;
            const Vector3d& v1 = mesh->vertex(tri.v1).position;
            const Vector3d& v2 = mesh->vertex(tri.v2).position;

            Vector3d e1 = v1 - v0;
            Vector3d e2 = v2 - v0;
            Vector3d h = rayDir.cross(e2);
            double a = e1.dot(h);

            if (std::abs(a) < 1e-30) continue;

            double f = 1.0 / a;
            Vector3d s = rayOrigin - v0;
            double u = f * s.dot(h);

            if (u < 0.0 || u > 1.0) continue;

            Vector3d q = s.cross(e1);
            double v = f * rayDir.dot(q);

            if (v < 0.0 || u + v > 1.0) continue;

            double t = f * e2.dot(q);
            if (t > 1e-12) {
                ++hits;
            }
        }

        if (hits & 1) {
            return static_cast<int>(gi);
        }
    }

    return -1;
}

// ---------------------------------------------------------------------------
// fixExteriorFaces
// ---------------------------------------------------------------------------

void LithologyModel::fixExteriorFaces() {
    const int nM = groupMeshCount();
    if (nM < 1) return;

    // 1. Hash all triangles across all meshes by quantized vertex positions.
    //    A face is "shared" (interior) if the same triangle appears in ≥2 meshes.
    struct TriKey {
        int64_t x0, y0, z0, x1, y1, z1, x2, y2, z2;
        bool operator<(const TriKey& o) const {
            if (x0 != o.x0) return x0 < o.x0;
            if (y0 != o.y0) return y0 < o.y0;
            if (z0 != o.z0) return z0 < o.z0;
            if (x1 != o.x1) return x1 < o.x1;
            if (y1 != o.y1) return y1 < o.y1;
            if (z1 != o.z1) return z1 < o.z1;
            if (x2 != o.x2) return x2 < o.x2;
            if (y2 != o.y2) return y2 < o.y2;
            return z2 < o.z2;
        }
    };

    auto quantize = [](double v) -> int64_t {
        return static_cast<int64_t>(std::round(v * 1e6));
    };

    // Build canonically-ordered triangle key (sort vertex indices by quantized pos)
    auto makeTriKey = [&](const SurfaceMesh* m, uint32_t ti) -> TriKey {
        const auto& t = m->triangle(ti);
        int64_t qx[3] = { quantize(m->vertex(t.v0).position.x()),
                          quantize(m->vertex(t.v1).position.x()),
                          quantize(m->vertex(t.v2).position.x()) };
        int64_t qy[3] = { quantize(m->vertex(t.v0).position.y()),
                          quantize(m->vertex(t.v1).position.y()),
                          quantize(m->vertex(t.v2).position.y()) };
        int64_t qz[3] = { quantize(m->vertex(t.v0).position.z()),
                          quantize(m->vertex(t.v1).position.z()),
                          quantize(m->vertex(t.v2).position.z()) };
        // Canonical ordering: sort by (z, y, x) so the key is independent
        // of triangle winding and vertex ordering.
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                if (qz[i] > qz[j] || (qz[i] == qz[j] && qy[i] > qy[j])
                    || (qz[i] == qz[j] && qy[i] == qy[j] && qx[i] > qx[j])) {
                    std::swap(qx[i], qx[j]);
                    std::swap(qy[i], qy[j]);
                    std::swap(qz[i], qz[j]);
                }
            }
        }
        return {qx[0], qy[0], qz[0], qx[1], qy[1], qz[1], qx[2], qy[2], qz[2]};
    };

    std::map<TriKey, int> triCount;
    std::vector<std::vector<TriKey>> meshTriKeys(nM);

    for (int mi = 0; mi < nM; ++mi) {
        const SurfaceMesh* mesh = m_groupMeshes[mi].get();
        if (!mesh) continue;
        meshTriKeys[mi].reserve(mesh->triangleCount());
        for (uint32_t ti = 0; ti < mesh->triangleCount(); ++ti) {
            TriKey key = makeTriKey(mesh, ti);
            triCount[key]++;
            meshTriKeys[mi].push_back(key);
        }
    }

    // 2. For each mesh, find exterior faces and classify their vertices.
    //    Also count how many exterior faces of each normal class each vertex
    //    belongs to — this lets us distinguish faces, edges, and corners.
    struct VertexFaceClass {
        int topBottom = 0;   // ±Z dominant
        int sideNorth = 0;   // +Y dominant
        int sideSouth = 0;   // -Y dominant
        int sideEast = 0;    // +X dominant
        int sideWest = 0;    // -X dominant
    };

    for (int mi = 0; mi < nM; ++mi) {
        SurfaceMesh* mesh = m_groupMeshes[mi].get();
        if (!mesh) continue;

        uint32_t nV = mesh->vertexCount();
        std::vector<VertexFaceClass> vfc(nV);

        for (uint32_t ti = 0; ti < mesh->triangleCount(); ++ti) {
            const TriKey& key = meshTriKeys[mi][ti];
            if (triCount[key] > 1) continue;  // interior face — skip

            // Exterior face: classify by dominant normal
            const auto& t = mesh->triangle(ti);
            Vector3d e1 = mesh->vertex(t.v1).position - mesh->vertex(t.v0).position;
            Vector3d e2 = mesh->vertex(t.v2).position - mesh->vertex(t.v0).position;
            Vector3d n = e1.cross(e2);
            double len = n.norm();
            if (len < 1e-30) continue;
            n /= len;

            double absX = std::abs(n.x());
            double absY = std::abs(n.y());
            double absZ = std::abs(n.z());

            auto updateVFC = [&](uint32_t vi) {
                if (absZ >= absX && absZ >= absY) {
                    vfc[vi].topBottom++;
                } else if (absY >= absX) {
                    if (n.y() > 0) vfc[vi].sideNorth++;
                    else           vfc[vi].sideSouth++;
                } else {
                    if (n.x() > 0) vfc[vi].sideEast++;
                    else           vfc[vi].sideWest++;
                }
            };
            updateVFC(t.v0);
            updateVFC(t.v1);
            updateVFC(t.v2);
        }

        // 3. Assign vertex freedoms based on face normal classes.
        //    Count normal *classes* (not individual faces) so that a vertex
        //    shared by 4 triangles on the same side wall is still treated as
        //    belonging to a single side-face class.
        for (uint32_t vi = 0; vi < nV; ++vi) {
            const auto& fc = vfc[vi];

            bool hasTopBottom = (fc.topBottom > 0);
            bool hasNS = (fc.sideNorth > 0 || fc.sideSouth > 0);
            bool hasEW = (fc.sideEast > 0 || fc.sideWest > 0);
            int nClass = (hasTopBottom ? 1 : 0) + (hasNS ? 1 : 0) + (hasEW ? 1 : 0);

            if (nClass == 0) continue;  // interior — keep current

            if (nClass >= 3) {
                // Corner: belongs to all three normal classes → FIXED
                mesh->setVertexFreedom(vi, VertexFreedom::FIXED);
            } else if (nClass == 2) {
                // Edge: belongs to two normal classes → 1 DOF along edge
                if (hasTopBottom && hasNS) {
                    mesh->setVertexFreedom(vi, VertexFreedom::X_ONLY);
                } else if (hasTopBottom && hasEW) {
                    mesh->setVertexFreedom(vi, VertexFreedom::Y_ONLY);
                } else {
                    mesh->setVertexFreedom(vi, VertexFreedom::Z_ONLY);
                }
            } else {
                // nClass == 1: single normal class.
                // Side-wall vertices: 1 DOF along the boundary to preserve
                // the outer envelope shape. Top/bottom: XY slide.
                if (hasTopBottom) {
                    mesh->setVertexFreedom(vi, VertexFreedom::XY_FREE);
                } else if (hasNS) {
                    mesh->setVertexFreedom(vi, VertexFreedom::X_ONLY);
                } else {
                    mesh->setVertexFreedom(vi, VertexFreedom::Y_ONLY);
                }
            }
        }

        // 3b. Correction pass: vertices incident to interior faces (shared
        //     with other meshes) are on contact surfaces.  Step 3 may have
        //     over-restricted them based solely on adjacent exterior faces.
        //     Restore Z freedom so the contact surface can move vertically
        //     during inversion.
        std::vector<bool> touchesInterior(nV, false);
        for (uint32_t ti = 0; ti < mesh->triangleCount(); ++ti) {
            if (triCount[meshTriKeys[mi][ti]] > 1) {
                const auto& t = mesh->triangle(ti);
                touchesInterior[t.v0] = true;
                touchesInterior[t.v1] = true;
                touchesInterior[t.v2] = true;
            }
        }
        for (uint32_t vi = 0; vi < nV; ++vi) {
            if (!touchesInterior[vi]) continue;

            // Vertices that touch the top/bottom model envelope must stay
            // locked in Z — the Earth's surface and model base are hard
            // constraints. Only vertices on side walls (X_ONLY, Y_ONLY, FIXED
            // from side+side corners) get Z freedom added.
            bool hasTopBottom = (vfc[vi].topBottom > 0);

            auto f = mesh->vertex(vi).freedom;
            switch (f) {
            case VertexFreedom::FIXED:
                if (!hasTopBottom)
                    mesh->setVertexFreedom(vi, VertexFreedom::XYZ_FREE);
                break;
            case VertexFreedom::X_ONLY:
                if (!hasTopBottom)
                    mesh->setVertexFreedom(vi, VertexFreedom::XZ_FREE);
                break;
            case VertexFreedom::Y_ONLY:
                if (!hasTopBottom)
                    mesh->setVertexFreedom(vi, VertexFreedom::YZ_FREE);
                break;
            case VertexFreedom::XY_FREE:
                if (!hasTopBottom)
                    mesh->setVertexFreedom(vi, VertexFreedom::XYZ_FREE);
                break;
            default:
                break;  // already has Z freedom (Z_ONLY, XZ_FREE, YZ_FREE, XYZ_FREE)
            }
        }
    }

    // Rebuild global DOFs since we changed vertex freedoms
    m_globalDofsDirty = true;

    // Build shared-vertex sync list so applyParameterVector keeps contact
    // surfaces gap-free even without control-point downsampling.
    buildSharedVertexSync();
}

// ---------------------------------------------------------------------------
// buildSharedVertexSync
// ---------------------------------------------------------------------------

void LithologyModel::buildSharedVertexSync() {
    const int nM = groupMeshCount();
    if (nM < 2) return;

    // Quantize all vertex positions across all meshes and group by position.
    struct PosKey {
        int64_t x, y, z;
        bool operator<(const PosKey& o) const {
            if (x != o.x) return x < o.x;
            if (y != o.y) return y < o.y;
            return z < o.z;
        }
    };
    auto quantize = [](double v) -> int64_t {
        return static_cast<int64_t>(std::round(v * 1e6));
    };

    std::map<PosKey, std::vector<SharedVertexRef>> posToVertices;
    for (int mi = 0; mi < nM; ++mi) {
        if (!m_groupMeshes[mi]) continue;
        const auto& verts = m_groupMeshes[mi]->vertices();
        for (uint32_t vi = 0; vi < static_cast<uint32_t>(verts.size()); ++vi) {
            const Vector3d& p = verts[vi].position;
            PosKey key{quantize(p.x()), quantize(p.y()), quantize(p.z())};
            posToVertices[key].push_back({static_cast<uint32_t>(mi), vi});
        }
    }

    m_sharedVertexSync.clear();
    for (auto& [key, refs] : posToVertices) {
        if (refs.size() < 2) continue;
        m_sharedVertexSync.push_back(std::move(refs));
    }
}

// ---------------------------------------------------------------------------
// validate — model quality checks
// ---------------------------------------------------------------------------

ValidationResult LithologyModel::validate(size_t overlapBaseline,
                                         size_t adjacentOverlapBaseline) const {
    ValidationResult result;
    const int nM = groupMeshCount();
    if (nM < 2) return result;

    // Compute model AABB and characteristic dimension (reused by several checks)
    double xmin = 1e30, xmax = -1e30;
    double ymin = 1e30, ymax = -1e30;
    double zmin = 1e30, zmax = -1e30;
    for (int g = 0; g < nM; ++g) {
        const auto* mesh = groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const auto& p = mesh->vertex(vi).position;
            if (p.x() < xmin) xmin = p.x();
            if (p.x() > xmax) xmax = p.x();
            if (p.y() < ymin) ymin = p.y();
            if (p.y() > ymax) ymax = p.y();
            if (p.z() < zmin) zmin = p.z();
            if (p.z() > zmax) zmax = p.z();
        }
    }
    double modelDim = std::max({xmax - xmin, ymax - ymin, zmax - zmin, 1.0});
    double aabbMargin = modelDim * 0.02;

    // --- Check 1: Shared vertex gaps ---
    for (const auto& group : m_sharedVertexSync) {
        if (group.size() < 2) continue;
        const auto& refPos = m_groupMeshes[group[0].meshIndex]
                                  ->vertex(group[0].vertexIndex).position;
        for (size_t i = 1; i < group.size(); ++i) {
            const auto& pos = m_groupMeshes[group[i].meshIndex]
                                  ->vertex(group[i].vertexIndex).position;
            double gap = (pos - refPos).norm();
            if (gap > 0.1) {
                result.passed = false;
                result.worstViolation = std::max(result.worstViolation, gap);
                result.failureReason = "Shared vertex gap: " + std::to_string(gap) + "m";
                return result;
            }
        }
    }

    // --- Check 2: Adjacent overlap (vertex of group N inside group N±1) ---
    size_t adjOverlaps = 0;
    for (int g = 0; g < nM; ++g) {
        const auto* mesh = groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            int classified = classifyPoint(mesh->vertex(vi).position);
            if (classified < 0 || classified == g) continue;
            if (std::abs(classified - g) == 1) {
                ++adjOverlaps;
            }
        }
    }
    size_t adjLimit = adjacentOverlapBaseline + 20;
    if (adjOverlaps > adjLimit) {
        result.passed = false;
        result.worstViolation = static_cast<double>(adjOverlaps - adjacentOverlapBaseline);
        result.failureReason = "Adjacent overlap: " + std::to_string(adjOverlaps)
            + " (baseline=" + std::to_string(adjacentOverlapBaseline)
            + ", limit=" + std::to_string(adjLimit) + ")";
        return result;
    }

    // --- Check 3: Non-adjacent overlap ---
    size_t nonAdjOverlaps = 0;
    for (int g = 0; g < nM; ++g) {
        const auto* mesh = groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            int classified = classifyPoint(mesh->vertex(vi).position);
            if (classified < 0 || classified == g) continue;
            if (std::abs(classified - g) > 1) {
                ++nonAdjOverlaps;
            }
        }
    }
    size_t overlapLimit = overlapBaseline + 120;
    if (nonAdjOverlaps > overlapLimit) {
        result.passed = false;
        result.worstViolation = static_cast<double>(nonAdjOverlaps - overlapBaseline);
        result.failureReason = "Non-adjacent overlap: " + std::to_string(nonAdjOverlaps)
            + " (baseline=" + std::to_string(overlapBaseline) + ", limit="
            + std::to_string(overlapLimit) + ")";
        return result;
    }

    // --- Check 4: Edge tearing (boundary vertices at shared positions) ---
    // Quantize positions of boundary vertices (freedom != XYZ_FREE) and
    // verify co-location across groups.
    struct PosKey {
        int64_t x, y, z;
        bool operator<(const PosKey& o) const {
            if (x != o.x) return x < o.x;
            if (y != o.y) return y < o.y;
            return z < o.z;
        }
    };
    auto quantize = [](double v) -> int64_t {
        return static_cast<int64_t>(std::round(v * 1e4));  // 0.1mm
    };
    std::map<PosKey, std::vector<Vector3d>> boundaryPositions;
    for (int g = 0; g < nM; ++g) {
        const auto* mesh = groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            if (mesh->vertex(vi).freedom == VertexFreedom::XYZ_FREE) continue;
            const auto& p = mesh->vertex(vi).position;
            PosKey key{quantize(p.x()), quantize(p.y()), quantize(p.z())};
            boundaryPositions[key].push_back(p);
        }
    }
    for (const auto& [key, positions] : boundaryPositions) {
        if (positions.size() < 2) continue;
        const auto& ref = positions[0];
        for (size_t i = 1; i < positions.size(); ++i) {
            double gap = (positions[i] - ref).norm();
            if (gap > 0.1) {
                result.passed = false;
                result.worstViolation = std::max(result.worstViolation, gap);
                result.failureReason = "Edge tearing: boundary vertices diverged by "
                    + std::to_string(gap) + "m";
                return result;
            }
        }
    }

    // --- Check 5: AABB containment ---
    for (int g = 0; g < nM; ++g) {
        const auto* mesh = groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const auto& p = mesh->vertex(vi).position;
            if (p.x() < xmin - aabbMargin || p.x() > xmax + aabbMargin ||
                p.y() < ymin - aabbMargin || p.y() > ymax + aabbMargin ||
                p.z() < zmin - aabbMargin || p.z() > zmax + aabbMargin) {
                result.passed = false;
                result.worstViolation = std::max({
                    xmin - aabbMargin - p.x(), p.x() - xmax - aabbMargin,
                    ymin - aabbMargin - p.y(), p.y() - ymax - aabbMargin,
                    zmin - aabbMargin - p.z(), p.z() - zmax - aabbMargin, 0.0});
                result.failureReason = "AABB containment: vertex in group "
                    + std::to_string(g) + " outside model envelope";
                return result;
            }
        }
    }

    // --- Check 6: Spike detection (Laplacian magnitude > modelDim/5) ---
    double spikeThreshold = modelDim / 5.0;
    for (int g = 0; g < nM; ++g) {
        const auto* mesh = groupMesh(g);
        if (!mesh) continue;
        // Neighbors must be built for Laplacian computation
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const auto& nb = mesh->neighborVertices(vi);
            if (nb.empty()) continue;
            Vector3d avg = Vector3d::Zero();
            for (uint32_t nvi : nb)
                avg += mesh->vertex(nvi).position;
            avg /= static_cast<double>(nb.size());
            Vector3d laplacian = mesh->vertex(vi).position - avg;
            double mag = laplacian.norm();
            if (mag > spikeThreshold) {
                result.passed = false;
                result.worstViolation = mag;
                result.failureReason = "Spike detected: vertex in group "
                    + std::to_string(g) + " vi=" + std::to_string(vi)
                    + " Laplacian=" + std::to_string(mag) + "m (threshold="
                    + std::to_string(spikeThreshold) + "m)";
                return result;
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// applyModelBounds
// ---------------------------------------------------------------------------

void LithologyModel::applyModelBounds(double marginMultiplier) {
    const int nM = groupMeshCount();
    if (nM < 1) return;

    // 1. Compute model AABB from all vertices in all group meshes
    double minX = 1e30, maxX = -1e30;
    double minY = 1e30, maxY = -1e30;
    double minZ = 1e30, maxZ = -1e30;

    for (int mi = 0; mi < nM; ++mi) {
        const SurfaceMesh* mesh = m_groupMeshes[mi].get();
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const Vector3d& p = mesh->vertex(vi).position;
            if (p.x() < minX) minX = p.x();
            if (p.x() > maxX) maxX = p.x();
            if (p.y() < minY) minY = p.y();
            if (p.y() > maxY) maxY = p.y();
            if (p.z() < minZ) minZ = p.z();
            if (p.z() > maxZ) maxZ = p.z();
        }
    }

    // 2. Compute margin: 1% of max dimension
    double dimX = maxX - minX;
    double dimY = maxY - minY;
    double dimZ = maxZ - minZ;
    double margin = std::max({dimX, dimY, dimZ, 1.0}) * 0.01 * marginMultiplier;

    m_modelBoundsMin = Vector3d(minX - margin, minY - margin, minZ - margin);
    m_modelBoundsMax = Vector3d(maxX + margin, maxY + margin, maxZ + margin);
    m_modelBoundsMargin = margin;

    // 3. Set per-mesh Z bounds
    for (int mi = 0; mi < nM; ++mi) {
        SurfaceMesh* mesh = m_groupMeshes[mi].get();
        if (!mesh) continue;
        mesh->setBounds(minZ - margin, maxZ + margin);
    }

    m_globalDofsDirty = true;
}

// ---------------------------------------------------------------------------
// Depth extent
// ---------------------------------------------------------------------------

void LithologyModel::setBottomDepth(double depth) {
    m_bottomDepth = depth;
}

double LithologyModel::bottomDepth() const {
    return m_bottomDepth;
}

uint32_t LithologyModel::setControlPointStride(int stride) {
    m_controlPointStride = stride;
    m_globalDofsDirty = true;

    if (stride <= 0) {
        for (auto& mesh : m_groupMeshes) {
            if (mesh) mesh->setControlPointStride(0);
        }
        return totalDofCount();
    }

    // --- Mesh-independent control-point selection ---
    // FPS (used by the unstructured path in SurfaceMesh) selects CPs
    // independently per mesh, so shared contact vertices are rarely both
    // chosen.  We instead collect all vertex positions across all meshes,
    // build a regular 3D grid on the common voxel lattice, and select CPs
    // at every stride-th grid position.  Every mesh that owns a vertex at
    // a CP position gets that CP, keeping contact faces in lock-step.

    // 1. Collect unique quantized positions across all meshes.
    //    Map key (qx, qy, qz) → list of (meshIndex, vertexIndex).
    struct PosKey {
        int64_t x, y, z;
        bool operator<(const PosKey& o) const {
            if (x != o.x) return x < o.x;
            if (y != o.y) return y < o.y;
            return z < o.z;
        }
    };
    std::map<PosKey, std::vector<std::pair<uint32_t, uint32_t>>> posToVertices;

    for (size_t mi = 0; mi < m_groupMeshes.size(); ++mi) {
        if (!m_groupMeshes[mi]) continue;
        const auto& verts = m_groupMeshes[mi]->vertices();
        for (uint32_t vi = 0; vi < static_cast<uint32_t>(verts.size()); ++vi) {
            const Vector3d& p = verts[vi].position;
            PosKey key{
                static_cast<int64_t>(std::round(p.x() * 1e6)),
                static_cast<int64_t>(std::round(p.y() * 1e6)),
                static_cast<int64_t>(std::round(p.z() * 1e6))};
            posToVertices[key].push_back({static_cast<uint32_t>(mi), vi});
        }
    }

    if (posToVertices.size() < 4) {
        m_controlPointStride = 0;
        return 0;
    }

    // 2. Estimate the voxel grid spacing from the minimum non-zero
    //    coordinate difference along each axis.  We sample the first
    //    few hundred unique positions to keep this cheap.
    double minDx = 1e12, minDy = 1e12, minDz = 1e12;
    {
        std::vector<PosKey> keys;
        keys.reserve(posToVertices.size());
        for (const auto& [k, _] : posToVertices)
            keys.push_back(k);
        size_t nSample = std::min(keys.size(), size_t(300));
        for (size_t i = 0; i < nSample; ++i) {
            for (size_t j = i + 1; j < nSample; ++j) {
                double dx = std::abs(static_cast<double>(keys[i].x - keys[j].x)) * 1e-6;
                double dy = std::abs(static_cast<double>(keys[i].y - keys[j].y)) * 1e-6;
                double dz = std::abs(static_cast<double>(keys[i].z - keys[j].z)) * 1e-6;
                if (dx > 0 && dx < minDx) minDx = dx;
                if (dy > 0 && dy < minDy) minDy = dy;
                if (dz > 0 && dz < minDz) minDz = dz;
            }
        }
    }
    double cellSize = std::max({minDx, minDy, minDz, 1.0});     // ≥ 1 m
    double cellSizeQ = cellSize * 1e6;                           // in µm

    // 3. Determine grid extents (all-mesh bounding box).
    int64_t minQX = std::numeric_limits<int64_t>::max();
    int64_t maxQX = std::numeric_limits<int64_t>::min();
    int64_t minQY = minQX, maxQY = maxQX;
    int64_t minQZ = minQX, maxQZ = maxQX;
    for (const auto& [k, _] : posToVertices) {
        if (k.x < minQX) minQX = k.x; if (k.x > maxQX) maxQX = k.x;
        if (k.y < minQY) minQY = k.y; if (k.y > maxQY) maxQY = k.y;
        if (k.z < minQZ) minQZ = k.z; if (k.z > maxQZ) maxQZ = k.z;
    }

    // 4. Unified grid sampling: all positions (shared and non-shared) are
    //    assigned to grid cells.  One CP per cell, chosen closest to the cell
    //    centre.  When the selected position is shared (appears in ≥2 meshes),
    //    ALL its mesh instances become CPs, so contact vertices stay in lock-step.
    int64_t strideQ = static_cast<int64_t>(stride) * static_cast<int64_t>(cellSizeQ);
    if (strideQ < 1) strideQ = 1;

    auto gridCell = [strideQ](int64_t q) -> int64_t {
        if (q >= 0) return q / strideQ;
        return (q - strideQ + 1) / strideQ;
    };

    struct CellKey {
        int64_t cx, cy, cz;
        bool operator<(const CellKey& o) const {
            if (cx != o.cx) return cx < o.cx;
            if (cy != o.cy) return cy < o.cy;
            return cz < o.cz;
        }
    };

    struct CellCandidate {
        PosKey pos;
        int64_t distSq;
    };
    std::map<CellKey, CellCandidate> cellToBest;

    for (const auto& [k, verts] : posToVertices) {
        CellKey ck{gridCell(k.x), gridCell(k.y), gridCell(k.z)};
        int64_t cxCenter = ck.cx * strideQ + strideQ / 2;
        int64_t cyCenter = ck.cy * strideQ + strideQ / 2;
        int64_t czCenter = ck.cz * strideQ + strideQ / 2;
        int64_t dSq = (k.x - cxCenter) * (k.x - cxCenter)
                    + (k.y - cyCenter) * (k.y - cyCenter)
                    + (k.z - czCenter) * (k.z - czCenter);

        auto it = cellToBest.find(ck);
        if (it == cellToBest.end()) {
            cellToBest[ck] = {k, dSq};
        } else if (dSq < it->second.distSq) {
            it->second = {k, dSq};
        }
    }

    std::set<PosKey> cpPositions;
    size_t sharedCount = 0;
    size_t nonsharedCount = 0;
    for (const auto& [ck, candidate] : cellToBest) {
        cpPositions.insert(candidate.pos);
        auto it = posToVertices.find(candidate.pos);
        if (it != posToVertices.end() && it->second.size() > 1)
            ++sharedCount;
        else
            ++nonsharedCount;
    }

    std::cout << "[DOF CP] " << posToVertices.size() << " unique positions, "
              << "cellSize=" << cellSize << "m, "
              << cellToBest.size() << " populated cells, "
              << sharedCount << " shared + "
              << nonsharedCount << " non-shared = "
              << cpPositions.size() << " total CP positions" << std::endl;

    // 6. Distribute CP indices to each mesh.
    std::vector<std::vector<uint32_t>> meshCPs(m_groupMeshes.size());
    for (const auto& pk : cpPositions) {
        auto it = posToVertices.find(pk);
        if (it == posToVertices.end()) continue;
        for (const auto& [mi, vi] : it->second)
            meshCPs[mi].push_back(vi);
    }

    // 7. Apply per-mesh CP sets.
    for (size_t mi = 0; mi < m_groupMeshes.size(); ++mi) {
        if (!m_groupMeshes[mi]) continue;
        if (meshCPs[mi].empty()) continue;
        m_groupMeshes[mi]->setControlPointsExplicit(meshCPs[mi]);
    }

    // 8. Build shared-vertex sync list via the shared method (also called
    //    unconditionally from fixExteriorFaces).
    buildSharedVertexSync();

    std::cout << "[DOF CP] " << posToVertices.size() << " unique positions, "
              << "cellSize=" << cellSize << "m, "
              << cellToBest.size() << " populated cells, "
              << sharedCount << " shared + "
              << nonsharedCount << " non-shared = "
              << cpPositions.size() << " total CP positions, "
              << m_sharedVertexSync.size() << " shared sync groups" << std::endl;

    return totalDofCount();
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

bool LithologyModel::isValid() const {
    if (m_groups.size() < 2) return false;
    if (m_groupMeshes.size() != m_groups.size()) return false;

    for (const auto& mesh : m_groupMeshes) {
        if (!mesh || !mesh->isValid()) return false;
    }
    for (const auto& g : m_groups) {
        if (!g.isValid()) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Group proposal hook
// ---------------------------------------------------------------------------

void LithologyModel::setGroupProposal(std::shared_ptr<GroupProposal> proposal) {
    m_groupProposal = std::move(proposal);
}

} // namespace litho_invert

