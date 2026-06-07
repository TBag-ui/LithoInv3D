#include <litho_invert/surface/surface_mesh.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>

namespace litho_invert {

// ==================== Vertex management ====================

uint32_t SurfaceMesh::addVertex(const Vector3d& pos, VertexFreedom freedom) {
    m_vertices.emplace_back(pos, freedom);
    m_dofsDirty = true;
    return static_cast<uint32_t>(m_vertices.size() - 1);
}

uint32_t SurfaceMesh::addVertex(const Vector3d& pos, VertexFreedom freedom, const Vector3d& moveVector) {
    m_vertices.emplace_back(pos, freedom, moveVector);
    m_dofsDirty = true;
    return static_cast<uint32_t>(m_vertices.size() - 1);
}

uint32_t SurfaceMesh::addVertex(double x, double y, double z, VertexFreedom freedom) {
    m_vertices.emplace_back(Vector3d(x, y, z), freedom);
    m_dofsDirty = true;
    return static_cast<uint32_t>(m_vertices.size() - 1);
}

const Vertex& SurfaceMesh::vertex(uint32_t i) const {
    return m_vertices[i];
}

Vertex& SurfaceMesh::vertex(uint32_t i) {
    return m_vertices[i];
}

uint32_t SurfaceMesh::vertexCount() const {
    return static_cast<uint32_t>(m_vertices.size());
}

const std::vector<Vertex>& SurfaceMesh::vertices() const {
    return m_vertices;
}

void SurfaceMesh::setAllFreedom(VertexFreedom f) {
    for (auto& v : m_vertices) {
        v.freedom = f;
    }
    rebuildDofMappings();
}

void SurfaceMesh::setVertexFreedom(uint32_t i, VertexFreedom f, const Vector3d& moveVec) {
    m_vertices[i].freedom = f;
    m_vertices[i].moveVector = moveVec;
    m_dofsDirty = true;
}

// ==================== Triangle management ====================

uint32_t SurfaceMesh::addTriangle(uint32_t v0, uint32_t v1, uint32_t v2) {
    m_triangles.emplace_back(v0, v1, v2);
    return static_cast<uint32_t>(m_triangles.size() - 1);
}

const Triangle& SurfaceMesh::triangle(uint32_t i) const {
    return m_triangles[i];
}

uint32_t SurfaceMesh::triangleCount() const {
    return static_cast<uint32_t>(m_triangles.size());
}

const std::vector<Triangle>& SurfaceMesh::triangles() const {
    return m_triangles;
}

// ==================== DOF / Parameter vector ====================

uint32_t SurfaceMesh::dofCount() const {
    ensureDofsBuilt();
    return static_cast<uint32_t>(m_dofMappings.size());
}

const std::vector<SurfaceMesh::DofMapping>& SurfaceMesh::dofMappings() const {
    ensureDofsBuilt();
    return m_dofMappings;
}

void SurfaceMesh::rebuildDofMappings() {
    m_dofMappings.clear();
    uint32_t paddingSkipped = 0;
    uint32_t cpSkipped = 0;
    uint32_t fixedSkipped = 0;
    uint32_t dofsAdded = 0;
    for (uint32_t vi = 0; vi < m_vertices.size(); ++vi) {
        // Padding vertices never contribute DOFs — they follow the interior
        if (isPaddingVertex(vi)) {
            ++paddingSkipped;
            continue;
        }
        // When control-point stride is active, only control points get DOFs
        if (m_controlPointStride > 0 && !isControlPoint(vi)) {
            ++cpSkipped;
            continue;
        }
        const Vertex& v = m_vertices[vi];
        switch (v.freedom) {
            case VertexFreedom::FIXED:
                ++fixedSkipped;
                break;
            case VertexFreedom::Z_ONLY:
                m_dofMappings.push_back({vi, 2});
                ++dofsAdded;
                break;
            case VertexFreedom::X_ONLY:
                m_dofMappings.push_back({vi, 0});
                ++dofsAdded;
                break;
            case VertexFreedom::Y_ONLY:
                m_dofMappings.push_back({vi, 1});
                ++dofsAdded;
                break;
            case VertexFreedom::ALONG_VECTOR:
                m_dofMappings.push_back({vi, 0});
                ++dofsAdded;
                break;
            case VertexFreedom::XY_FREE:
                m_dofMappings.push_back({vi, 0});
                m_dofMappings.push_back({vi, 1});
                dofsAdded += 2;
                break;
            case VertexFreedom::XZ_FREE:
                m_dofMappings.push_back({vi, 0});
                m_dofMappings.push_back({vi, 2});
                dofsAdded += 2;
                break;
            case VertexFreedom::YZ_FREE:
                m_dofMappings.push_back({vi, 1});
                m_dofMappings.push_back({vi, 2});
                dofsAdded += 2;
                break;
            case VertexFreedom::XYZ_FREE:
                m_dofMappings.push_back({vi, 0});
                m_dofMappings.push_back({vi, 1});
                m_dofMappings.push_back({vi, 2});
                dofsAdded += 3;
                break;
        }
    }
#ifdef LITHO_INVERT_DEBUG
    std::cerr << "[DEBUG rebuildDofMappings: vertices=" << m_vertices.size()
              << " stride=" << m_controlPointStride
              << " paddingRings=" << m_paddingRings
              << " paddingSkipped=" << paddingSkipped
              << " cpSkipped=" << cpSkipped
              << " fixedSkipped=" << fixedSkipped
              << " dofsAdded=" << dofsAdded
              << "]" << std::endl;
    int shown = 0;
    for (uint32_t vi = 0; vi < m_vertices.size() && shown < 5; ++vi) {
        if (m_controlPointStride > 0 && !isControlPoint(vi)) continue;
        if (isPaddingVertex(vi)) continue;
        int f = static_cast<int>(m_vertices[vi].freedom);
        std::cerr << "[DEBUG   cp vi=" << vi << " freedom=" << f
                  << " (0=FIXED 1=Z_ONLY 2=X_ONLY 3=Y_ONLY 4=VEC 5=XY 6=XZ 7=YZ 8=XYZ)]"
                  << std::endl;
        ++shown;
    }
#endif
}

void SurfaceMesh::ensureDofsBuilt() const {
    if (!m_dofsDirty) return;
    const_cast<SurfaceMesh*>(this)->rebuildDofMappings();
    m_dofsDirty = false;
}

void SurfaceMesh::setAxisScale(int axis, double scale) {
    if (axis >= 0 && axis < 3) m_axisScales[axis] = scale;
}

double SurfaceMesh::axisScale(int axis) const {
    return (axis >= 0 && axis < 3) ? m_axisScales[axis] : 1.0;
}

void SurfaceMesh::applyParams(const VectorXd& params, uint32_t offset) {
    ensureDofsBuilt();
    for (uint32_t i = 0; i < m_dofMappings.size(); ++i) {
        const DofMapping& mapping = m_dofMappings[i];
        Vertex& v = m_vertices[mapping.vertexIndex];

        if (v.freedom == VertexFreedom::ALONG_VECTOR) {
            Vector3d mv = v.moveVector.normalized();
            double proj = v.position.dot(mv);
            double newProj = params[offset + i];
            v.position += mv * (newProj - proj);
        } else {
            double s = m_axisScales[mapping.axis];
            v.position[mapping.axis] = params[offset + i] * s;
        }
    }

    // After setting control-point positions, interpolate non-control vertices
    if (m_controlPointStride > 0) {
        interpolateFromControlPoints();
    }
}

void SurfaceMesh::extractParams(VectorXd& target, uint32_t offset) const {
    ensureDofsBuilt();
    for (uint32_t i = 0; i < m_dofMappings.size(); ++i) {
        const DofMapping& mapping = m_dofMappings[i];
        const Vertex& v = m_vertices[mapping.vertexIndex];

        if (v.freedom == VertexFreedom::ALONG_VECTOR) {
            Vector3d mv = v.moveVector.normalized();
            target[offset + i] = v.position.dot(mv);
        } else {
            double s = m_axisScales[mapping.axis];
            target[offset + i] = v.position[mapping.axis] / s;
        }
    }
}

void SurfaceMesh::getBounds(VectorXd& lower, VectorXd& upper, uint32_t offset) const {
    ensureDofsBuilt();
    for (uint32_t i = 0; i < m_dofMappings.size(); ++i) {
        const DofMapping& mapping = m_dofMappings[i];
        const Vertex& v = m_vertices[mapping.vertexIndex];

        double s = m_axisScales[mapping.axis];
        if (v.freedom == VertexFreedom::ALONG_VECTOR) {
            // Basically unbounded for movement along the vector
            lower[offset + i] = -1e6;
            upper[offset + i] = 1e6;
        } else if (mapping.axis == 2) {
            // Z axis: use depth bounds (scaled)
            lower[offset + i] = m_minDepth / s;
            upper[offset + i] = m_maxDepth / s;
        } else {
            // X or Y axis: basically unbounded (scaled)
            lower[offset + i] = -1e6 / s;
            upper[offset + i] = 1e6 / s;
        }
    }
}

void SurfaceMesh::setBounds(double minZ, double maxZ) {
    m_minDepth = minZ;
    m_maxDepth = maxZ;
}

double SurfaceMesh::minDepth() const {
    return m_minDepth;
}

double SurfaceMesh::maxDepth() const {
    return m_maxDepth;
}

// ==================== Neighbor relationships ====================

void SurfaceMesh::buildNeighbors() {
    m_neighbors.clear();
    m_neighbors.resize(m_vertices.size());

    // Gather all unique neighbor pairs from triangles
    for (const auto& t : m_triangles) {
        // Each pair within a triangle are neighbors
        m_neighbors[t.v0].push_back(t.v1);
        m_neighbors[t.v0].push_back(t.v2);
        m_neighbors[t.v1].push_back(t.v0);
        m_neighbors[t.v1].push_back(t.v2);
        m_neighbors[t.v2].push_back(t.v0);
        m_neighbors[t.v2].push_back(t.v1);
    }

    // Deduplicate each vertex's neighbor list
    for (auto& nbrs : m_neighbors) {
        std::sort(nbrs.begin(), nbrs.end());
        auto last = std::unique(nbrs.begin(), nbrs.end());
        nbrs.erase(last, nbrs.end());
    }
}

const std::vector<uint32_t>& SurfaceMesh::neighborVertices(uint32_t vertexIndex) const {
    return m_neighbors[vertexIndex];
}

// ==================== Mesh validation ====================

bool SurfaceMesh::isValid() const {
    if (m_vertices.size() < 3) return false;
    if (m_triangles.empty()) return false;
    for (const auto& t : m_triangles) {
        if (t.v0 >= m_vertices.size() ||
            t.v1 >= m_vertices.size() ||
            t.v2 >= m_vertices.size()) {
            return false;
        }
    }
    return true;
}

// ==================== Naming ====================

void SurfaceMesh::setName(const std::string& name) {
    m_name = name;
}

const std::string& SurfaceMesh::name() const {
    return m_name;
}

// ==================== Control-point downsampling ====================

uint32_t SurfaceMesh::setControlPointStride(int stride) {
    m_controlPointStride = stride;
    m_controlPointIndices.clear();
    m_vertexToControlWeights.clear();

    if (stride <= 0) {
        m_dofsDirty = true;
        return 0;
    }

    uint32_t nVerts = static_cast<uint32_t>(m_vertices.size());
    int fullDim = static_cast<int>(std::sqrt(static_cast<double>(nVerts)));

    // Only support regular grids (N×N vertices)
    if (fullDim * fullDim != static_cast<int>(nVerts) || fullDim < 2) {
        // Unstructured mesh fallback (e.g. closed volume meshes from voxel extraction)
        return setControlPointStrideUnstructured(stride);
    }

    m_gridDim = fullDim;

    // Interior region in grid index space
    int i0 = m_paddingRings;
    int i1 = fullDim - m_paddingRings - 1;

    // Only interior vertices can be control points.
    // Padding vertices are always interpolated from interior CPs.
    std::vector<bool> isCP(nVerts, false);
    for (int r = i0; r <= i1; ++r) {
        int ri = ((r - i0) / stride) * stride + i0;
        if (ri > i1) ri = i1;
        if (r != ri && r != i0 && r != i1) continue;
        for (int c = i0; c <= i1; ++c) {
            int ci = ((c - i0) / stride) * stride + i0;
            if (ci > i1) ci = i1;
            if (c != ci && c != i0 && c != i1) continue;
            uint32_t vi = static_cast<uint32_t>(r * fullDim + c);
            isCP[vi] = true;
        }
    }

    // Build control point index list and mapping from vertex index to CP index
    std::vector<int> vertexToCPIndex(nVerts, -1);
    for (uint32_t vi = 0; vi < nVerts; ++vi) {
        if (isCP[vi]) {
            vertexToCPIndex[vi] = static_cast<int>(m_controlPointIndices.size());
            m_controlPointIndices.push_back(vi);
        }
    }

    // Build interpolation weights for all vertices
    m_vertexToControlWeights.resize(nVerts);

    for (uint32_t vi = 0; vi < nVerts; ++vi) {
        if (isCP[vi]) {
            m_vertexToControlWeights[vi].push_back({
                static_cast<uint32_t>(vertexToCPIndex[vi]), 1.0});
            continue;
        }

        int r = static_cast<int>(vi) / fullDim;
        int c = static_cast<int>(vi) % fullDim;

        // Find surrounding control points — clamp to interior region
        int r0, r1, c0, c1;

        if (r < i0) {
            r0 = r1 = i0;
        } else if (r > i1) {
            r0 = r1 = i1;
        } else {
            r0 = ((r - i0) / stride) * stride + i0;
            r1 = r0 + stride;
            if (r1 > i1) r1 = i1;
        }

        if (c < i0) {
            c0 = c1 = i0;
        } else if (c > i1) {
            c0 = c1 = i1;
        } else {
            c0 = ((c - i0) / stride) * stride + i0;
            c1 = c0 + stride;
            if (c1 > i1) c1 = i1;
        }

        double dr = (r1 > r0) ? static_cast<double>(r - r0) / static_cast<double>(r1 - r0) : 0.0;
        double dc = (c1 > c0) ? static_cast<double>(c - c0) / static_cast<double>(c1 - c0) : 0.0;

        double w00 = (1.0 - dr) * (1.0 - dc);
        double w01 = (1.0 - dr) * dc;
        double w10 = dr * (1.0 - dc);
        double w11 = dr * dc;

        auto addWeight = [&](int rr, int cc, double w) {
            if (w < 1e-15) return;
            uint32_t idx = static_cast<uint32_t>(rr * fullDim + cc);
            int cpIdx = vertexToCPIndex[idx];
            if (cpIdx >= 0) {
                m_vertexToControlWeights[vi].push_back({
                    static_cast<uint32_t>(cpIdx), w});
            }
        };

        addWeight(r0, c0, w00);
        addWeight(r0, c1, w01);
        addWeight(r1, c0, w10);
        addWeight(r1, c1, w11);
    }

    m_dofsDirty = true;
    return static_cast<uint32_t>(m_controlPointIndices.size());
}

uint32_t SurfaceMesh::setControlPointStrideUnstructured(int stride) {
    m_controlPointStride = stride;
    m_controlPointIndices.clear();
    m_vertexToControlWeights.clear();
    m_gridDim = 0;

    uint32_t nVerts = static_cast<uint32_t>(m_vertices.size());
    if (nVerts < 4 || stride < 2) {
        m_controlPointStride = 0;
        return 0;
    }

    // Target CP count: roughly nVerts / (stride^2), matching the reduction
    // factor of the regular-grid path. Clamp to [4, nVerts-1].
    uint32_t targetCPs = std::max(4u, nVerts / static_cast<uint32_t>(stride * stride));
    targetCPs = std::min(targetCPs, nVerts - 1);
    if (targetCPs >= nVerts) {
        m_controlPointStride = 0;
        return 0;
    }

    // --- Farthest-point sampling on Euclidean positions ---
    // Start with vertex 0, then iteratively pick the vertex farthest from
    // all existing CPs. This gives good spatial coverage on any mesh topology.
    std::vector<double> minDistToCP(nVerts, std::numeric_limits<double>::max());

    m_controlPointIndices.push_back(0);
    for (uint32_t vi = 0; vi < nVerts; ++vi) {
        double d = (m_vertices[vi].position - m_vertices[0].position).squaredNorm();
        minDistToCP[vi] = d;
    }

    while (m_controlPointIndices.size() < targetCPs) {
        uint32_t farthest = 0;
        double farthestDist = -1.0;
        for (uint32_t vi = 0; vi < nVerts; ++vi) {
            if (minDistToCP[vi] > farthestDist) {
                farthestDist = minDistToCP[vi];
                farthest = vi;
            }
        }
        if (farthestDist < 1e-20) break;  // degenerate mesh

        m_controlPointIndices.push_back(farthest);

        // Update minimum distances to the new CP
        const Vector3d& cpPos = m_vertices[farthest].position;
        for (uint32_t vi = 0; vi < nVerts; ++vi) {
            double d = (m_vertices[vi].position - cpPos).squaredNorm();
            if (d < minDistToCP[vi]) minDistToCP[vi] = d;
        }
    }

    uint32_t nCPs = static_cast<uint32_t>(m_controlPointIndices.size());
    if (nCPs < 4) {
        m_controlPointStride = 0;
        m_controlPointIndices.clear();
        return 0;
    }

    // Build CP vertex → CP index lookup
    std::unordered_map<uint32_t, uint32_t> cpVertexToIndex;
    for (uint32_t ci = 0; ci < nCPs; ++ci) {
        cpVertexToIndex[m_controlPointIndices[ci]] = ci;
    }

    // --- Interpolation weights: K nearest CPs per vertex ---
    // Non-CP vertices use inverse-distance-weighted interpolation from their
    // K nearest control points. CP vertices get weight 1.0 on themselves.
    constexpr int K = 4;
    m_vertexToControlWeights.resize(nVerts);

    for (uint32_t vi = 0; vi < nVerts; ++vi) {
        auto cpIt = cpVertexToIndex.find(vi);
        if (cpIt != cpVertexToIndex.end()) {
            m_vertexToControlWeights[vi].push_back({cpIt->second, 1.0});
            continue;
        }

        // Find K nearest CPs by Euclidean distance
        std::vector<std::pair<double, uint32_t>> nearest;  // (dist, cpIdx)
        nearest.reserve(nCPs);
        for (uint32_t ci = 0; ci < nCPs; ++ci) {
            uint32_t cpVi = m_controlPointIndices[ci];
            double d = (m_vertices[cpVi].position - m_vertices[vi].position).norm();
            nearest.push_back({d, ci});
        }

        // Partial sort: only need top K
        std::partial_sort(nearest.begin(),
                          nearest.begin() + std::min(static_cast<size_t>(K), nearest.size()),
                          nearest.end());

        // Inverse-distance weighting (Shepard interpolation)
        double totalW = 0.0;
        const double eps = 1e-12;
        size_t keep = std::min(static_cast<size_t>(K), nearest.size());
        for (size_t j = 0; j < keep; ++j) {
            totalW += 1.0 / std::max(nearest[j].first, eps);
        }

        for (size_t j = 0; j < keep; ++j) {
            double w = (1.0 / std::max(nearest[j].first, eps)) / totalW;
            m_vertexToControlWeights[vi].push_back({nearest[j].second, w});
        }
    }

    m_dofsDirty = true;
    return nCPs;
}

uint32_t SurfaceMesh::setControlPointsExplicit(const std::vector<uint32_t>& cpVertexIndices) {
    m_controlPointStride = 1;  // non-zero to activate CP mode
    m_controlPointIndices = cpVertexIndices;
    m_vertexToControlWeights.clear();
    m_gridDim = 0;

    uint32_t nVerts = static_cast<uint32_t>(m_vertices.size());
    uint32_t nCPs = static_cast<uint32_t>(m_controlPointIndices.size());
    if (nCPs < 2 || nVerts < 2) {
        m_controlPointStride = 0;
        m_controlPointIndices.clear();
        return 0;
    }

    // Build CP vertex → CP index lookup
    std::unordered_map<uint32_t, uint32_t> cpVertexToIndex;
    for (uint32_t ci = 0; ci < nCPs; ++ci) {
        cpVertexToIndex[m_controlPointIndices[ci]] = ci;
    }

    // Interpolation weights: K nearest CPs per vertex (same as unstructured path)
    constexpr int K = 4;
    m_vertexToControlWeights.resize(nVerts);

    for (uint32_t vi = 0; vi < nVerts; ++vi) {
        auto cpIt = cpVertexToIndex.find(vi);
        if (cpIt != cpVertexToIndex.end()) {
            m_vertexToControlWeights[vi].push_back({cpIt->second, 1.0});
            continue;
        }

        std::vector<std::pair<double, uint32_t>> nearest;
        nearest.reserve(nCPs);
        for (uint32_t ci = 0; ci < nCPs; ++ci) {
            uint32_t cpVi = m_controlPointIndices[ci];
            double d = (m_vertices[cpVi].position - m_vertices[vi].position).norm();
            nearest.push_back({d, ci});
        }

        std::partial_sort(nearest.begin(),
                          nearest.begin() + std::min(static_cast<size_t>(K), nearest.size()),
                          nearest.end());

        double totalW = 0.0;
        const double eps = 1e-12;
        size_t keep = std::min(static_cast<size_t>(K), nearest.size());
        for (size_t j = 0; j < keep; ++j)
            totalW += 1.0 / std::max(nearest[j].first, eps);

        for (size_t j = 0; j < keep; ++j) {
            double w = (1.0 / std::max(nearest[j].first, eps)) / totalW;
            m_vertexToControlWeights[vi].push_back({nearest[j].second, w});
        }
    }

    m_dofsDirty = true;
    return nCPs;
}

bool SurfaceMesh::isControlPoint(uint32_t vertexIndex) const {
    if (m_controlPointStride <= 0) return true;  // all vertices are CPs
    if (m_vertexToControlWeights.empty()) return true;
    return m_vertexToControlWeights[vertexIndex].size() == 1
        && m_vertexToControlWeights[vertexIndex][0].second == 1.0;
}

void SurfaceMesh::interpolateFromControlPoints() {
    if (m_controlPointStride <= 0) return;

    for (uint32_t vi = 0; vi < static_cast<uint32_t>(m_vertices.size()); ++vi) {
        if (isControlPoint(vi)) continue;

        Vector3d pos = Vector3d::Zero();
        for (const auto& [cpIndex, weight] : m_vertexToControlWeights[vi]) {
            uint32_t cpVertex = m_controlPointIndices[cpIndex];
            pos += weight * m_vertices[cpVertex].position;
        }
        m_vertices[vi].position = pos;
    }
}

void SurfaceMesh::downsampleVertexGradient(
    const std::vector<double>& fullVertexGrad,
    VectorXd& coarseGrad, uint32_t coarseOffset,
    int axis) const
{
    ensureDofsBuilt();

    // Build vertex → DOF index lookup, filtered by the requested axis.
    std::unordered_map<uint32_t, uint32_t> vertexToDof;
    for (uint32_t d = 0; d < static_cast<uint32_t>(m_dofMappings.size()); ++d) {
        if (static_cast<int>(m_dofMappings[d].axis) == axis) {
            vertexToDof[m_dofMappings[d].vertexIndex] = d;
        }
    }

    if (m_controlPointStride <= 0) {
        // No downsampling: directly map vertex gradients to DOFs
        for (const auto& [vertIdx, dofIdx] : vertexToDof) {
            if (vertIdx < static_cast<uint32_t>(fullVertexGrad.size())) {
                coarseGrad[coarseOffset + dofIdx] = fullVertexGrad[vertIdx];
            }
        }
        return;
    }

    for (uint32_t k = 0; k < static_cast<uint32_t>(m_controlPointIndices.size()); ++k) {
        uint32_t cpVertex = m_controlPointIndices[k];
        auto it = vertexToDof.find(cpVertex);
        if (it == vertexToDof.end()) continue;

        uint32_t dofIdx = it->second;
        double g = 0.0;
        for (uint32_t vi = 0; vi < static_cast<uint32_t>(m_vertices.size()); ++vi) {
            for (const auto& [cpIdx, weight] : m_vertexToControlWeights[vi]) {
                if (cpIdx == k) {
                    g += weight * fullVertexGrad[vi];
                    break;
                }
            }
        }
        coarseGrad[coarseOffset + dofIdx] = g;
    }
}

// ==================== Lateral padding ====================

void SurfaceMesh::setInteriorGrid(int interiorDim, int paddingRings) {
    m_interiorGridDim = interiorDim;
    m_paddingRings = paddingRings;
}

bool SurfaceMesh::isPaddingVertex(uint32_t vertexIndex) const {
    if (m_paddingRings <= 0) return false;
    int fullDim = m_interiorGridDim + 2 * m_paddingRings;
    int r = static_cast<int>(vertexIndex) / fullDim;
    int c = static_cast<int>(vertexIndex) % fullDim;
    return r < m_paddingRings
        || r >= m_paddingRings + m_interiorGridDim
        || c < m_paddingRings
        || c >= m_paddingRings + m_interiorGridDim;
}

bool SurfaceMesh::isInteriorVertex(uint32_t vertexIndex) const {
    if (m_paddingRings <= 0) return true;
    return !isPaddingVertex(vertexIndex);
}

void SurfaceMesh::extrapolatePadding(const std::vector<double>& upperZ,
                                      const std::vector<double>& lowerZ) {
    if (m_paddingRings <= 0) return;

    int fullDim = m_interiorGridDim + 2 * m_paddingRings;
    int p = m_paddingRings;
    int i0 = p;
    int i1 = p + m_interiorGridDim - 1;

    for (uint32_t vi = 0; vi < static_cast<uint32_t>(m_vertices.size()); ++vi) {
        if (!isPaddingVertex(vi)) continue;

        int r = static_cast<int>(vi) / fullDim;
        int c = static_cast<int>(vi) % fullDim;

        // Determine which side(s) and compute extrapolated z
        double zExtrap = 0.0;
        int contributing = 0;

        // Vertical extrapolation (top / bottom)
        if (r < p) {
            // Top padding: extrapolate upward using dip from interior boundary
            uint32_t bnd = static_cast<uint32_t>(i0 * fullDim + c);
            uint32_t inr = static_cast<uint32_t>((i0 + 1) * fullDim + c);
            double dz = m_vertices[bnd].position.z() - m_vertices[inr].position.z();
            double dist = static_cast<double>(p - r);
            zExtrap += m_vertices[bnd].position.z() + dist * dz;
            ++contributing;
        } else if (r > i1) {
            // Bottom padding: extrapolate downward
            uint32_t bnd = static_cast<uint32_t>(i1 * fullDim + c);
            uint32_t inr = static_cast<uint32_t>((i1 - 1) * fullDim + c);
            double dz = m_vertices[bnd].position.z() - m_vertices[inr].position.z();
            double dist = static_cast<double>(r - i1);
            zExtrap += m_vertices[bnd].position.z() + dist * dz;
            ++contributing;
        }

        // Horizontal extrapolation (left / right)
        if (c < p) {
            uint32_t bnd = static_cast<uint32_t>(r * fullDim + i0);
            uint32_t inr = static_cast<uint32_t>(r * fullDim + i0 + 1);
            double dz = m_vertices[bnd].position.z() - m_vertices[inr].position.z();
            double dist = static_cast<double>(p - c);
            zExtrap += m_vertices[bnd].position.z() + dist * dz;
            ++contributing;
        } else if (c > i1) {
            uint32_t bnd = static_cast<uint32_t>(r * fullDim + i1);
            uint32_t inr = static_cast<uint32_t>(r * fullDim + i1 - 1);
            double dz = m_vertices[bnd].position.z() - m_vertices[inr].position.z();
            double dist = static_cast<double>(c - i1);
            zExtrap += m_vertices[bnd].position.z() + dist * dz;
            ++contributing;
        }

        // For interior rows/cols (sides-only), the if-else above sets zExtrap once
        // For corners, average the two directional extrapolations
        if (contributing > 1) {
            zExtrap /= static_cast<double>(contributing);
        }

        // Apply pinching bounds if provided
        if (!upperZ.empty() && static_cast<size_t>(vi) < upperZ.size()) {
            zExtrap = std::min(zExtrap, upperZ[vi]);
        }
        if (!lowerZ.empty() && static_cast<size_t>(vi) < lowerZ.size()) {
            zExtrap = std::max(zExtrap, lowerZ[vi]);
        }

        m_vertices[vi].position.z() = zExtrap;
    }
}

// ==================== Face normals ====================

void SurfaceMesh::ensureFaceNormals() const {
    if (m_normalsComputed) return;
    const uint32_t nTris = triangleCount();
    m_faceNormals.resize(nTris);
    for (uint32_t ti = 0; ti < nTris; ++ti) {
        const Triangle& tri = triangle(ti);
        const Vector3d& p0 = vertex(tri.v0).position;
        const Vector3d& p1 = vertex(tri.v1).position;
        const Vector3d& p2 = vertex(tri.v2).position;
        Vector3d e1 = p1 - p0;
        Vector3d e2 = p2 - p0;
        Vector3d n = e1.cross(e2);
        double len = n.norm();
        if (len > 1e-30) {
            n /= len;
        }
        m_faceNormals[ti] = n;
    }
    m_normalsComputed = true;
}

// ==================== 3D side classification ====================

int SurfaceMesh::classifySide(const Vector3d& point) const {
    ensureFaceNormals();

    const uint32_t nTris = triangleCount();
    const uint32_t nVerts = vertexCount();
    if (nTris == 0 || nVerts == 0) return 0;

    // Step 1: try barycentric XY projection (fast path — handles most cases).
    // For height-field-like surfaces, the XY projection of the point lands on
    // the correct triangle.  We reuse the existing barycentric logic but use
    // the 3D triangle normal for side determination instead of Z comparison.
    for (uint32_t t = 0; t < nTris; ++t) {
        const Triangle& tri = triangle(t);
        const Vector3d& p0 = vertex(tri.v0).position;
        const Vector3d& p1 = vertex(tri.v1).position;
        const Vector3d& p2 = vertex(tri.v2).position;

        // Barycentric coordinates in XY plane
        double v0x = p2.x() - p0.x();
        double v0y = p2.y() - p0.y();
        double v1x = p1.x() - p0.x();
        double v1y = p1.y() - p0.y();
        double v2x = point.x() - p0.x();
        double v2y = point.y() - p0.y();

        double dot00 = v0x * v0x + v0y * v0y;
        double dot01 = v0x * v1x + v0y * v1y;
        double dot02 = v0x * v2x + v0y * v2y;
        double dot11 = v1x * v1x + v1y * v1y;
        double dot12 = v1x * v2x + v1y * v2y;

        double denom = dot00 * dot11 - dot01 * dot01;
        if (std::abs(denom) < 1e-30) continue;

        double invDenom = 1.0 / denom;
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        if (u >= -1e-12 && v >= -1e-12 && (u + v) <= 1.0 + 1e-12) {
            // Point projects onto this triangle in XY.
            // Use the 3D triangle normal for side determination.
            double w = 1.0 - u - v;
            Vector3d surfPt = w * p0 + v * p2 + u * p1;
            Vector3d toPoint = point - surfPt;
            double side = toPoint.dot(m_faceNormals[t]);
            if (side > 1e-9) return 1;
            if (side < -1e-9) return -1;
            return 0;
        }
    }

    // Step 2: point outside the XY extent — find the closest triangle in 3D
    // and use its normal for side determination.
    uint32_t bestTri = 0;
    double bestDist2 = std::numeric_limits<double>::max();
    Vector3d bestClosestPt;

    for (uint32_t ti = 0; ti < nTris; ++ti) {
        const Triangle& tri = triangle(ti);
        const Vector3d& p0 = vertex(tri.v0).position;
        const Vector3d& p1 = vertex(tri.v1).position;
        const Vector3d& p2 = vertex(tri.v2).position;

        Vector3d e0 = p1 - p0;
        Vector3d e1 = p2 - p0;
        Vector3d d = p0 - point;
        double a = e0.dot(e0);
        double b = e0.dot(e1);
        double c = e1.dot(e1);
        double det = a * c - b * b;
        if (std::abs(det) < 1e-30) continue;

        double u = b * e1.dot(d) - c * e0.dot(d);
        double v = b * e0.dot(d) - a * e1.dot(d);
        u /= det;
        v /= det;

        // Clamp to triangle
        if (u < 0.0) u = 0.0;
        if (v < 0.0) v = 0.0;
        if (u + v > 1.0) {
            double sum = u + v;
            u /= sum;
            v /= sum;
        }

        Vector3d closest = p0 + u * e0 + v * e1;
        double d2 = (point - closest).squaredNorm();
        if (d2 < bestDist2) {
            bestDist2 = d2;
            bestTri = ti;
            bestClosestPt = closest;
        }
    }

    Vector3d toPoint = point - bestClosestPt;
    double side = toPoint.dot(m_faceNormals[bestTri]);
    if (side > 1e-9) return 1;
    if (side < -1e-9) return -1;
    return 0;
}

} // namespace litho_invert
