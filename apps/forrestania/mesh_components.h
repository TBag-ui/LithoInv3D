#pragma once
#include <litho_invert/surface/surface_mesh.h>
#include <vector>
#include <memory>
#include <unordered_map>
#include <queue>

namespace litho_invert {

// Find connected components of a triangle mesh based on shared vertices.
// Returns one component label per triangle (0, 1, 2, ...).
inline std::vector<int> findConnectedComponents(const SurfaceMesh& mesh) {
    uint32_t nTri = mesh.triangleCount();
    std::vector<int> comp(nTri, -1);
    if (nTri == 0) return comp;

    // Build vertex→triangle adjacency
    std::vector<std::vector<uint32_t>> vertTris(mesh.vertexCount());
    for (uint32_t ti = 0; ti < nTri; ++ti) {
        const auto& t = mesh.triangle(ti);
        vertTris[t.v0].push_back(ti);
        vertTris[t.v1].push_back(ti);
        vertTris[t.v2].push_back(ti);
    }

    int curComp = 0;
    for (uint32_t seed = 0; seed < nTri; ++seed) {
        if (comp[seed] >= 0) continue;
        // BFS from seed triangle
        std::queue<uint32_t> q;
        q.push(seed);
        comp[seed] = curComp;
        while (!q.empty()) {
            uint32_t ti = q.front(); q.pop();
            const auto& t = mesh.triangle(ti);
            // Visit all triangles sharing any vertex with this triangle
            for (uint32_t vi : {t.v0, t.v1, t.v2}) {
                for (uint32_t nbr : vertTris[vi]) {
                    if (comp[nbr] < 0) {
                        comp[nbr] = curComp;
                        q.push(nbr);
                    }
                }
            }
        }
        ++curComp;
    }
    return comp;
}

// Extract a connected component into its own SurfaceMesh.
// compLabels: per-triangle component labels from findConnectedComponents().
// compIdx: which component to extract.
inline std::shared_ptr<SurfaceMesh> extractComponentMesh(
    const SurfaceMesh& src,
    const std::vector<int>& compLabels,
    int compIdx) {

    uint32_t nTri = src.triangleCount();
    uint32_t nVert = src.vertexCount();

    // Collect triangles belonging to this component
    std::vector<uint32_t> compTris;
    for (uint32_t ti = 0; ti < nTri; ++ti) {
        if (compLabels[ti] == compIdx) compTris.push_back(ti);
    }
    if (compTris.empty()) return nullptr;

    // Map old vertex indices → new vertex indices
    std::vector<int> oldToNew(nVert, -1);
    std::vector<uint32_t> newToOld; // for vertex property copying
    for (uint32_t ti : compTris) {
        const auto& t = src.triangle(ti);
        for (uint32_t vi : {t.v0, t.v1, t.v2}) {
            if (oldToNew[vi] < 0) {
                oldToNew[vi] = static_cast<int>(newToOld.size());
                newToOld.push_back(vi);
            }
        }
    }

    auto out = std::make_shared<SurfaceMesh>();
    out->setName(src.name());

    // Copy vertices
    for (uint32_t oldVi : newToOld) {
        const auto& v = src.vertex(oldVi);
        out->addVertex(v.position, v.freedom, v.moveVector);
    }

    // Copy triangles with reindexed vertices
    for (uint32_t ti : compTris) {
        const auto& t = src.triangle(ti);
        out->addTriangle(
            oldToNew[t.v0],
            oldToNew[t.v1],
            oldToNew[t.v2]);
    }

    // Copy per-axis scales
    for (int a = 0; a < 3; ++a) {
        out->setAxisScale(a, src.axisScale(a));
    }
    // Copy depth bounds
    out->setBounds(src.minDepth(), src.maxDepth());

    out->buildNeighbors();
    return out;
}

// Split a mesh into its disconnected components.
// Returns a vector of meshes, one per connected component.
inline std::vector<std::shared_ptr<SurfaceMesh>> splitDisconnectedComponents(
    const SurfaceMesh& mesh) {
    auto comps = findConnectedComponents(mesh);
    if (comps.empty()) return {};

    // Count components
    int numComps = 0;
    for (int c : comps) numComps = std::max(numComps, c + 1);

    std::vector<std::shared_ptr<SurfaceMesh>> result;
    for (int c = 0; c < numComps; ++c) {
        auto m = extractComponentMesh(mesh, comps, c);
        if (m) result.push_back(m);
    }
    return result;
}

} // namespace litho_invert

