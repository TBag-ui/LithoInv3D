#include <litho_invert/surface/boundary_loop.h>
#include <litho_invert/surface/surface_mesh.h>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace litho_invert {

using Edge = std::pair<uint32_t, uint32_t>;

static Edge makeEdge(uint32_t a, uint32_t b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

std::vector<BoundaryEdge> findBoundaryEdges(const SurfaceMesh& mesh) {
    // Count how many triangles share each edge. Boundary edges appear once.
    std::map<Edge, int> edgeCount;
    std::map<Edge, bool> edgeIsAB; // true if v0→v1 is CCW in the mesh

    for (uint32_t t = 0; t < mesh.triangleCount(); ++t) {
        const Triangle& tri = mesh.triangle(t);

        auto processEdge = [&](uint32_t v0, uint32_t v1) {
            Edge key = makeEdge(v0, v1);
            if (edgeCount.find(key) == edgeCount.end()) {
                // First time seeing this edge — record CCW direction.
                // In a CCW-wound triangle, edges go v0→v1, v1→v2, v2→v0.
                // The sorted key (lo, hi) matches v0→v1 when v0 < v1.
                edgeIsAB[key] = (v0 < v1);
            }
            edgeCount[key]++;
        };

        processEdge(tri.v0, tri.v1);
        processEdge(tri.v1, tri.v2);
        processEdge(tri.v2, tri.v0);
    }

    std::vector<BoundaryEdge> result;
    for (const auto& [e, count] : edgeCount) {
        if (count != 1) continue; // interior edge

        uint32_t a = e.first;
        uint32_t b = e.second;
        bool ab = edgeIsAB[e];
        // CCW perimeter: start is the smaller index when edge was recorded a→b
        result.push_back({ab ? a : b, ab ? b : a});
    }
    return result;
}

std::vector<std::vector<uint32_t>> extractBoundaryLoops(const SurfaceMesh& mesh) {
    // Count edge occurrences
    std::map<Edge, int> edgeCount;
    for (uint32_t t = 0; t < mesh.triangleCount(); ++t) {
        const auto& tri = mesh.triangle(t);
        edgeCount[makeEdge(tri.v0, tri.v1)]++;
        edgeCount[makeEdge(tri.v1, tri.v2)]++;
        edgeCount[makeEdge(tri.v2, tri.v0)]++;
    }

    // Build adjacency from boundary edges only
    std::map<uint32_t, std::vector<uint32_t>> adj;
    for (const auto& [e, count] : edgeCount) {
        if (count == 1) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
    }

    // Walk connected components of the boundary graph into ordered loops
    std::set<uint32_t> visited;
    std::vector<std::vector<uint32_t>> loops;

    for (const auto& [v, nbs] : adj) {
        if (visited.count(v) || nbs.size() != 2) continue;

        std::vector<uint32_t> loop;
        uint32_t curr = v;
        uint32_t prev = UINT32_MAX;

        while (!visited.count(curr)) {
            loop.push_back(curr);
            visited.insert(curr);

            uint32_t next = UINT32_MAX;
            for (uint32_t nb : adj[curr]) {
                if (nb != prev) { next = nb; break; }
            }
            if (next == UINT32_MAX) break;
            prev = curr;
            curr = next;
        }

        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }
    return loops;
}

} // namespace litho_invert
