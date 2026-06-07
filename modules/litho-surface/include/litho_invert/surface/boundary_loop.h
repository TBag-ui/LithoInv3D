#pragma once

#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace litho_invert {

// Forward declaration — full definition in surface_mesh.h
class SurfaceMesh;

/// A directed boundary edge in CCW perimeter order.
struct BoundaryEdge {
    uint32_t start;   // CCW start vertex
    uint32_t end;     // CCW end vertex
};

/// Find all boundary edges of a surface mesh.
///
/// A boundary edge appears in exactly one triangle. Interior edges appear
/// in two triangles (one from each adjacent quad). The returned edges are
/// directed in CCW order around the perimeter, so they can be used directly
/// for side-wall triangulation.
///
/// Used by: GravityForward, MagneticForward (polyhedron side walls)
std::vector<BoundaryEdge> findBoundaryEdges(const SurfaceMesh& mesh);

/// Walk boundary edges into ordered connected-component loops.
///
/// Each loop is a closed cycle of vertex indices. Multiple loops occur
/// when the surface has holes or disconnected components. Loops are
/// ordered CCW.
///
/// Used by: InversionExporter (closed-volume TS export)
std::vector<std::vector<uint32_t>> extractBoundaryLoops(const SurfaceMesh& mesh);

} // namespace litho_invert
