#pragma once
#include <litho_invert/surface/surface_mesh.h>
#include <memory>
#include <string>

namespace litho_invert {

class OBJSurfaceLoader {
public:
    // Read a surface mesh from a Wavefront OBJ file.
    // Supports v (vertex) and f (face/triangle) lines.
    // OBJ indices are 1-based; converted to 0-based internally.
    // Returns empty shared_ptr on parse failure.
    static std::shared_ptr<SurfaceMesh> load(const std::string& filepath);
};

class TSurfaceLoader {
public:
    // Read a surface mesh from a GOCAD TSurf (.ts) file.
    // Supports VRTX/PVRTX (vertex) and TRGL (triangle) lines.
    // Accepts both TFACE and TRIANGLES geometry keywords.
    // GOCAD indices are 1-based; converted to 0-based internally.
    // Returns empty shared_ptr on parse failure.
    static std::shared_ptr<SurfaceMesh> load(const std::string& filepath);
};

// Auto-detect loader by file extension: .ts -> TSurfaceLoader, .obj -> OBJSurfaceLoader.
// Returns nullptr for unrecognized extensions.
std::shared_ptr<SurfaceMesh> loadSurfaceMesh(const std::string& filepath);

} // namespace litho_invert

