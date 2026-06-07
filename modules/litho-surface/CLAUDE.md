# litho-surface — Triangulated Surface Mesh with DOF System

## Purpose

The fundamental geometric primitive of the inversion. Each boundary between adjacent
litho groups is a SurfaceMesh — a regular grid of vertices and triangles with a
configurable DOF system, control-point downsampling, and padding support.

## API

```cpp
#include <litho_invert/surface/surface_mesh.h>

// Vertex freedom levels
enum class VertexFreedom : uint8_t {
    FIXED, Z_ONLY, X_ONLY, Y_ONLY, ALONG_VECTOR, XY_FREE, XYZ_FREE
};

struct Vertex {
    Vector3d position;
    VertexFreedom freedom = VertexFreedom::FIXED;  // DEFAULT: locked
    Vector3d moveVector = Vector3d::UnitZ();
};

struct Triangle {
    uint32_t v0, v1, v2;  // CCW when viewed from OUTSIDE
};

class SurfaceMesh {
    // === Construction ===
    uint32_t addVertex(pos, freedom, moveVector);
    uint32_t addVertex(x, y, z, freedom);
    uint32_t addTriangle(v0, v1, v2);

    // === Access ===
    const Vertex& vertex(i) const;
    Vertex& vertex(i);
    const Triangle& triangle(i) const;
    uint32_t vertexCount() const;
    uint32_t triangleCount() const;
    const std::vector<Vertex>& vertices() const;

    // === Freedom ===
    void setAllFreedom(VertexFreedom f);
    void setVertexFreedom(i, f, moveVec);

    // === DOF System (respects controlPointStride) ===
    uint32_t dofCount() const;         // free DOFs on this surface
    void applyParams(params, offset);  // write params → vertex positions
    void extractParams(target, offset); // read vertex positions → params
    void getBounds(lower, upper, offset);

    // === Control Points ===
    uint32_t setControlPointStride(stride);
    int controlPointStride() const;
    bool isControlPoint(vi) const;
    void downsampleVertexGradient(fullGrad, coarse, offset, axis);
    void interpolateFromControlPoints();

    // === Neighbors (for regularization) ===
    void buildNeighbors();
    const vector<uint32_t>& neighborVertices(vi) const;

    // === Padding ===
    void setInteriorGrid(interiorDim, paddingRings);
    int interiorGridDim() const;     // N (invertible interior)
    int paddingRings() const;        // P (padding rings beyond interior)
    int fullGridDim() const;         // N + 2*P
    bool isPaddingVertex(vi) const;
    bool isInteriorVertex(vi) const;
    void extrapolatePadding(upperZ, lowerZ);  // dip continuation + pinching

    // === Naming ===
    void setName(name);
    const string& name() const;
    bool isValid() const;
};
```

## Design Notes

- **controlPointStride=0** means all vertices are free DOFs (use direct fast-path)
- **Non-control vertices** are bilinearly interpolated from surrounding control points
- **Padding vertices** are extrapolated from interior boundary, optionally clamped by upper/lower Z arrays
- **Neighbors** are the 4-8 adjacent grid vertices (used by Laplacian smoothness)
- Grid is structured quads triangulated into 2 triangles per cell

## boundary_loop.h — Shared Boundary Utility

```cpp
#include <litho_invert/surface/boundary_loop.h>

struct BoundaryEdge { uint32_t start, end; };  // CCW-directed boundary edge
std::vector<BoundaryEdge> findBoundaryEdges(const SurfaceMesh& mesh);
std::vector<std::vector<uint32_t>> extractBoundaryLoops(const SurfaceMesh& mesh);
```

Used by GravityForward, MagneticForward (polyhedron side walls), and InversionExporter (closed-volume export).
Extracted from 6 duplicate copies across gravity_forward and magnetic_forward.

## Build

```powershell
cd modules/litho-surface
qmake litho-surface.pro
nmake release
```

## Dependencies

- litho-core (Eigen types)

## Tests

```powershell
cd modules/litho-surface/tests
qmake tests.pro && nmake release && release\tests.exe
```

Tests: DOF counting, param apply/extract, neighbor building, control-point downsampling,
padding extrapolation, interior/padding vertex classification, bounds.
