#pragma once
#include <litho_invert/core/common.h>
#include <vector>
#include <string>
#include <cstdint>

namespace litho_invert {

enum class VertexFreedom : uint8_t {
    FIXED = 0,      // No DOF contribution — vertex doesn't move
    Z_ONLY,         // 1 DOF: depth only
    X_ONLY,         // 1 DOF: x only
    Y_ONLY,         // 1 DOF: y only
    ALONG_VECTOR,   // 1 DOF: along user-specified unit vector
    XY_FREE,        // 2 DOF: free in horizontal plane (X+Y)
    XZ_FREE,        // 2 DOF: free in vertical N/S plane (X+Z, Y fixed)
    YZ_FREE,        // 2 DOF: free in vertical E/W plane (Y+Z, X fixed)
    XYZ_FREE        // 3 DOF: full 3D freedom
};

struct Vertex {
    Vector3d position;
    VertexFreedom freedom = VertexFreedom::FIXED;
    Vector3d moveVector = Vector3d::UnitZ();  // only used for ALONG_VECTOR

    Vertex() = default;
    Vertex(const Vector3d& pos, VertexFreedom f = VertexFreedom::FIXED,
           const Vector3d& mv = Vector3d::UnitZ())
        : position(pos), freedom(f), moveVector(mv) {}
    Vertex(double x, double y, double z, VertexFreedom f = VertexFreedom::FIXED)
        : position(x, y, z), freedom(f) {}
};

struct Triangle {
    uint32_t v0 = 0, v1 = 0, v2 = 0;
    Triangle() = default;
    Triangle(uint32_t a, uint32_t b, uint32_t c) : v0(a), v1(b), v2(c) {}
};

class SurfaceMesh {
public:
    SurfaceMesh() = default;

    // --- Vertex management ---
    uint32_t addVertex(const Vector3d& pos, VertexFreedom freedom = VertexFreedom::FIXED);
    uint32_t addVertex(const Vector3d& pos, VertexFreedom freedom, const Vector3d& moveVector);
    uint32_t addVertex(double x, double y, double z, VertexFreedom freedom = VertexFreedom::FIXED);

    const Vertex& vertex(uint32_t i) const;
    Vertex& vertex(uint32_t i);
    uint32_t vertexCount() const;
    const std::vector<Vertex>& vertices() const;

    // Set freedom for all vertices in one call (convenience)
    void setAllFreedom(VertexFreedom f);
    // Set freedom + move vector for a single vertex
    void setVertexFreedom(uint32_t i, VertexFreedom f, const Vector3d& moveVec = Vector3d::UnitZ());

    // --- Triangle management ---
    uint32_t addTriangle(uint32_t v0, uint32_t v1, uint32_t v2);
    const Triangle& triangle(uint32_t i) const;
    uint32_t triangleCount() const;
    const std::vector<Triangle>& triangles() const;

    // --- DOF / Parameter vector ---
    // Returns how many DOFs this surface contributes to inversion
    uint32_t dofCount() const;

    // Maps each DOF index to (vertex, axis). Axis: 0=X, 1=Y, 2=Z.
    struct DofMapping {
        uint32_t vertexIndex = 0;
        uint8_t axis = 0;  // 0=X, 1=Y, 2=Z

        DofMapping() = default;
        DofMapping(uint32_t vi, uint8_t ax) : vertexIndex(vi), axis(ax) {}
    };
    const std::vector<DofMapping>& dofMappings() const;

    // Apply parameter values from the full parameter vector to vertices.
    // 'params' is the full parameter vector; 'offset' is where this surface's DOFs start.
    void applyParams(const VectorXd& params, uint32_t offset);

    // Write current DOF values into a parameter vector at the given offset.
    void extractParams(VectorXd& target, uint32_t offset) const;

    // Get lower and upper bounds for each DOF (written into vectors at offset).
    void getBounds(VectorXd& lower, VectorXd& upper, uint32_t offset) const;

    // Set per-surface depth bounds (used for Z_ONLY and XYZ_FREE DOFs)
    void setBounds(double minZ, double maxZ);
    double minDepth() const;
    double maxDepth() const;

    // --- Parameter scaling ---
    // Each axis can have a scale factor for the parameter vector.
    // applyParams multiplies by scale, extractParams divides by scale.
    // This keeps the optimizer working in a well-conditioned space
    // where 1 unit ≈ similar physical significance across axes.
    void setAxisScale(int axis, double scale);
    double axisScale(int axis) const;

    // --- Control-point downsampling ---
    // When stride > 0, only every stride-th vertex (in grid index space)
    // gets a free DOF. Non-control-point vertices are bilinearly interpolated
    // from surrounding control points. This reduces the parameter count for
    // the gravity inversion while preserving the full mesh for other methods.
    // Returns the number of control points, or 0 if the mesh is not a regular grid.
    uint32_t setControlPointStride(int stride);
    // Set control points from externally-provided vertex indices.
    // Used by LithologyModel for mesh-independent CP selection across
    // group volume meshes that share contact vertices.
    uint32_t setControlPointsExplicit(const std::vector<uint32_t>& cpVertexIndices);
    int controlPointStride() const { return m_controlPointStride; }
    bool isControlPoint(uint32_t vertexIndex) const;
    uint32_t controlPointCount() const { return static_cast<uint32_t>(m_controlPointIndices.size()); }

    // Downsample a full-resolution per-vertex gradient (one scalar per vertex
    // for a given DOF axis) to the control-point DOF space using the chain rule
    // through bilinear interpolation weights.
    // axis: 0=X, 1=Y, 2=Z — filters DOF mappings to match the given axis.
    void downsampleVertexGradient(const std::vector<double>& fullVertexGrad,
                                  VectorXd& coarseGrad, uint32_t coarseOffset,
                                  int axis) const;

    // Interpolate non-control-point vertices from control point z-values.
    // Called automatically by applyParams() when stride is active.
    void interpolateFromControlPoints();

    // --- Neighbor relationships (for regularization) ---
    void buildNeighbors();
    const std::vector<uint32_t>& neighborVertices(uint32_t vertexIndex) const;

    // --- Mesh validation ---
    bool isValid() const;

    // --- Lateral padding ---
    // Set the interior (invertible) region for a regular grid with padding rings.
    // interiorDim = N (the N×N interior), paddingRings = P (extra rings each side).
    // Total grid = (N + 2P) × (N + 2P).
    void setInteriorGrid(int interiorDim, int paddingRings);
    int interiorGridDim() const { return m_interiorGridDim; }
    int paddingRings() const { return m_paddingRings; }
    int fullGridDim() const { return m_interiorGridDim + 2 * m_paddingRings; }
    bool isPaddingVertex(uint32_t vertexIndex) const;
    bool isInteriorVertex(uint32_t vertexIndex) const;

    // Extrapolate Z values from interior to padding region.
    // Each padding vertex continues the dip from the nearest interior boundary.
    // Optional upperZ/lowerZ (per-vertex, full grid) clamp for pinching.
    void extrapolatePadding(const std::vector<double>& upperZ = {},
                            const std::vector<double>& lowerZ = {});

    // --- Naming ---
    void setName(const std::string& name);
    const std::string& name() const;

    // --- 3D side classification ---
    // Returns +1 if the point is on the same side as the surface normal
    // (the "above" / group_above side for marching-cubes contact surfaces),
    // -1 if on the opposite side, 0 if exactly on the surface.
    // Uses per-triangle normals oriented via marching-cubes convention.
    int classifySide(const Vector3d& point) const;

private:
    std::vector<Vertex> m_vertices;
    std::vector<Triangle> m_triangles;
    std::vector<DofMapping> m_dofMappings;
    std::vector<std::vector<uint32_t>> m_neighbors;  // adjacency per vertex
    std::string m_name;
    double m_minDepth = -10000.0;  // positive down, so -10000 = 10km max depth
    double m_maxDepth = 0.0;       // default: surface at 0

    // Lateral padding
    int m_interiorGridDim = 0;
    int m_paddingRings = 0;

    // Control-point downsampling
    int m_controlPointStride = 0;  // 0 = disabled
    int m_gridDim = 0;             // sqrt(vertexCount) for regular grid
    std::vector<uint32_t> m_controlPointIndices;  // vertex indices of control points
    // For each vertex: list of (controlPointIndex, bilinearWeight) pairs
    std::vector<std::vector<std::pair<uint32_t, double>>> m_vertexToControlWeights;

    // Recompute m_dofMappings from current vertex freedoms
    void rebuildDofMappings();
    void ensureDofsBuilt() const;

    // Unstructured-mesh control-point downsampling via farthest-point sampling
    // and inverse-distance-weighted interpolation. Used when the mesh is not
    // an N×N regular grid.
    uint32_t setControlPointStrideUnstructured(int stride);

    mutable bool m_dofsDirty = false;

    // Parameter scaling per axis (default 1.0 = no scaling)
    double m_axisScales[3] = {1.0, 1.0, 1.0};

    // Cached per-triangle normals for 3D side classification.
    // Normals follow marching-cubes convention: point outward from the
    // binary volume (i.e., toward group_above for contact surfaces).
    mutable std::vector<Vector3d> m_faceNormals;
    mutable bool m_normalsComputed = false;
    void ensureFaceNormals() const;
};

} // namespace litho_invert
