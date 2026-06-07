#pragma once
#include <litho_invert/litho/litho_group.h>
#include <litho_invert/surface/surface_mesh.h>
#include <litho_invert/core/common.h>
#include <vector>
#include <memory>

namespace litho_invert {

class GroupProposal;

struct ValidationResult {
    bool passed = true;
    std::string failureReason;
    double worstViolation = 0.0;  // magnitude of worst violation (metres or count)
};

class LithologyModel {
public:
    LithologyModel() = default;

    // --- Group management ---
    int addGroup(const LithoGroup& group);
    void removeGroup(int groupId);

    const std::vector<LithoGroup>& groups() const;
    const LithoGroup& group(int index) const;
    int groupCount() const;
    LithoGroup& group(int index);

    // --- Group mesh management ---
    // Each group has one closed triangulated boundary mesh.
    // For N groups you need N meshes. Groups can be adjacent in any direction.
    void setGroupMesh(int groupIndex, std::shared_ptr<SurfaceMesh> mesh);
    SurfaceMesh* groupMesh(int index);
    const SurfaceMesh* groupMesh(int index) const;
    int groupMeshCount() const;

    // --- Label grid for fast point classification ---
    // The pipeline exports a 3D integer array where each voxel holds
    // its group ID. This provides O(1) classification without ray casting.
    void setLabelGrid(const std::vector<int>& labels, int nx, int ny, int nz,
                      double x0, double y0, double z0, double dx, double dy, double dz);
    bool hasLabelGrid() const;

    // --- DOF / Parameter vector ---
    // Global DOF management with deduplication: vertices at the same
    // spatial position on the same axis share a single DOF, keeping
    // adjacent group meshes in contact during inversion.
    uint32_t totalDofCount() const;
    VectorXd assembleParameterVector() const;
    void applyParameterVector(const VectorXd& params);
    void getBounds(VectorXd& lower, VectorXd& upper) const;

    // Look up the global DOF index for a (meshIndex, localDofIndex) pair.
    // Returns -1 if not found. Used by regularization to map per-mesh
    // DOF contributions to the global gradient vector.
    int globalDofIndex(int meshIndex, int localDofIndex) const;

    // --- Point classification ---
    // Returns which litho group is at a given 3D point.
    // Primary path: label grid O(1) lookup. Fallback: ray casting.
    // Returns -1 if point is outside all groups or below bottom depth.
    int classifyPoint(const Vector3d& point) const;

    // --- Control-point downsampling ---
    uint32_t setControlPointStride(int stride);
    int controlPointStride() const { return m_controlPointStride; }

    // --- Exterior face constraints ---
    // Detect faces appearing in only one group mesh (model boundary), classify
    // by dominant normal direction, and assign per-vertex freedom so boundary
    // vertices slide in-plane but cannot leave the model envelope.
    void fixExteriorFaces();

    // --- Model bounds enforcement ---
    // Compute AABB of all group mesh vertices, then set per-mesh depth bounds
    // and tighten global DOF bounds to the model envelope plus a margin.
    // marginMultiplier scales the margin (1.0 = 1% of max dimension).
    void applyModelBounds(double marginMultiplier = 1.0);

    // --- Depth extent ---
    void setBottomDepth(double depth);
    double bottomDepth() const;

    // --- Shared-vertex synchronization ---
    // Builds m_sharedVertexSync from all mesh vertex positions, identifying
    // vertices shared between ≥2 meshes. Called by fixExteriorFaces() and
    // setControlPointStride(). Ensures applyParameterVector keeps contact
    // surfaces gap-free even without control-point downsampling.
    void buildSharedVertexSync();

    // --- Validation ---
    bool isValid() const;

    // Quality validation: checks shared vertex gaps, adjacent/non-adjacent
    // overlaps, edge tearing, AABB containment, spike detection, and face
    // orientation. Returns pass/fail with the worst violation magnitude.
    // overlapBaseline is the starting non-adjacent overlap count (pre-inversion).
    // adjacentOverlapBaseline is the starting adjacent overlap count — both
    // are subtracted so discretization noise doesn't cause false positives.
    ValidationResult validate(size_t overlapBaseline = 0,
                              size_t adjacentOverlapBaseline = 0) const;

    // --- Group proposal hook ---
    void setGroupProposal(std::shared_ptr<GroupProposal> proposal);

private:
    std::vector<LithoGroup> m_groups;
    std::vector<std::shared_ptr<SurfaceMesh>> m_groupMeshes;

    // Global DOF deduplication — contact vertices shared between adjacent groups
    struct GlobalDofRef {
        uint32_t meshIndex = 0;
        uint32_t localDofIndex = 0;
    };
    // m_globalDofRefs[g] = all mesh-local DOFs controlled by global DOF g.
    // Typically 1 entry for interior vertices, 2+ for shared contact vertices.
    mutable std::vector<std::vector<GlobalDofRef>> m_globalDofRefs;
    mutable bool m_globalDofsDirty = true;
    void rebuildGlobalDofs() const;

    // Label grid for fast point classification
    std::vector<int> m_labelGrid;
    int m_labelNx = 0, m_labelNy = 0, m_labelNz = 0;
    double m_labelX0 = 0, m_labelY0 = 0, m_labelZ0 = 0;
    double m_labelDx = 0, m_labelDy = 0, m_labelDz = 0;

    // Shared-vertex synchronization: for every shared spatial position
    // (vertex appearing in ≥2 group meshes), we store the list of
    // (meshIndex, vertexIndex) so that applyParameterVector can enforce
    // identical positions after interpolation.  This eliminates voids
    // at non-CP shared vertices that would otherwise interpolate
    // differently in each mesh.
    struct SharedVertexRef {
        uint32_t meshIndex = 0;
        uint32_t vertexIndex = 0;
    };
    // Each entry groups all mesh/vertex instances of one shared spatial position.
    // Built by setControlPointStride, enforced by applyParameterVector.
    std::vector<std::vector<SharedVertexRef>> m_sharedVertexSync;

    double m_bottomDepth = -5000.0;
    std::shared_ptr<GroupProposal> m_groupProposal;
    int m_controlPointStride = 0;

    // Cached model AABB for applyModelBounds
    Vector3d m_modelBoundsMin = Vector3d(-1e12, -1e12, -1e12);
    Vector3d m_modelBoundsMax = Vector3d(1e12, 1e12, 1e12);
    double m_modelBoundsMargin = 0.0;
};

} // namespace litho_invert

