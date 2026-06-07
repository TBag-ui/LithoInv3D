# litho-model — Lithology Groups and Model Container

## Purpose

Manages N lithology groups with N-1 boundary surfaces. The central container that
all forward models and the inversion runner operate on.

**Critical invariant:** `groups.size() == surfaces.size() + 1` (always).

## Architecture

```
Group 0 (top)     ← Topography / z=0 (implicit)
  ↓
Surface 0         ← Boundary between groups 0 and 1
  ↓
Group 1
  ↓
Surface 1         ← Boundary between groups 1 and 2
  ↓
Group 2 (bottom)  ← Extends to bottomDepth
```

## API

```cpp
#include <litho_invert/litho/litho_group.h>
#include <litho_invert/litho/lithology_model.h>

struct LithoGroup {
    int id = -1;
    std::string name;
    double density = 0.0;           // g/cm³ — 0 triggers RUNTIME WARNING
    double susceptibility = 0.0;    // SI
    double conductivity = 1e-4;     // S/m
    // Remanence
    double remanence_magnitude = 0.0;
    double remanence_inclination = 0.0;
    double remanence_declination = 0.0;
    Vector3d magnetization;
};

class LithologyModel {
    // Group management
    int addGroup(const LithoGroup& group);
    void removeGroup(int groupId);
    const vector<LithoGroup>& groups() const;
    const LithoGroup& group(int index) const;
    LithoGroup& group(int index);       // non-const for property updates
    int groupCount() const;

    // Surface management
    void addSurface(shared_ptr<SurfaceMesh> surface);
    SurfaceMesh* surface(int index);
    const SurfaceMesh* surface(int index) const;
    int surfaceCount() const;
    const vector<shared_ptr<SurfaceMesh>>& surfaces() const;

    // DOF system (delegates to all surfaces)
    uint32_t totalDofCount() const;
    VectorXd assembleParameterVector() const;
    void applyParameterVector(const VectorXd& params);
    void getBounds(VectorXd& lower, VectorXd& upper) const;

    // Point-in-polyhedron classification
    int classifyPoint(const Vector3d& point) const;  // returns group index, -1 outside

    // Control-point downsampling
    uint32_t setControlPointStride(int stride);
    int controlPointStride() const;

    void setBottomDepth(double depth);
    double bottomDepth() const;
    bool isValid() const;  // asserts invariant
};
```

## Density Safety

When `LithoGroup::density == 0.0`, the forward models emit a runtime warning:
`[GRAVITY] Group N 'name' has density=0.0 — uninitialized?`

This catches the common bug of forgetting to set density after constructing groups.

## Build

```powershell
cd modules/litho-model
qmake litho-model.pro
nmake release
```

## Dependencies

- litho-core (Eigen types)
- litho-surface (SurfaceMesh)

## Tests

```powershell
cd modules/litho-model/tests
qmake tests.pro && nmake release && release\tests.exe
```

Tests: group creation, surface management, invariant enforcement, point classification,
DOF assembly/disassembly, bounds propagation, stride delegation.
