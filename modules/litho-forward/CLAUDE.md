# litho-forward — Analytical Gravity and Magnetic Forward Models

## Purpose

Computes predicted geophysical responses from a LithologyModel. Uses analytical
polyhedron formulas: Nagy (2000) for gravity, Okabe/Plouff for magnetics.

Each litho group is a closed polyhedron: top surface (CCW) + bottom surface (reversed)
+ side walls stitched from boundary edges. Reuses `findBoundaryEdges()` from litho-core.

## API

```cpp
#include <litho_invert/forward/forward_model.h>
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>

// === Abstract base ===
class ForwardModel {
public:
    virtual VectorXd compute(const VectorXd& params) = 0;
    virtual MatrixXd computeJacobian(const VectorXd& params);  // default: FD
    virtual void computeBoth(const VectorXd& params, VectorXd& predicted, MatrixXd& jacobian);
    virtual size_t dataCount() const = 0;
    virtual size_t parameterCount() const = 0;  // MUST return m_model->totalDofCount()
};

// === GravityForward (Nagy 2000) ===
GravityForward(shared_ptr<LithologyModel>, const GravityData&);

VectorXd compute(params);                           // → nData gravity (mGal)
MatrixXd computeJacobian(params);                   // analytical Z, FD XY
VectorXd computeGroupUnitResponse(groupIndex);      // per unit density

// Static helpers
static double gravityClosedMesh(obs, closedMesh, density_g_per_cm3);
static double gravityFacet(obs, a, b, c);           // dimensionless geometric factor

// === MagneticForward (Okabe/Plouff) ===
MagneticForward(shared_ptr<LithologyModel>, const MagneticData&, inc, dec, field_nT);

VectorXd compute(params);                           // → nData TMI (nT)
MatrixXd computeJacobian(params);
VectorXd computeGroupUnitResponse(groupIndex);      // per unit susceptibility
VectorXd computeGroupUnitResponseDirection(groupIndex, magDir);
VectorXd computeGroupUnitRemanenceResponse(groupIndex);

// Static helpers
static double magneticClosedMesh(obs, closedMesh, chi, fieldUnitVector, field_nT);
static double magneticFacet(obs, a, b, c, magDir, fieldUnitVector);
```

## Constants

```cpp
static constexpr double G_SI = 6.67430e-11;
static constexpr double M_S2_TO_MGAL = 1e5;
static constexpr double DENSITY_SCALE = 1000.0;  // g/cm³ → kg/m³
```

## Design Notes

- **parameterCount()** MUST return `m_model->totalDofCount()` — never `vertexCount()`
- **Density safety**: emits `[GRAVITY]` warning when group density is 0.0
- **Polyhedron construction** uses `findBoundaryEdges()` from boundary_loop.h for side walls
- **FD Jacobian step**: h=0.1 for geometry DOFs, h=1e-4 for facet-level gradient

## Build

```powershell
cd modules/litho-forward
qmake litho-forward.pro
nmake release
```

## Dependencies

- litho-core (geometry, boundary_loop, types)
- litho-surface (SurfaceMesh)
- litho-model (LithologyModel)

## Tests

```powershell
cd modules/litho-forward/tests
qmake tests.pro && nmake release && release\tests.exe
```

Tests: Nagy vs analytical sphere, Nagy vs infinite slab, unit response scaling,
padding response, Bouguer correction, FD Jacobian validation.
