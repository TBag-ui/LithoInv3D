# LithoInvert3D — Developer Guide

**Level 2 documentation** — architecture, build system, APIs, and extension points.

## Architecture

```
InversionRunner (modules/litho-inversion)
  ├── JointObjective
  │   ├── ObjectiveFunction (gravity core + reg + constraints)
  │   ├── MagneticForward (optional)
  │   ├── EMActiveForward (optional)
  │   └── EMMTForward (optional)
  ├── LBFGSBOptimizer (bound-constrained L-BFGS)
  ├── LithologyModel (modules/litho-model)
  │   └── N × SurfaceMesh (modules/litho-surface)
  └── Forward models (modules/litho-forward, modules/litho-em)
      └── Geometry primitives (modules/litho-core)
```

### Data Flow

```
INI config → InversionConfig → InversionRunner::setup()
  ├── Creates forward models from config
  ├── Wires up objective + regularization + constraints
  └── InversionRunner::run()
      ├── Phase A: Geometry optimization (L-BFGS-B on surface DOFs)
      ├── Phase B: Property inversion (least-squares on densities/susc)
      └── Export via InversionExporter (modules/litho-io)
```

## Build System

### C++ Static Library + Executable

Each module builds as a static library. The top-level executable links all modules.

**Build order** (dependency chain):
1. `litho-core` — no deps
2. `litho-surface` — depends on core
3. `litho-model` — depends on core, surface
4. `litho-forward` — depends on core, surface, model
5. `litho-em` — depends on core, surface, model
6. `litho-io` — depends on core, surface, model
7. `litho-regularization` — depends on core, surface, model
8. `litho-inversion` — depends on all above

Each module's `.pro` file:
```qmake
TEMPLATE = lib
CONFIG += staticlib c++17
INCLUDEPATH += ../../vendor/eigen
# plus module-specific HEADERS/SOURCES
```

### Building

```powershell
# From project root:
.\build\build.bat   # calls qmake + nmake for all modules, then links executable

# Build a single module:
cd modules\litho-core
qmake litho-core.pro
nmake release
```

### Python

```bash
pip install -e python/cluster_api
pip install -e python/lithoseed
```

## Module APIs

### litho-core
- `Vector3d`, `VectorXd`, `MatrixXd` — Eigen typedefs (`common.h`)
- `GravityPoint`, `MagneticPoint`, `Constraint` — data structs (`common.h`)
- `TopographyConfig`, `TopographyMode` — topography settings (`common.h`)
- `solidAngle()`, `lineIntegralTerm()`, `surfaceIntegralTerm()` — geometry (`geometry.h`)
- `BoundaryLoop`, `findBoundaryEdges()` — shared boundary utility (`boundary_loop.h`)
- `computeDurbinWatson()` — spatial autocorrelation (`stats.h`)

### litho-surface
- `Vertex`, `Triangle`, `VertexFreedom` — surface primitives
- `SurfaceMesh` — full mesh with DOF system, padding, control points
  - `addVertex()`, `addTriangle()`, `buildNeighbors()`
  - `dofCount()`, `applyParams()`, `extractParams()`, `getBounds()`
  - `setControlPointStride()`, `interpolateFromControlPoints()`
  - `extrapolatePadding()`, `isInteriorVertex()`, `isPaddingVertex()`

### litho-model
- `LithoGroup` — density, susceptibility, conductivity, remanence
- `LithologyModel` — N groups + N-1 surfaces
  - Invariant: `groups.size() == surfaces.size() + 1`
  - `totalDofCount()`, `assembleParameterVector()`, `applyParameterVector()`
  - `classifyPoint()` — ray-cast point-in-polyhedron test
  - `setControlPointStride()` — delegates to all surfaces

### litho-forward
- `ForwardModel` — abstract base: `compute()`, `computeJacobian()`, `parameterCount()`
- `GravityForward` — Nagy polyhedron formula
  - `compute()`, `computeJacobian()`, `computeGroupUnitResponse()`
  - `gravityClosedMesh()`, `gravityFacet()` (static helpers)
- `MagneticForward` — Okabe/Plouff polyhedron formula
  - `compute()`, `computeJacobian()`, `computeGroupUnitResponse()`
  - `magneticClosedMesh()`, `magneticFacet()` (static helpers)

### litho-em
- `EMActiveForward` — airborne/large-loop TEM (fallback dipole approx)
- `EMMTForward` — magnetotelluric (fallback 1D half-space)
- `EMSolver` — pluggable solver interface (IE solver stubbed)
- `EMConfig`, `EMSource`, `EMReceiver`, etc. — data types

### litho-inversion
- `Optimizer` → `LBFGSBOptimizer` — bound-constrained L-BFGS
- `ObjectiveFunction` — gravity misfit + regularization + constraints
- `JointObjective` — weighted multi-physics objective
- `ConstraintHandler` — drillhole penalty
- `InversionRunner` + `InversionConfig` + `InversionResult`

### litho-io
- Loaders: `CSVGravityLoader`, `OBJSurfaceLoader`, `CSVConstraintLoader`, `loadDEM()`
- `IniConfig` — simple INI parser
- `InversionExporter` — OBJ, TS (GOCAD TSurf), UBC-GIF, CSV, closed volumes

### litho-regularization
- `Regularization` — abstract base
- `SurfaceSmoothness` — Laplacian curvature penalty
- `ReferenceModelRegularization` — deviation from starting model

## Extension Points

| Feature | Where to Hook |
|---------|---------------|
| New forward model | Extend `ForwardModel` in `litho-forward/` |
| New EM solver | Implement `EMSolver` in `litho-em/` |
| New optimizer | Implement `Optimizer` in `litho-inversion/` |
| New regularizer | Extend `Regularization` in `litho-regularization/` |
| New data type | Add struct + forward + wire into `JointObjective` |
| New export format | Add method to `InversionExporter` in `litho-io/` |

## Coordinate Conventions

- **Vertices**: Z positive UP (surface at z=0, deep at z=-5000)
- **Constraint depths**: POSITIVE DOWN (z_top=100 means top at 100m depth)
- **DEM elevations**: positive UP (matching vertices)
- **UBC-GIF export**: Z positive DOWN (converted at export time)
- **Gravity**: positive DOWN (normal gravity attraction)

## Known Limitations

- **IE solver stubbed**: EM forward models use fallback approximations until full IE solver is implemented
- **No parallelization**: Forward model evaluations are single-threaded
- **Padding group must be last**: The deep half-space padding group must be group N-1
