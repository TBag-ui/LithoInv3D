# LithoInvert3D — Complete LLM Context Document

**Level 4 documentation** — designed to fully re-instantiate an AI coding session without
needing to search or read any other file. Load this first, then start working.

**Date:** 2026-06-02
**Language:** C++17 + Python 3.9+
**Build:** qmake + nmake (MSVC 2022) or MinGW-w64 11.2
**Linear algebra:** Eigen 3.4.0 (vendored in `vendor/eigen/`)

---

## 1. What This Project Is

LithoInvert3D recovers 3D subsurface lithological boundary surfaces from gravity,
magnetic, and EM data. Instead of inverting for voxel property values, it inverts
for the positions of explicit triangulated surfaces that separate lithological units.

**Key features:**
- Analytical forward models: Nagy polyhedron (gravity), Okabe/Plouff (magnetics)
- L-BFGS-B optimization with bound constraints
- Control-point DOF downsampling for performance
- Borehole constraints, topography support, lateral padding
- XYZ_FREE vertex freedom (full 3-DOF surface inversion)
- Reference model regularization, hard depth bounds
- Python starting-model pipeline: SimPEG → GMM clustering → voxel-face volume extraction + contact surfaces

---

## 2. Project Structure

```
LithoInvert3D/
├── README.md                     # Level 1: Overview
├── CLAUDE.md                     # Master Claude Code instructions (THIS FILE's summary)
├── DUPLICATE_TRACKING.txt        # Records of duplicate code across original codebase
│
├── docs/
│   ├── developer.md              # Level 2: Architecture + API reference
│   ├── scientific_reference.md   # Level 3: Mathematical formulation
│   └── llm_context.md            # Level 4: THIS FILE
│
├── modules/                      # C++ submodules (each independently buildable)
│   ├── litho-core/               # common.h, geometry.h, stats.h
│   ├── litho-surface/            # SurfaceMesh, Vertex, Triangle, VertexFreedom
│   ├── litho-model/              # LithoGroup, LithologyModel
│   ├── litho-forward/            # GravityForward, MagneticForward (analytical)
│   ├── litho-em/                 # EMActiveForward, EMMTForward (stubs)
│   ├── litho-inversion/          # LBFGSB, Objective, JointObjective, Runner
│   ├── litho-io/                 # Loaders, InversionExporter, IniConfig
│   └── litho-regularization/     # Smoothness, ReferenceModel
│
├── python/                       # Python packages
│   ├── lithoseed/                # SimPEG → clustering → surfaces → INI
│   └── cluster_api/              # GMM clustering + mesh I/O (NO duplicates)
│
├── apps/                         # End-to-end inversion applications
│   ├── voiseys_bay/              # Synthetic Ni-Cu-Co test case
│   ├── forrestania/              # Regional gravity-magnetic inversion
│   └── quest_south/              # Regional application
│
├── datasets/                     # Test data, INI configs, synthetic models
├── scripts/                      # batch_deposits.py, build helpers
├── vendor/eigen/                 # Eigen 3.4.0 (header-only)
└── build/                        # Build output + build.bat
```

### Module Dependency Order

```
litho-core  (NO deps)
  ↓
litho-surface  (depends: core)
  ↓
litho-model  (depends: core, surface)
  ↓
litho-forward  litho-em  litho-io  litho-regularization  (depend: core, surface, model)
  ↓
litho-inversion  (depends: ALL above)
```

Each module has its own `CLAUDE.md` in `modules/<module-name>/CLAUDE.md` with
API reference, build instructions, and test locations.

---

## 3. Complete C++ Data Types

### 3.1 Core Types (`litho-core/include/litho_invert/core/common.h`)

```cpp
namespace litho_invert {

using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using Index    = Eigen::Index;  // ptrdiff_t

struct GravityPoint {
    Vector3d position;      // (x, y, z) — z positive UP
    double g_obs = 0.0;     // observed gravity (mGal)
    double g_std = 0.0;     // std dev (0 = unweighted, w=1.0)
};
using GravityData = std::vector<GravityPoint>;

struct MagneticPoint {
    Vector3d position;       // (x, y, z)
    double t_obs = 0.0;      // total-field anomaly (nT)
    double t_std = 0.0;
};
using MagneticData = std::vector<MagneticPoint>;

struct Constraint {
    Vector3d position;           // (x, y) center
    double z_top = 0.0;          // top depth — POSITIVE DOWN
    double z_bottom = 0.0;       // bottom depth — POSITIVE DOWN
    int litho_group_id = -1;     // known lithology in this interval
};

enum class TopographyMode { None, Raw, TerrainCorrected };

struct TopographyConfig {
    TopographyMode mode = TopographyMode::None;
    std::string demFile;
    double datumElevation = 0.0;
    double bouguerDensity = 2.67;  // g/cm³
    int paddingRings = 0;
    double paddingCellSize = 0.0;
    bool invertHalfspaceProperties = false;

    // Controls how closure surfaces (top, bottom, deep) are built
    // when padding rings extend beyond the survey area:
    //   IndependentFlat       — independent horizontal grid from bounding box (default)
    //   ProjectedBeyondSurvey — copy x,y topology from reference surface
    PaddingProjection paddingProjection = PaddingProjection::IndependentFlat;
};
```

**CRITICAL COORDINATE CONVENTIONS:**
- Vertices: z positive UP (surface at z=0, deep at z=-5000)
- Constraint depths: POSITIVE DOWN (z_top=100 means 100m below surface)
- UBC-GIF export: z positive DOWN (converted at export time)

### 3.2 Boundary Loop Utility (`litho-surface/include/litho_invert/surface/boundary_loop.h`)

Shared utility extracted from 6 duplicate copies in gravity_forward and magnetic_forward.
Lives in litho-surface (not litho-core) because it depends on `SurfaceMesh`.

```cpp
struct BoundaryEdge {
    uint32_t start;   // CCW start vertex
    uint32_t end;     // CCW end vertex
};

// Find all boundary edges of a surface mesh (edges appearing in exactly one triangle).
// Returns edges with CCW winding (following the surface perimeter direction).
std::vector<BoundaryEdge> findBoundaryEdges(const SurfaceMesh& mesh);

// Walk boundary edges into ordered loops (connected components of the boundary graph).
std::vector<std::vector<uint32_t>> extractBoundaryLoops(const SurfaceMesh& mesh);
```

Used by: GravityForward, MagneticForward, and InversionExporter (closed-volume export).

### 3.3 Surface Types (`litho-surface/include/litho_invert/surface/surface_mesh.h`)

```cpp
enum class VertexFreedom : uint8_t {
    FIXED = 0,      // 0 DOF
    Z_ONLY,         // 1 DOF — vertical only
    X_ONLY, Y_ONLY, // 1 DOF
    ALONG_VECTOR,   // 1 DOF — along Vertex::moveVector
    XY_FREE,        // 2 DOF — horizontal plane
    XYZ_FREE        // 3 DOF — full 3D
};

struct Vertex {
    Vector3d position;
    VertexFreedom freedom = VertexFreedom::FIXED;  // DEFAULT IS FIXED
    Vector3d moveVector = Vector3d::UnitZ();       // for ALONG_VECTOR
};

struct Triangle {
    uint32_t v0, v1, v2;  // CCW winding when viewed from OUTSIDE
};

class SurfaceMesh {
public:
    // Vertices
    uint32_t addVertex(const Vector3d& pos, VertexFreedom f = VertexFreedom::FIXED);
    uint32_t addVertex(double x, double y, double z, VertexFreedom f = VertexFreedom::FIXED);
    const Vertex& vertex(uint32_t i) const;
    Vertex& vertex(uint32_t i);
    uint32_t vertexCount() const;

    // Triangles
    uint32_t addTriangle(uint32_t v0, uint32_t v1, uint32_t v2);
    uint32_t triangleCount() const;

    // DOF system (respects controlPointStride)
    uint32_t dofCount() const;
    void applyParams(const VectorXd& params, uint32_t offset);
    void extractParams(VectorXd& target, uint32_t offset) const;
    void getBounds(VectorXd& lower, VectorXd& upper, uint32_t offset) const;
    void setBounds(double minZ, double maxZ);

    // Control-point downsampling
    uint32_t setControlPointStride(int stride);
    int controlPointStride() const;
    bool isControlPoint(uint32_t vi) const;
    void downsampleVertexGradient(const std::vector<double>& fullGrad,
                                  VectorXd& coarse, uint32_t offset, int axis) const;
    void interpolateFromControlPoints();

    // Neighbors (for regularization)
    void buildNeighbors();
    const std::vector<uint32_t>& neighborVertices(uint32_t vi) const;

    // Padding
    void setInteriorGrid(int interiorDim, int paddingRings);
    int interiorGridDim() const;
    int paddingRings() const;
    int fullGridDim() const;
    bool isPaddingVertex(uint32_t vi) const;
    bool isInteriorVertex(uint32_t vi) const;
    void extrapolatePadding(const std::vector<double>& upperZ = {},
                            const std::vector<double>& lowerZ = {});
    // ... plus name, isValid, etc.
};
```

### 3.4 Lithology Types (`litho-model/`)

```cpp
struct LithoGroup {
    int id = -1;
    std::string name;
    double density = 0.0;           // g/cm³ (0 = UNINITIALIZED — triggers runtime warning)
    double susceptibility = 0.0;    // SI
    double conductivity = 1e-4;     // S/m
    double remanence_magnitude = 0.0;     // A/m
    double remanence_inclination = 0.0;  // deg
    double remanence_declination = 0.0;  // deg
    Vector3d magnetization;               // A/m total remanence vector
};

class LithologyModel {
public:
    int addGroup(const LithoGroup& group);
    void addSurface(std::shared_ptr<SurfaceMesh> surface);
    const std::vector<LithoGroup>& groups() const;
    const LithoGroup& group(int index) const;
    LithoGroup& group(int index);
    int groupCount() const;
    SurfaceMesh* surface(int index);
    int surfaceCount() const;

    // DOF system — delegates to all surfaces
    uint32_t totalDofCount() const;
    VectorXd assembleParameterVector() const;
    void applyParameterVector(const VectorXd& params);
    void getBounds(VectorXd& lower, VectorXd& upper) const;

    // Point classification
    int classifyPoint(const Vector3d& point) const;

    // Control-point downsampling
    uint32_t setControlPointStride(int stride);
    int controlPointStride() const;

    bool isValid() const;  // groups.size() == surfaces.size() + 1
};

// INVARIANT: groups.size() == surfaces.size() + 1 (always)
// Group 0 → Surface 0 → Group 1 → Surface 1 → Group 2 (bottom)
```

---

## 4. Forward Models — Complete API

### 4.1 ForwardModel Base (`litho-forward/`)

```cpp
class ForwardModel {
public:
    virtual ~ForwardModel() = default;
    virtual VectorXd compute(const VectorXd& params) = 0;
    virtual size_t dataCount() const = 0;
    virtual size_t parameterCount() const = 0;  // MUST return m_model->totalDofCount()
    virtual MatrixXd computeJacobian(const VectorXd& params);  // default: FD
    virtual void computeBoth(const VectorXd& params, VectorXd& predicted, MatrixXd& jacobian);
};
```

### 4.2 GravityForward

Constructor: `GravityForward(shared_ptr<LithologyModel>, const GravityData&)`

Key methods:
- `VectorXd compute(const VectorXd& params)` → nData gravity (mGal)
- `MatrixXd computeJacobian(const VectorXd& params)` → analytical for Z, FD for XY
- `VectorXd computeGroupUnitResponse(int groupIndex)` → per unit density
- `static double gravityClosedMesh(obs, closedMesh, density)` → Nagy polyhedron gravity for a closed mesh
- `static double gravityFacet(obs, a, b, c)` → dimensionless geometric factor for one facet

Constants:
```cpp
static constexpr double G_SI = 6.67430e-11;
static constexpr double M_S2_TO_MGAL = 1e5;
static constexpr double DENSITY_SCALE = 1000.0;  // g/cm³ → kg/m³
```

Polyhedron construction reuses `findBoundaryEdges()` from `boundary_loop.h` for side walls.
Padding, topography, and Bouguer correction are handled at the `InversionRunner`/`ObjectiveFunction`
level, not on `GravityForward` directly.

### 4.3 MagneticForward

Constructor: `MagneticForward(shared_ptr<LithologyModel>, const MagneticData&, double inc, double dec, double field_nT)`

Same interface pattern as GravityForward. Uses Okabe/Plouff polyhedron formula.

### 4.4 EM Forward (stubs in `litho-em/`)

- `EMActiveForward` — airborne/large-loop TEM, trust-region subsetting, dipole fallback
- `EMMTForward` — MT, 1D half-space fallback
- `EMSolver` — pluggable solver interface (factory returns nullptr)
- `EMConfig` — solver method, time gates, subsetting params, physics flags

---

## 5. Inversion Framework

### 5.1 ObjectiveFunction (`litho-inversion/`)

```cpp
class ObjectiveFunction {
public:
    ObjectiveFunction(shared_ptr<ForwardModel>, const GravityData&);
    void addRegularization(shared_ptr<Regularization>);
    void setConstraintHandler(shared_ptr<ConstraintHandler>);
    double evaluate(const VectorXd& params);
    VectorXd gradient(const VectorXd& params);
    VectorXd residuals(const VectorXd& params);
    Components { double dataMisfit, regularization, constraintPenalty, total; }
    Components evaluateComponents(const VectorXd& params);
};
```

### 5.2 JointObjective

Weighted registry of sub-objectives:
```cpp
JointObjective(gravityObj, magForward, magData, alpha_mag);
void addActiveEM(activeEMForward, data, alpha);
void addMT(mtForward, data, alpha);

double evaluate(const VectorXd& params);
VectorXd gradient(const VectorXd& params);
JointComponents evaluateComponents(const VectorXd& params);
// Individual misfits: gravityMisfit(), magneticMisfit(), activeEMMisfit(), mtMisfit()
```

### 5.3 LBFGSBOptimizer

```cpp
class LBFGSBOptimizer : public Optimizer {
    void setHistorySize(int m);      // default 10
    void setClearOnMinimize(bool);   // false = warm-start
    OptimizerResult minimize(obj, grad, x0, lower, upper) override;
};

struct OptimizerResult {
    VectorXd params; int iterations; bool converged; double value;
};
```

### 5.4 InversionRunner

Main loop:
```
While totalIters < maxIterations AND not converged:
├─ Phase A: Geometry Optimization (L-BFGS-B on surface DOFs)
│   ├─ geoObj → JointObjective::evaluateComponents + DW computation
│   ├─ geoGrad → JointObjective::gradient
│   └─ Apply params → extrapolate padding
├─ Phase B: Property Inversion (if enabled)
│   ├─ Compute unit-response matrices U
│   └─ Bounded least-squares for densities/susc/conductivity
└─ Repeat
```

### 5.5 InversionConfig (key fields)

```cpp
struct InversionConfig {
    shared_ptr<LithologyModel> model;
    GravityData observedData;
    vector<Constraint> constraints;
    double lambda = 1.0;          // regularization weight
    double omega = 1e6;           // constraint penalty weight
    int maxIterations = 500;
    int controlPointStride = 0;   // 0 = all free
    bool enablePropertyInversion = false;
    int propertyInversionInterval = 50;
    // Magnetic (empty = skip)
    MagneticData magneticData; double magneticWeight = 1.0;
    double magneticInclination = 75.0, magneticDeclination = -20.0, magneticField_nT = 55000.0;
    // EM (empty = skip)
    ActiveEMData activeEMData; EMConfig emConfig; double activeEMWeight = 1.0;
    MTData mtData; vector<MTStation> mtStations; double mtWeight = 1.0;
    // Padding, topography, reference model, depth bounds, remanence...
    TopographyConfig topography;
};
```

### 5.6 InversionResult

```cpp
struct InversionResult {
    shared_ptr<LithologyModel> finalModel;
    vector<InversionIteration> history;
    VectorXd predictedData;
    bool converged; int totalIterations; double finalMisfit, finalRMS;
    VectorXd finalDensities, finalSusceptibilities, finalConductivities;
    double finalPaddingDensity, finalPaddingSusceptibility, finalPaddingConductivity;
    // Closure surfaces for export
    shared_ptr<SurfaceMesh> closureTop, closureBottom, closureDeepBottom;
};
```

---

## 6. I/O System (`litho-io/`)

### InversionExporter

```cpp
class InversionExporter {
public:
    InversionExporter(const string& outputDir, const string& baseName);
    void setSubfolder(const string& sub);
    void setGroupNaming(const vector<string>& groupExportNames);

    // Individual exports
    void exportTS(const SurfaceMesh& mesh, const string& suffix);       // GOCAD TSurf
    void exportOBJ(const SurfaceMesh& mesh, const string& suffix);
    void exportRawVertices(const SurfaceMesh&, const string&);           // .txt
    void exportRawTriangles(const SurfaceMesh&, const string&);

    // Closed volumes (per litho group)
    void exportClosedVolume(const SurfaceMesh* top, const SurfaceMesh* bot,
                            const string& suffix, double bottomDepth);

    // Starting model export (all contacts + closed volumes)
    void exportStartingModel(const LithologyModel& model,
                             const SurfaceMesh* flatTop, const SurfaceMesh* flatBottom,
                             double bottomDepth);

    // Interior-only (survey-area, truncates padding)
    void exportInteriorTS(const SurfaceMesh&, const string&);
    void exportInteriorOBJ(const SurfaceMesh&, const string&);

    // Bulk exports
    void exportAll(const InversionResult& result, const GravityData& observed,
                   double xmin, double xmax, double ymin, double ymax,
                   double zmin, double zmax, double cellSize);

    void exportUBCGIF(const LithologyModel&, ...);   // UBC-GIF mesh
    void exportLithoCSV(const LithologyModel&, ...); // per-cell CSV
    void exportPredictedCSV(const GravityData&, const VectorXd& predicted);
    void exportLog(const InversionResult&);
};
```

**Closed volume TS format**: Top surface (CCW), bottom surface (reversed), side walls
stitched via `extractBoundaryLoops()` from `boundary_loop.h`.

**Dual export (padding mode)**: When padding is active, `exportAll()` produces:
- `csv_inv_<name>.ts` — interior-only (survey area)
- `csv_inv_<name>_full.ts` — full model including padding (debug)

### Loaders

- `CSVGravityLoader`: CSV `x,y,z,g_obs,g_std` → `GravityData`
- `OBJSurfaceLoader`: Wavefront OBJ → `shared_ptr<SurfaceMesh>`
- `CSVConstraintLoader`: CSV `x,y,z_top,z_bottom,litho_id` → `vector<Constraint>`
- `JSONLithoConfigLoader`: JSON → `LithologyModel`
- `loadDEM(path, xVec, yVec)`: XYZ grid → `vector<double>` (bilinear interpolation)

### IniConfig

Simple INI parser (no external dependencies):
```cpp
class IniConfig {
    bool load(const string& path);
    string getString(section, key, default);
    double getDouble(section, key, default);
    int getInt(section, key, default);
    bool getBool(section, key, default);
};
```

Config wiring: `InversionConfig::fromIni(const IniConfig&)` in `litho-inversion`.

---

## 7. Regularization (`litho-regularization/`)

```cpp
class SurfaceSmoothness : public Regularization {
    explicit SurfaceSmoothness(shared_ptr<LithologyModel>);
    // Laplacian curvature: ½ Σ ||v − mean(neighbors(v))||²
};

class ReferenceModelRegularization : public Regularization {
    explicit ReferenceModelRegularization(shared_ptr<LithologyModel>);
    void captureReference(const VectorXd& params);  // snapshot starting state
    // Penalty: ½·λ_ref·||m − m_start||²
};
```

---

## 8. Python Packages

### 8.1 cluster_api (`python/cluster_api/`)

```python
# _cluster.py
from cluster_api import cluster_lithology
labels, summary = cluster_lithology(density, susceptibility, n_clusters=4)
# summary: LithologySummary with density_mean, density_std, susc_mean, susc_std, counts

# _io.py
from cluster_api import load_inversion, MeshData
mesh, density, susceptibility = load_inversion("model.msh", "density.mod", "susc.mod")
# mesh.to_3d_grid(values, fill=-1) → dense [nx, ny, nz] numpy array

# _io.py — MeshData dataclass, load_inversion(), to_3d_grid()
# _cluster.py — cluster_lithology(), cluster_intersection() (1D GMM density × 1D GMM susceptibility)
```

Dependencies: numpy, scikit-learn, scikit-image

### 8.2 lithoseed (`python/lithoseed/`)

Full pipeline orchestrator:

```python
# _surfaces.py — extract_contact_surfaces(), extract_group_volumes(), cleanup_label_fragments()
#                decimate_mesh(), enforce_stratigraphic_ordering(), separate_connected_components()
# _export.py — write_contact_ts(), write_contact_obj(), write_volume_ts(), write_volume_obj(),
#              write_cluster_csv(), etc.
# _ini.py — write_ini_config() → resolved_config.ini for the C++ inversion executable
# _pipeline.py — run_extract(), run_extract_from_labels(), run_extract_group_volumes()
```

Pipeline flow:
```
SimPEG .msh + .mod files
  → cluster_intersection() [1D GMM density × 1D GMM susceptibility]
  → extract_contact_surfaces() [marching cubes + neighbor filter, for INI contacts]
  → run_extract_group_volumes() [voxel-face boundary extraction, non-overlapping closed volumes]
  → cleanup_label_fragments() [reassign small disconnected label fragments]
  → write_contact_ts/obj() + write_volume_ts/obj() [GOCAD + OBJ export]
  → write_ini_config() [C++ INI config]
  → forrestania_invert.exe resolved_config.ini
```

Volume extraction uses voxel-face boundary extraction: every face in the 3D label grid
is iterated once; shared quads are assigned to both adjacent groups with opposite winding.
This guarantees zero-gap, shared-vertex, perfectly squared-off closed surfaces.
No decimation is applied — shared boundary vertices remain identical between neighbouring volumes.

### 8.3 Standalone Python Tools

**`python/gate_clusters.py`** — Manual property-range gating (no GMM clustering):
```bash
python gate_clusters.py \
    --mesh simpeg_mesh.msh \
    --density simpeg_density.mod \
    --susc simpeg_susceptibility.mod \
    --clusters my_gates.csv \
    --background-density 2.67 \
    --output-dir gated_output/
```
Cluster CSV format:
```csv
cluster_name,density_min,density_max,susc_min,susc_max
Mafic_Upper,2.70,2.85,0.010,0.040
Granite_Gneiss,2.60,2.72,0.000,0.003
```
First-match wins when ranges overlap. Unmatched cells get label -1 (background).
Uses `run_extract_group_volumes()` from lithoseed for volume extraction.

**`python/split_mesh_volumes.py`** — Detect and split disconnected components in .ts meshes:
```bash
python split_mesh_volumes.py volume_group_3.ts --csv cluster_properties.csv
python split_mesh_volumes.py meshes/ --all --csv cluster_properties.csv --formats ts,obj
```
Uses `separate_connected_components()` from lithoseed. Updates cluster_properties.csv
with copied rows for split components (new sequential cluster_ids).

---

## 9. Key INI Config Sections

```ini
[inversion]    max_iterations, tolerance, lambda, omega, control_point_stride
[gravity]      density_min, density_max
[magnetic]     magnetic_weight, inclination, declination, field_nT
[topography]   mode (none|raw|terrain_corrected), dem_file, datum_elevation
[padding]      enable_padding_group, padding_depth, padding_density_*
[regularization] enable_reference_model, lambda_ref
[bounds]       enable_depth_bounds, depth_bound_margin
[data]         cluster_csv, gravity_data, magnetic_data, constraints
[output]       output_dir
```

---

## 10. Build

```powershell
# Full build (all modules + executable)
.\build\build.bat

# Single module
cd modules\litho-core
qmake litho-core.pro && nmake release
```

Requirements: Visual Studio 2022, C++17, Eigen 3.4.0 (vendor/eigen/), Qt 6.2+.

---

## 11. Code Conventions

- **Namespace**: Everything in `litho_invert`
- **Shared pointers**: `LithologyModel` and `SurfaceMesh` via `shared_ptr`
- **Forward model pattern**: compute(params) → model->applyParameterVector(params) → physics
- **parameterCount()**: MUST return `m_model->totalDofCount()` (NOT vertexCount())
- **Debug prints**: Guard with `#ifdef LITHO_INVERT_DEBUG`
- **Index types**: `Eigen::Index` for Eigen, `size_t` for STL, cast explicitly
- **Density safety**: Forward models warn when `LithoGroup::density == 0.0` (uninitialized)

---

## 12. Key Bug Fixes

### From original codebase

1. **MT misfit weights** — fixed: weight now applied to residual in both evaluateComponents() and mtMisfit()
2. **Double gravity evaluation** — fixed: evaluateComponents() reuses single compute() call
3. **Constraint sign bug** — fixed: generateClusterIdConstraints used negative model coords instead of positive-down depths
4. **downsampleVertexGradient stride=0** — fixed: stride=0 now maps vertices directly to DOFs
5. **EM parameterCount()** — fixed: returns totalDofCount() instead of vertexCount() (was SIGSEGV)
6. **L-BFGS-B convergence** — fixed: uses ||x|| per Byrd et al. 1995
7. **Side-wall boundary detection** — extracted into shared boundary_loop.h (was 6 duplicate copies)

### 2026-05-27 to 2026-06-02

8. **`.mod` writer (2026-05-27)** — `run_forrestania_e2e.py` writes one value per line instead of all on one line. Geoscience Analyst compatibility.
9. **`.mod` reader (2026-05-27)** — `cluster_api/_io.py` `_read_property()` detects and skips the 3-integer header (`nx ny nz`) that UBC-format `.mod` files prepend.
10. **`cleanup_label_fragments()` sentinel leak (2026-05-29)** — Sentinel label (-1) in inactive cells was winning neighbor majority votes at the bottom boundary of the active region, causing scattered holes in volume meshes. Fixed by checking `nb_lab >= 0` before considering a neighbor label. Critical for zero-gap volume meshes.

---

## 13. When Working on This Codebase

1. Read this file first to get complete context
2. Read the specific module's CLAUDE.md for its API and build commands
3. Check `DUPLICATE_TRACKING.txt` before modifying any code that might have copies
4. The Forrestania `final/` output files are the gold standard for exporter correctness
5. Each module can be worked on independently — don't load all CLAUDE.md files at once

