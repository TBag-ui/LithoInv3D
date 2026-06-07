# cluster_api — GMM Clustering and Mesh I/O for SimPEG Outputs

## Purpose

Reads SimPEG inversion output (tensor mesh + property models), clusters recovered
physical properties into lithology groups via Gaussian Mixture Models, and provides
utilities for working with the resulting clustered volumes.

This is the **authoritative copy** of the clustering code. The old
`project_Forrestania/api/` was a byte-identical duplicate — deleted during
restructuring. All regional applications should import from this package.

## API

```python
from cluster_api import load_inversion, MeshData
from cluster_api import cluster_lithology, LithologySummary

# === I/O ===
mesh, density, susceptibility = load_inversion(
    "model.msh",           # SimPEG tensor mesh
    "density.mod",         # density per active cell (.mod or .xyz)
    "susc.mod",            # susceptibility (optional — returns NaN array if None)
)

# MeshData attributes
mesh.nx, mesh.ny, mesh.nz         # grid dimensions
mesh.x0, mesh.y0, mesh.z0         # origin (z0 = top surface)
mesh.dx, mesh.dy, mesh.dz         # cell sizes
mesh.n_active                     # number of active cells
mesh.ix, mesh.iy, mesh.iz         # 1-based indices of active cells [n_active]
mesh.x_center, y_center, z_center # cell-centre coordinates [n_active]

# Map per-active-cell values to dense 3D grid
grid = mesh.to_3d_grid(values, fill=-1)  # → (nx, ny, nz) numpy array

# Node coordinate arrays for marching cubes
xn, yn, zn = mesh.node_arrays()

# === Clustering ===
labels, summary = cluster_lithology(
    density,              # (N,) array
    susceptibility,       # (N,) array
    n_clusters=4,         # number of lithology groups
    sentinel=-9999.0,     # values below sentinel+1 are excluded (label=-1)
    random_state=42,
)

# labels: (N,) int array — -1=sentinel, 1..n_clusters (ordered by density)
# summary: LithologySummary dataclass
#   summary.labels       — [1, 2, ..., n_clusters]
#   summary.counts       — cells per group
#   summary.density_mean — mean density per group [g/cc]
#   summary.density_std  — std density per group
#   summary.susc_mean    — mean susceptibility per group [SI]
#   summary.susc_std
#   summary.n_groups
#   summary.to_csv(path) — export cluster properties CSV
```

## Algorithm

1. Filter: cells with density < sentinel+1 are excluded (outside ROI)
2. Standardize: (density, susceptibility) → StandardScaler
3. Cluster: GaussianMixture (full covariance, n_init=10, max_iter=500)
4. Order: clusters sorted by mean density (lowest → group 1, highest → group N)
5. Remap: labels reassigned to reflect ordered groups

## Dependencies

- numpy
- scikit-learn (GaussianMixture, StandardScaler)

## Tests

```bash
# Tests not yet implemented — tests/ directory exists but is empty
python -m pytest python/cluster_api/tests/ -v
```

## Build (editable install)

```bash
pip install -e python/cluster_api
```

## Recent Fixes

- **`.mod` reader (2026-05-27)**: `_read_property()` now detects the 3-integer header (`nx ny nz`) that UBC-format `.mod` files prepend, and skips it. Previously assumed all `.mod` files had no header.

## C++ Header

`cluster_api/surface_export.h` provides standalone C++ TS/OBJ writers for use
in scripts that can't link the full litho-io module. Uses raw Eigen vectors
(Vector3d + Vector3i) rather than SurfaceMesh.

**Note**: This is a convenience header for standalone use. The authoritative
TS/OBJ writer is `InversionExporter` in `modules/litho-io`.

