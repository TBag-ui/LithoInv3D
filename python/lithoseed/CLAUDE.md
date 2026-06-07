# lithoseed — Starting-Model Pipeline for LithoInvert3D

## Purpose

Generates 3D starting models for the C++ litho-constrained inversion from SimPEG
voxel inversions. Two output types:

1. **Contact surfaces** (marching cubes + neighbor filter) — interfaces between
   adjacent lithology groups, used for the C++ inversion INI config.
2. **Closed volumes** (voxel-face boundary extraction) — non-overlapping, zero-gap
   per-group volumes with shared vertices at contacts. No decimation applied.

## Pipeline

```
Raw data (CSV)
    → [1] SimPEG gravity + magnetic inversions (tensor mesh)
    → [2] Intersection clustering (1D GMM density × 1D GMM susceptibility)
    → [3] cleanup_label_fragments() — reassign small disconnected label fragments
    → [4] extract_contact_surfaces() — marching cubes + neighbor filter (for INI)
    → [5] extract_group_volumes() — voxel-face boundary extraction (closed volumes)
    → [6] Export: GOCAD TSurf + OBJ + CSVs + INI config
    → [7] forrestania_invert.exe resolved_config.ini
```

## Quick Start

```bash
# Run full Forrestania e2e pipeline
python run_forrestania_e2e.py
```

## Module Map

| File | Purpose |
|------|---------|
| `_surfaces.py` | `extract_contact_surfaces()`, `extract_group_volumes()`, `cleanup_label_fragments()`, `decimate_mesh()`, `enforce_stratigraphic_ordering()`, `separate_connected_components()` |
| `_export.py` | `write_contact_ts()`, `write_contact_obj()`, `write_volume_ts()`, `write_volume_obj()`, `write_cluster_csv()`, etc. |
| `_ini.py` | `write_ini_config()` → resolved_config.ini |
| `_pipeline.py` | `run_extract()`, `run_extract_from_labels()`, `run_extract_group_volumes()` |

## Standalone Tools (import from lithoseed but don't modify it)

| Script | Purpose |
|--------|---------|
| `../gate_clusters.py` | Manual property-range gating — define per-cluster density/susc windows, extract closed volumes |
| `../split_mesh_volumes.py` | Detect & split disconnected components in .ts meshes, update cluster_properties.csv |
| `../run_forrestania_e2e.py` | Full Forrestania pipeline end-to-end |

## Recent Bug Fixes

- **`cleanup_label_fragments()` (2026-05-29)**: Sentinel label (-1) in inactive cells was winning neighbor majority votes at the bottom boundary, causing holes in volume meshes. Fixed by checking `nb_lab >= 0` before considering a neighbor label.
- **`.mod` reader (2026-05-27)**: `cluster_api/_io.py` now detects and skips the 3-integer header (nx ny nz) that UBC-format `.mod` files prepend.
- **`.mod` writer (2026-05-27)**: `run_forrestania_e2e.py` now writes one value per line instead of all values on one line, for Geoscience Analyst compatibility.

## API

```python
from lithoseed import SimPEGConfig, build_parser
from lithoseed import ContactSurface, extract_contact_surfaces, decimate_mesh
from lithoseed import write_contact_ts, write_contact_obj, write_cluster_csv
from lithoseed import write_ini_config
from lithoseed._pipeline import run_extract, PipelineResult

# === Config ===
config = SimPEGConfig(
    mode="joint",
    mesh_type="octree",
    base_cell_size=30.0,
    depth_core=1000.0,
    gravity_csv="gravity.csv",
    magnetic_csv="magnetic.csv",
    n_clusters=4,
    target_vertices_per_contact=500,
    output_dir="lithoseed_output",
)
config.validate()  # raises ValueError on invalid combinations

# === Contact Extraction ===
contacts = extract_contact_surfaces(
    mesh, labels, group_order=[1, 2, 3, 4],
    local_origin=(x0, y0, z0),
)
# contacts: list of ContactSurface
#   cs.group_above, cs.group_below — int
#   cs.vertices — (N, 3) numpy array, world coordinates
#   cs.faces — (M, 3) numpy array, 0-based indices
#   cs.median_depth — float

# === Decimation ===
dec_v, dec_f = decimate_mesh(vertices, faces, target_vertices=500)
# Vertex-clustering: bins bounding box into cubes, merges to centroids

# === Stratigraphic Ordering ===
contacts = enforce_stratigraphic_ordering(contacts, margin=5.0)
# Ensures surfaces don't cross in 3D

# === Pipeline ===
result = run_extract(
    mesh, density, susceptibility,
    n_clusters=4, target_vertices=500,
    local_origin=(0, 0, 0), z_datum=0,
    output_dir="output",
    export_formats=("ts", "obj"),
    gravity_xyz=..., gravity_obs=..., gravity_std=...,
    magnetic_xyz=..., magnetic_obs=..., magnetic_std=...,
    borehole_constraints=...,
)
# result.contacts, result.group_order, result.ini_path
```

## Output Files

```
lithoseed_output/
  meshes/contact_1_2.ts       (GOCAD TSurf per contact pair)
  meshes/contact_1_2.obj      (OBJ interchange)
  cluster_properties.csv      (group_id, density_mean/std, susc_mean/std)
  observed_gravity.csv        (x,y,z,g_obs)
  observed_magnetic.csv       (x,y,z,t_obs)
  borehole_constraints.csv    (x,y,z_top,z_bottom,litho_group_id)
  resolved_config.ini         (ready for forrestania_invert.exe)
```

## Coordinate Conventions

- **MeshData z0**: top surface, Z decreases downward in index space
- **World coordinates**: Z positive UP (matching C++ LithoInvert3D)
- **Borehole constraints CSV**: Z positive DOWN (matching C++ CSVConstraintLoader)

## Dependencies

- cluster_api (GMM clustering + mesh I/O)
- numpy, scipy, scikit-learn, scikit-image
- Optional: SimPEG + discretize (for `pipeline` and `invert` commands)

## Tests

```bash
# Unit tests (not yet implemented — tests/ directory exists, tests to be written)
python -m pytest python/lithoseed/tests/ -v

# Integration tests (not yet implemented)
python -m pytest python/lithoseed/tests/integration/voiseys_bay/ -v
```

## Build (editable install)

```bash
pip install -e python/lithoseed
```

