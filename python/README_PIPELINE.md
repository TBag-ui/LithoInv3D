# LithoInv — Python Pipeline Documentation

Complete reference for all Python scripts in the LithoInv project.
These scripts automate the full workflow from raw survey data to
presentation-ready inversion results.

---

## Pipeline Overview

```
Raw Survey Data (CSV)
    │
    ▼
┌─────────────────────────────────────────────┐
│  forrestania_master_pipeline.py  (PRE)       │
│  ═══════════════════════════════════════     │
│  Steps 1-17:                                 │
│    1.  Data preparation                      │
│    2.  Build tensor mesh                     │
│    3.  SimPEG gravity inversion              │
│    4.  SimPEG magnetic inversion             │
│    5.  Absolute density conversion           │
│    6.  GA-compatible .mod export             │
│    7.  UTM coordinate copies                 │
│    8.  MeshData conversion                   │
│    9.  GMM intersection clustering           │
│    10. Volume extraction (voxel-face)        │
│    11. QC plots (convergence, scatter, maps) │
│    12. Split disconnected volumes            │
│    13. Single-property clustering            │
│    13b. Perturbed/user clustering            │
│    13c. Split all volumes                    │
│    14. INI config generation                 │
│    15. Presentation assembly                 │
│    16. Data lineage                           │
│    17. Convenience copy                       │
│                                               │
│  Output: presentation/                        │
└──────────────┬──────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────┐
│  forrestania_post_inversion.py  (POST)       │
│  ═══════════════════════════════════════     │
│    1. Generate INI config                    │
│    2. Launch forrestania_invert.exe          │
│    3. Parse convergence → CSV                │
│    4. Translate final surfaces → UTM         │
│    5. Collect outputs                        │
│                                               │
│  Output: inversion_output/<cluster>/          │
└──────────────────────────────────────────────┘
```

---

## Quick Start

```powershell
# 1. Run the full pre-processing pipeline
cd python
python forrestania_master_pipeline.py

# 2. Run the C++ litho inversion on combination clusters
python forrestania_post_inversion.py combination

# 3. Run on user-defined clusters
python forrestania_post_inversion.py user

# 4. Run both in sequence
python forrestania_post_inversion.py both
```

### Output Structure

```
apps/forrestania/presentation/
├── simpeg/
│   ├── simpeg_mesh.msh                    # UBC-format tensor mesh (34x31x8)
│   ├── simpeg_density.mod                  # Density contrast (GA format, no header)
│   ├── simpeg_density_absolute.mod         # Absolute density = contrast + 2.67
│   ├── simpeg_susceptibility.mod           # Susceptibility (GA format)
│   └── data_coords/                        # UTM-coordinate copies
├── volumes/
│   ├── combination_clusters/               # GMM intersection (density x susc)
│   ├── density_clusters/                   # Density-only GMM
│   ├── susc_clusters/                      # Susceptibility-only GMM
│   └── user_clusters/                      # Manual lithology gates
├── volumes_split/                          # All meshes split by component
├── configs/
│   ├── baseline/resolved_config.ini
│   └── perturbed/resolved_config.ini
├── figures/                                # 8 QC plots
├── lineage/                                # 10-step data provenance
├── convenience/                            # Flat copy of all key outputs
└── inversion_output/
    ├── combination/
    │   ├── convergence.csv                 # 753 iteration lines
    │   ├── iterations/                     # 84 geometry iteration dirs
    │   ├── final_model_space/              # Final surfaces (local coords)
    │   └── final_utm/                      # Final surfaces (UTM)
    └── user/
        └── ...
```

---

## Script Reference

### `forrestania_master_pipeline.py`

**Purpose**: Complete pre-processing pipeline from raw data to C++-ready INI configs.

**Usage**: `python forrestania_master_pipeline.py`

**Key constants** (editable at top of file):
| Constant | Default | Description |
|----------|---------|-------------|
| `CELL_SIZE` | 200.0 m | Tensor mesh cell size |
| `DEPTH_CORE` | 1000.0 m | Depth below min topography |
| `MAX_ITER` | 100 | SimPEG inversion max iterations |
| `GRAV_UNCERTAINTY` | 0.005 mGal | Data uncertainty (tight = more iterations) |
| `MAG_UNCERTAINTY` | 5.0 nT | Data uncertainty |
| `N_DENSITY_CLUSTERS` | 4 | Density GMM bins |
| `N_SUSC_CLUSTERS` | 3 | Susceptibility GMM bins |
| `BACKGROUND_DENSITY` | 2.67 g/cm³ | Added to contrast for absolute density |
| `CLEAN_MIN_CELLS` | 25 | Min cells per fragment |

**Clustering modes** (controlled by `FORRESTANIA_GATES`):
- `None` — GMM density-only with 8 clusters (current default for user_clusters)
- List of tuples — manual lithology gates (name, d_min, d_max, s_min, s_max)

**Z-axis fix**: The `.mod` writer now correctly writes z top-to-bottom (UBC-GIF
spec) by flipping the SimPEG bottom-to-top ordering via `data_3d[:,:,::-1]`.

**Dependencies**: SimPEG, discretize, numpy, scipy, scikit-learn, scikit-image,
pandas, matplotlib, cluster_api, lithoseed.

---

### `forrestania_post_inversion.py`

**Purpose**: Launch the C++ litho-constrained inversion and collect outputs.

**Usage**: `python forrestania_post_inversion.py [combination|user|both]`

**Per-cluster settings** (in `CLUSTER_CONFIGS` dict):
| Key | Property Inversion |
|-----|--------------------|
| `combination` | ON (`enable_property_inversion = true`) |
| `user` | OFF (`enable_property_inversion = false`) |

**Automatic outputs**:
- `convergence.csv` — parsed from C++ stdout (iter, rms, objective, DW stats)
- `iterations/` — all per-iteration surface exports
- `final_utm/` — final surfaces translated to UTM MGA Zone 50
- `final_model_space/` — final surfaces in local coordinates

**Coordinate translation**: Adds `LOCAL_ORIGIN_X/Y` and `Z_DATUM` to model-space
vertices. Handles GOCAD TSurf VRTX format: `VRTX <id> <x> <y> <z>`.

---

### `validate_simpeg_mesh.py`

**Purpose**: Load the SimPEG mesh, place synthetic density/susceptibility bodies,
compute forward gravity + magnetic responses, and verify the mesh produces
physically plausible anomalies.

**Usage**: `python validate_simpeg_mesh.py`

**Output**: `presentation/mesh_validation/` — synthetic data CSVs + QC plots.

The synthetic model places four bodies:
1. Shallow felsic block (low density, -0.15 g/cm³)
2. Intermediate mafic lens (high density, +0.30 g/cm³)
3. Deep ultramafic body (very high density, +0.50 g/cm³)
4. Magnetic pipe (strong susceptibility, 0.15 SI)

---

### `gate_clusters.py`

**Purpose**: Manual property-range gating — define per-cluster density/susceptibility
windows and extract closed volumes. Replaces GMM clustering with user-specified ranges.

**Usage**:
```bash
python gate_clusters.py \
    --mesh simpeg_mesh.msh \
    --density simpeg_density.mod \
    --susc simpeg_susceptibility.mod \
    --clusters my_gates.csv \
    --background-density 2.67 \
    --output-dir gated_output/
```

**Cluster CSV format**:
```csv
cluster_name,density_min,density_max,susc_min,susc_max
Mafic_Upper,2.70,2.85,0.010,0.040
Granite_Gneiss,2.60,2.72,0.000,0.003
```

First-match wins when ranges overlap. Unmatched cells get label -1.
Optional: `--auto-cluster-remainder N` to GMM-cluster unassigned cells.

---

### `split_mesh_volumes.py`

**Purpose**: Detect and split disconnected components in .ts mesh files.
Each component becomes an independent volume file.

**Usage**:
```bash
python split_mesh_volumes.py volume_group_3.ts --csv cluster_properties.csv
python split_mesh_volumes.py meshes/ --all --csv cluster_properties.csv --formats ts,obj
```

Uses `separate_connected_components()` from lithoseed. Updates cluster_properties.csv
with copied rows for split components (new sequential cluster_ids).

---

### `run_forrestania_e2e.py`

**Purpose**: Original end-to-end pipeline (SimPEG → clustering → volumes → INI).
Largely superseded by `forrestania_master_pipeline.py` but kept for reference.

**Usage**: `python run_forrestania_e2e.py`

Uses the OLD mesh-building logic (with `+2*cell_size` overhead padding on z_top).
The new pipeline fixes this.

---

### `generate_single_property_volumes.py`

**Purpose**: Generate single-property (density-only or susceptibility-only)
cluster volumes using GMM, with scatter plots.

**Usage**: Import from this script or run standalone.

---

### `recluster_forrestania.py`

**Purpose**: Re-run clustering on existing SimPEG models with different parameters.

---

### `regenerate_volumetric_ini.py`

**Purpose**: Regenerate INI configs from existing cluster volumes.

---

### `run_diagnostic_cluster_volumes.py`

**Purpose**: Diagnostic script for inspecting cluster volume properties.

---

## Module Dependencies

```
forrestania_master_pipeline.py
├── cluster_api/
│   ├── _cluster.py    — cluster_intersection(), cluster_lithology()
│   ├── _io.py         — MeshData, load_inversion()
│   └── _surface.py    — ContactSurface
├── lithoseed/
│   ├── _pipeline.py   — run_extract_from_labels(), run_extract_group_volumes()
│   ├── _surfaces.py    — extract_contact_surfaces(), cleanup_label_fragments()
│   ├── _export.py     — write_contact_ts(), write_volume_ts(), write_cluster_csv()
│   └── _ini.py        — write_ini_config()
└── gate_clusters.py   — gate_cells(), build_summary()
```

---

## Environment

- **Python**: 3.9+
- **Key packages**: SimPEG 0.22+, discretize 0.10+, numpy, scipy, scikit-learn,
  scikit-image, pandas, matplotlib

---

## See Also

- `../README.md` — project overview, build instructions, citation
- `../docs/developer.md` — C++ module architecture
- `../docs/scientific_reference.md` — mathematical formulation
- `../docs/key_parameters.md` — inversion parameter guide
- `../CLAUDE.md` — AI coding assistant context
