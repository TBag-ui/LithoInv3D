# Forrestania — Regional Gravity-Magnetic Joint Inversion

## Purpose

Regional-scale litho-constrained joint inversion of the Forrestania greenstone belt,
Western Australia. Gravity + magnetic data, 4 lithology groups (Granite Gneiss,
Mafic Complex 1+2, Massive Sulfide). The reference implementation whose `final/`
output files define the gold standard for exporter correctness.

## Pipeline

```
Raw gravity CSV + magnetic CSV
  → [1] SimPEG joint inversion (run_joint_inversion.py)
  → [2] GMM clustering → lithology groups
  → [3] Volume extraction + fragment cleanup (via lithoseed)
  → [4] Mesh component splitting (split_mesh_volumes.py)
  → [5] INI config generation
  → [6] forrestania_invert.exe resolved_config.ini
  → [7] Export: TS, OBJ, UBC-GIF, closed volumes
```

## Quick Start

```bash
# Step 1: Run SimPEG joint inversion
python run_joint_inversion.py

# Step 2: Extract starting model with manual gating
python ../python/gate_clusters.py \
    --mesh inversion_output/forrestania_joint.msh \
    --density inversion_output/forrestania_density.mod \
    --susc inversion_output/forrestania_susceptibility.mod \
    --clusters my_gates.csv \
    --output-dir gated_output/

# Step 3: Split multi-component meshes
python ../python/split_mesh_volumes.py gated_output/meshes/ \
    --all --csv gated_output/cluster_properties.csv --formats ts,obj

# Step 4: Run C++ litho-constrained inversion
cd ../../build/release
./forrestania_invert.exe ../../apps/forrestania/resolved_config.ini
```

## Key Files

| File | Purpose |
|------|---------|
| `run_joint_inversion.py` | SimPEG gravity-magnetic cross-gradient inversion |
| `simpg_to_lithoinv.py` | Bridge: SimPEG output → cluster CSV + surface point clouds |
| `groom_data.py` | Data preparation: coordinate transforms, IGRF, QA |
| `build_tensor_mesh.py` | Tensor mesh construction from survey extents |

## Lithology Model

| Group | Density (g/cm³) | Susc (SI) | Notes |
|-------|-----------------|-----------|-------|
| Mafic Complex 1 | ~2.80 | ~0.02 | Upper mafic |
| Granite Gneiss | ~2.67 | ~0.001 | Felsic host |
| Mafic Complex 2 | ~3.00 | ~0.05 | Lower mafic |
| Massive Sulfide | ~3.50 | ~0.10 | Target |

## Inversion Output Directories

| Directory | Clusters | Status |
|-----------|----------|--------|
| `inversion_output/volumes_voxel/` | original | Raw volume extraction |
| `inversion_output/volumes_voxel_cleaned/` | 11 | After fragment cleanup, has bottom-hole bug |
| `inversion_output/volumes_voxel_cleaned_v2/` | 26 | **Current** — sentinel fix applied + mesh splitting |

## IGRF (epoch 2026.0)

- F = 58874 nT
- I = -66.2° (upward, southern hemisphere)
- D = -0.1°

## Gold Standard

The output files in the original `project_Forrestania/Inversion/final/` directory
are the reference for C++ exporter correctness. Any changes to `InversionExporter`
should be verified against these outputs.

## Dependencies

- Python: SimPEG, discretize, numpy, scipy, pandas, matplotlib
- Python: cluster_api, lithoseed (local packages)
- C++: All litho-* modules

## Data Sources

- Gravity: `Forrestania_Gravity_Station_trim_.csv`
- Magnetics: `Forrestania_mag_at_gravity_stations.csv`
- DEM: `Forrestania_SRTM1 Australia_MGA50` (ERS format)
- IGRF: `geomagnetic_field.txt`

