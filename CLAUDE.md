# CLAUDE.md — LithoInvert3D Master Project

## Quick Start

To start working on this codebase in a fresh Claude Code session:

> Read `docs/llm_context.md` for complete project context. Then read the CLAUDE.md of the specific module you'll be working on.

## Build

```powershell
.\build\build.bat
```
Produces `build/release/forrestania_invert.exe`. Requires: Visual Studio 2022, C++17, Eigen 3.4.0 (vendored), Qt 6.2+ with qmake.

## Run an Inversion

```powershell
cd build\release
.\forrestania_invert.exe <path-to-config.ini>
```

Key test configs:
- `datasets/Voiseys_Bay/SyntheticDomeTest/dome_test.ini` — gravity-only Z_ONLY baseline
- `datasets/Voiseys_Bay/SyntheticDomeTest/dome_test_combined.ini` — all features

## Module Map

```
modules/
├── litho-core/          — common types, geometry, stats (NO deps)
├── litho-surface/       — SurfaceMesh, vertices, DOFs (depends: core)
├── litho-model/         — LithoGroup, LithologyModel (depends: core, surface)
├── litho-forward/       — gravity + magnetic forward (depends: core, surface, model)
├── litho-em/            — EM forward stubs (depends: core, surface, model)
├── litho-inversion/     — optimizer, objective, runner (depends: all above)
├── litho-io/            — loaders, exporters, INI parser (depends: core, surface, model)
└── litho-regularization/ — smoothness, reference model (depends: core, surface, model)
```

Dependency order for building: core → surface → model → forward/em/io/regularization → inversion

## Python

```bash
# Run all tests
python -m pytest python/lithoseed/tests/ python/lithoseed/tests/integration/ -v

# Extract contacts from SimPEG results
python -m lithoseed extract --mesh model.msh --density density.mod --n-clusters 4
```

### Standalone Tools

| Script | Purpose |
|--------|---------|
| `python/gate_clusters.py` | Extract volumes by manual property-range gating (no GMM clustering) |
| `python/split_mesh_volumes.py` | Split multi-component .ts meshes into separate files + update CSV |
| `python/run_forrestania_e2e.py` | Full Forrestania pipeline: SimPEG → clustering → volumes → INI |

See `python/lithoseed/CLAUDE.md` and `python/cluster_api/CLAUDE.md` for details.

## Working on a Specific Module

Each module has its own CLAUDE.md. When working on `litho-forward`:
> Read `modules/litho-forward/CLAUDE.md`, then make changes.

The master CLAUDE.md (this file) only needs to route to submodules. Do not load all module CLAUDE.md files at once — pick the one you need.

## Project Conventions

- **Namespace**: `litho_invert` for all C++ code
- **Shared pointers**: `LithologyModel` and `SurfaceMesh` via `std::shared_ptr`
- **Error handling**: `std::runtime_error` for fatal errors, `std::cerr` for diagnostics
- **Debug prints**: Guard with `#ifdef LITHO_INVERT_DEBUG`
- **Index types**: `Eigen::Index` for Eigen, `size_t` for STL, cast explicitly
- **Coordinate convention**: Z positive UP (surface z=0, deep z=-5000)
- **Constraint depths**: POSITIVE DOWN (z_top=100 means 100m below surface)

