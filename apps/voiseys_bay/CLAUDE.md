# Voisey's Bay — Synthetic Ni-Cu-Co Test Case

## Purpose

Synthetic inversion test case based on the Voisey's Bay Ni-Cu-Co deposit, Labrador.
Exercises all features of the inversion framework: gravity + magnetics + EM joint
inversion, property inversion, topography, padding, constraints, XYZ_FREE vertices,
reference model regularization, and depth bounds.

## True Model

| Group | Density (g/cm³) | Susc (SI) | Conductivity (S/m) | Geometry |
|-------|-----------------|-----------|---------------------|----------|
| Anorthosite (host) | 2.70 | 0.001 | 1e-4 | Regional background |
| Troctolite (envelope) | 2.95 | 0.020 | 1e-2 | Ellipsoidal lens, ~250m thick |
| Massive Sulfide (ore) | 4.00 | 0.300 | 1.0 | Nested pod, ~150m thick |

## Survey Geometry

- 500m × 500m area, 11×11=121 gravity + magnetic points on 50m grid
- 1 airborne TEM source at (0, 0, 30), 8 time gates (50µs–6.4ms)
- 4 MT stations at corners, 6 frequencies (0.001–100 Hz)
- Initial model: 2 flat surfaces at -300m and -400m

## Quick Start

```powershell
# Build
cd ..\..\build
.\build.bat

# Run baseline test (Z_ONLY, gravity only)
cd build\release
.\csv_invert_fixed.exe ..\..\datasets\Voiseys_Bay\SyntheticDomeTest\dome_test.ini

# Run full-feature test (XYZ_FREE + reference model + depth bounds)
.\csv_invert_fixed.exe ..\..\datasets\Voiseys_Bay\SyntheticDomeTest\dome_test_combined.ini

# Override max iterations for quick test
.\csv_invert_fixed.exe dome_test.ini --max-iter=5
```

## Test Configs

| Config | Features |
|--------|----------|
| `dome_test.ini` | Z_ONLY baseline |
| `dome_test_xyzfree.ini` | XYZ_FREE vertex freedom |
| `dome_test_refmodel.ini` | Reference model (λ=0.1) |
| `dome_test_strong_ref.ini` | Strong reference (λ=100) |
| `dome_test_tight_bounds.ini` | Hard depth bounds (±10m) |
| `dome_test_combined.ini` | All features combined |

## Key Source Files

| File | Purpose |
|------|---------|
| `main.cpp` | Entry point — INI or legacy config, run, export |
| `generate_synthetic.h/cpp` | True model + synthetic data generation |
| `cluster_loader.h/cpp` | Lithosplitter CSV cluster parser |

## Dependencies

- All litho-* C++ modules (linked via top-level build)
- Eigen 3.4.0

## Expected Results

- Gravity RMS should converge toward <1 mGal on the synthetic dome
- DW statistics should approach 2.0 (random residuals)
- Recovered surfaces should match the ellipsoidal true model geometry
