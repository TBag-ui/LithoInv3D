# litho-core — Common Types, Geometry, and Statistics

## Purpose

Foundation module for all other litho-* modules. Provides Eigen typedefs, data
structs, analytical geometry primitives, boundary-loop extraction, and statistical
diagnostics. Zero dependencies beyond Eigen.

## API

### common.h — Data Types
```cpp
#include <litho_invert/core/common.h>

// Eigen wrappers
using litho_invert::Vector3d;   // Eigen::Vector3d
using litho_invert::VectorXd;   // dynamic double vector
using litho_invert::MatrixXd;   // dynamic double matrix
using litho_invert::Index;      // Eigen::Index (ptrdiff_t)

// Observation data
struct GravityPoint { Vector3d position; double g_obs, g_std; };
struct MagneticPoint { Vector3d position; double t_obs, t_std; };

// Constraints (borehole intervals)
struct Constraint { Vector3d position; double z_top, z_bottom; int litho_group_id; };

// Topography
enum class TopographyMode { None, Raw, TerrainCorrected };
struct TopographyConfig { TopographyMode mode; std::string demFile; /*...*/ };
```

**Coordinate conventions:**
- Vertices: Z positive UP (surface z=0, deep z=-5000)
- Constraint depths: POSITIVE DOWN (z_top=100 means 100m below surface)

### geometry.h — Nagy/Okabe Primitives
```cpp
#include <litho_invert/core/geometry.h>

double solidAngle(a, b, c, obs);             // Ω subtended by triangle at obs
double lineIntegralTerm(a, b, obs);          // edge line integral L
double surfaceIntegralTerm(a, b, c, obs);    // face surface integral
double tripleProduct(a, b, c);              // a·(b×c)
double tetraVolume(a, b, c, d);             // signed tetrahedron volume
```

### stats.h — Durbin-Watson
```cpp
#include <litho_invert/core/stats.h>

DurbinWatsonResult { double dw_x, dw_y; };
DurbinWatsonResult computeDurbinWatson(positions, residuals);
// DW ≈ 2.0 → random, DW ≈ 0.0 → underfitting, DW ≈ 4.0 → oscillating
```

## Build

```powershell
cd modules/litho-core
qmake litho-core.pro
nmake release
```

## Dependencies

- Eigen 3.4.0 (vendor/eigen/) — header only
- C++17 standard library

## Tests

```powershell
cd modules/litho-core/tests
qmake tests.pro && nmake release && release\tests.exe
```

Tests: solidAngle, lineIntegralTerm, tripleProduct, tetraVolume, boundary edge detection,
Durbin-Watson on known patterns.
