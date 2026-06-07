# litho-io — Loaders, Exporters, and Config Parser

## Purpose

All file I/O for the inversion pipeline. Loads observations, surfaces, constraints;
exports results in multiple formats (OBJ, GOCAD TSurf, UBC-GIF, CSV); parses INI
config files.

## API

```cpp
#include <litho_invert/io/io.h>           // umbrella — all loaders + exporters

// === Loaders ===
GravityData CSVGravityLoader::load(path);          // CSV: x,y,z,g_obs,g_std
shared_ptr<SurfaceMesh> OBJSurfaceLoader::load(path); // OBJ: v/f
vector<Constraint> CSVConstraintLoader::load(path);   // CSV: x,y,z_top,z_bottom,litho_id
shared_ptr<LithologyModel> JSONLithoConfigLoader::load(path); // JSON groups
vector<double> loadDEM(path, xPositions, yPositions);  // XYZ grid → bilinear interp

// === Config ===
#include <litho_invert/io/ini_config.h>
class IniConfig {
    bool load(path); void save(path);
    string getString(section, key, default);
    double getDouble(section, key, default);
    int getInt(section, key, default);
    bool getBool(section, key, default);
};

// === Exporters ===
#include <litho_invert/io/exporters.h>
class InversionExporter {
    InversionExporter(outputDir, baseName);
    void setSubfolder(sub);
    void setGroupNaming(groupExportNames);  // e.g. {"Granite", "Basalt", "Sulfide"}

    // Single-surface exports
    void exportTS(mesh, suffix);             // GOCAD TSurf (.ts)
    void exportOBJ(mesh, suffix);            // Wavefront OBJ
    void exportRawVertices(mesh, suffix);    // .txt
    void exportRawTriangles(mesh, suffix);   // .txt
    void exportSurfaceCSV(mesh, suffix);     // x,y,z

    // Interior-only (truncates padding)
    void exportInteriorTS(mesh, suffix);
    void exportInteriorOBJ(mesh, suffix);

    // Closed volumes (top + bottom + side walls)
    void exportClosedVolume(topSurf, botSurf, suffix, bottomDepth);

    // Starting model (all contacts + closed volumes per group)
    void exportStartingModel(model, flatTop, flatBottom, bottomDepth);

    // Bulk
    void exportAll(result, observed, xmin, xmax, ymin, ymax, zmin, zmax, cellSize);

    // Specialized
    void exportUBCGIF(model, ...);           // UBC-GIF mesh (litho codes per cell)
    void exportLithoCSV(model, ...);         // per-cell CSV
    void exportPredictedCSV(observed, predicted);
    void exportLog(result);
};
```

## Export Conventions

### GOCAD TSurf (.ts)
```
GOCAD TSurf 1
HEADER { name: ... }
TFACE
VRTX 1 x y z    (scientific, 8 decimal places)
TRGL v1 v2 v3
END
```

### Closed Volume TS
Top surface (CCW) + bottom surface (reversed winding, vertex offset by nVertsTop)
+ side walls stitched via `extractBoundaryLoops()` from litho-core.

### Dual Export (padding mode)
When `paddingRings > 0`: `exportAll()` produces both survey-area (interior-only)
and full-model (with _full suffix) TS/OBJ sets.

## Gold Standard

The Forrestania `final/` output files in the original codebase were produced by
this exporter and serve as the reference for format correctness.

## Build

```powershell
cd modules/litho-io
qmake litho-io.pro
nmake release
```

## Dependencies

- litho-core (types, boundary_loop for closed-volume side walls)
- litho-surface (SurfaceMesh)
- litho-model (LithologyModel)
