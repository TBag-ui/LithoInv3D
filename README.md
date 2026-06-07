# LithoInv — Lithology-constrained 3D joint inversion of gravity and magnetics using closed surface volumes

From voxel models to geological contact surfaces.

**LithoInv** recovers 3D subsurface lithological boundary surfaces from gravity and
magnetic data. Instead of inverting for voxel property values, it inverts for the
positions of explicit triangulated surfaces that separate distinct lithological
units. The output is a geological contact model — surfaces that can be imported
directly into 3D geological modelling packages.

---

## Method Overview

The inversion pipeline has three stages:

1. **SimPEG voxel inversion** — Separate gravity and magnetic inversions on a
   shared tensor mesh produce 3D density-contrast and susceptibility models.

2. **GMM clustering + surface extraction** — Gaussian Mixture Model intersection
   clustering (density × susceptibility) partitions the voxel model into lithology
   groups.  Closed volumes are extracted via voxel-face boundary enumeration:
   every face in the 3D label grid is visited once; shared quads between adjacent
   groups are assigned to both with opposite winding, guaranteeing zero-gap,
   shared-vertex, perfectly edge-matched surfaces with no decimation.

3. **L-BFGS-B geometry inversion** — The contact surfaces are loaded into a C++
   engine that minimises a joint gravity + magnetic objective function using
   bound-constrained L-BFGS-B optimisation.  Control-point degree-of-freedom
   downsampling enables practical inversion of large surface meshes.  The forward
   model uses analytical polyhedron formulae (Nagy, 1966; Okabe/Plouff, 1976–1979)
   with full analytical Z-derivatives and finite-difference XY-derivatives.

Supporting features include borehole depth constraints, reference-model
regularisation, topography, lateral padding, property inversion (density and
susceptibility), and XYZ_FREE vertex freedom for full 3-DOF surface inversion.

---

## Quick Start

### Build

Requirements: Visual Studio 2022 or MinGW-w64 11.2, C++17, Eigen 3.4.0,
Qt 6.2+ (qmake).

```powershell
.\build\build.bat
```

### Run

```powershell
# Pre-processing (SimPEG + clustering + volumes)
cd python
python forrestania_master_pipeline.py

# C++ litho-constrained inversion
cd build\release
.\forrestania_invert.exe config.ini
```

---

## Repository

https://github.com/TBag-ui/LithoInv3D

---

## How to Cite

> Thomas Bagley. LithoInv: Lithology-constrained 3D inversion of gravity and
> magnetics, preprint (2026). https://eartharxiv.org/repository/view/XXXXX/

*(Update the EarthArXiv URL after the preprint is accepted.)*

---

## License

MIT License — see [LICENSE](LICENSE) for full terms.

Copyright (c) 2026 Thomas Bagley

---

## Trademark

**LithoInv**™ is a trademark of Thomas Bagley.

---

## Consulting

Consulting and custom development are available.
Contact: [bagleyt0@gmail.com](mailto:bagleyt0@gmail.com)
