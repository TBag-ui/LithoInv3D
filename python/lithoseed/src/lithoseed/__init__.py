"""
LithoSeed — generate 3D starting models for the C++ litho-constrained inversion
from SimPEG voxel inversions.

Pipeline:
  1. Run a configurable SimPEG gravity/magnetic inversion (tensor mesh)
  2. Cluster recovered properties into lithology groups (GMM or intersection clustering)
  3. Extract 3D contact surfaces as triangulated meshes (marching cubes + neighbor filter)
  4. Extract non-overlapping closed volumes via voxel-face boundary extraction
     (shared vertices, zero-gap, perfectly squared-off edges)
  5. Export GOCAD TSurf (.ts) and OBJ meshes, plus LithoInv CSV + INI

Usage:
  python run_forrestania_e2e.py  # Full pipeline
  from lithoseed._pipeline import run_extract, run_extract_group_volumes
"""

from ._version import __version__
from ._config import SimPEGConfig, build_parser
from ._surfaces import ContactSurface, extract_contact_surfaces, decimate_mesh
from ._surfaces import extract_group_volumes, cleanup_label_fragments
from ._export import (
    write_contact_ts,
    write_contact_obj,
    write_cluster_csv,
    write_observed_gravity_csv,
    write_observed_magnetic_csv,
    write_borehole_constraints_csv,
)
from ._ini import write_ini_config
from ._pipeline import run_extract, run_extract_from_labels, run_extract_group_volumes, PipelineResult

__all__ = [
    "__version__",
    "SimPEGConfig",
    "build_parser",
    "ContactSurface",
    "extract_contact_surfaces",
    "extract_group_volumes",
    "cleanup_label_fragments",
    "decimate_mesh",
    "write_contact_ts",
    "write_contact_obj",
    "write_cluster_csv",
    "write_observed_gravity_csv",
    "write_observed_magnetic_csv",
    "write_borehole_constraints_csv",
    "write_ini_config",
    "run_extract",
    "run_extract_from_labels",
    "run_extract_group_volumes",
    "PipelineResult",
]
