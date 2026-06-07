"""Configuration for LithoSeed: SimPEG inversion + clustering + surface extraction."""

import argparse
from dataclasses import dataclass, field
from typing import Optional, Tuple


@dataclass
class SimPEGConfig:
    """All parameters for running the LithoSeed pipeline.

    Covers SimPEG inversion, GMM clustering, mesh decimation, and export.
    Sensible defaults fix the known issues from run_joint_inversion.py.
    """

    # ---- Mode ----
    mode: str = "joint"  # "gravity" | "magnetic" | "joint"

    # ---- Mesh ----
    mesh_type: str = "octree"  # "octree" | "tensor"
    base_cell_size: float = 30.0  # metres (was 60m in old script)
    depth_core: float = 1000.0  # core mesh depth below surface
    pad_distance_xy: float = 500.0  # lateral padding beyond data extent
    extraction_cell_size: float = 30.0  # tensor grid resolution for surface extraction

    # ---- Data paths ----
    gravity_csv: Optional[str] = None
    magnetic_csv: Optional[str] = None
    dem_csv: Optional[str] = None
    borehole_csv: Optional[str] = None  # litho intervals for constraints

    # ---- Gravity ----
    gravity_uncertainty: float = 0.05  # mGal (was 0.015 — too tight)
    density_lower: float = -0.4  # g/cc
    density_upper: float = 0.6  # g/cc (was missing — critical fix)

    # ---- Magnetics ----
    magnetic_uncertainty: float = 5.0  # nT (was 60 — too loose as floor)
    susc_lower: float = 0.0  # SI
    susc_upper: float = 0.15  # SI

    # ---- IGRF ----
    igrf_f: float = 55000.0  # nT — total field strength
    igrf_i: float = -66.0  # deg — inclination (negative = upward, southern hemisphere)
    igrf_d: float = 0.0  # deg — declination

    # ---- Cross-gradient ----
    enable_cross_gradient: bool = True
    cross_grad_weight: float = -1.0  # auto: cell_size**4 / 16 when < 0

    # ---- Optimization ----
    max_iter: int = 20
    max_irls: int = 10
    beta0_ratio: float = 1.0

    # ---- Downsampling (magnetic data) ----
    mag_downsample: int = 2500  # target number of magnetic points (0 = no downsample)

    # ---- Clustering ----
    n_clusters: int = 4
    cluster_random_state: int = 42

    # ---- Mesh decimation ----
    target_vertices_per_contact: int = 500

    # ---- Depth weighting ----
    depth_weight_beta_grav: float = 2.0  # exponent for gravity depth weighting
    depth_weight_beta_mag: float = 3.0   # exponent for magnetic depth weighting
    depth_weight_z0: float = -1.0        # reference depth (< 0 = auto: base_cell_size)

    # ---- Borehole geophysics ----
    borehole_properties_csv: Optional[str] = None  # downhole density/susc measurements
    borehole_ref_model: bool = True                 # build reference model from property logs
    borehole_local_weight_radius: float = -1.0      # reduce regularization near boreholes (< 0 = auto)

    # ---- Coordinate transform ----
    local_origin_x: float = 0.0
    local_origin_y: float = 0.0
    z_datum: float = 0.0

    # ---- Output ----
    output_dir: str = "lithoseed_output"
    plot: bool = True
    export_formats: Tuple[str, ...] = field(default=("ts", "obj"))

    # ---- Topography ----
    bouguer_density: float = 2.67  # g/cc for terrain correction

    # ---- Borehole constraints (LithoInv side) ----
    # These flow through to the C++ INI config
    borehole_uncertainty: float = 25.0  # metres — tolerance for borehole intersection

    def validate(self) -> None:
        """Raise ValueError on invalid combinations."""
        if self.mode not in ("gravity", "magnetic", "joint"):
            raise ValueError(f"Invalid mode: {self.mode}. Use gravity|magnetic|joint.")
        if self.mesh_type not in ("octree", "tensor"):
            raise ValueError(f"Invalid mesh_type: {self.mesh_type}. Use octree|tensor.")
        if self.base_cell_size <= 0:
            raise ValueError(f"base_cell_size must be positive, got {self.base_cell_size}")
        if self.n_clusters < 2:
            raise ValueError(f"n_clusters must be >= 2, got {self.n_clusters}")
        if self.target_vertices_per_contact < 10:
            raise ValueError(f"target_vertices_per_contact must be >= 10")
        if self.density_lower >= self.density_upper:
            raise ValueError(
                f"density_lower ({self.density_lower}) must be < density_upper ({self.density_upper})"
            )
        if self.susc_lower >= self.susc_upper:
            raise ValueError(
                f"susc_lower ({self.susc_lower}) must be < susc_upper ({self.susc_upper})"
            )
        if self.mode in ("gravity", "joint") and self.gravity_csv is None:
            raise ValueError(f"Mode '{self.mode}' requires gravity_csv.")
        if self.mode in ("magnetic", "joint") and self.magnetic_csv is None:
            raise ValueError(f"Mode '{self.mode}' requires magnetic_csv.")

    @classmethod
    def from_argparse(cls, ns: argparse.Namespace) -> "SimPEGConfig":
        """Build config from parsed CLI arguments."""
        kwargs = {}
        for field_name in [
            "mode", "mesh_type", "base_cell_size", "depth_core",
            "pad_distance_xy", "extraction_cell_size",
            "gravity_csv", "magnetic_csv", "dem_csv", "borehole_csv",
            "gravity_uncertainty", "density_lower", "density_upper",
            "magnetic_uncertainty", "susc_lower", "susc_upper",
            "igrf_f", "igrf_i", "igrf_d",
            "enable_cross_gradient", "cross_grad_weight",
            "max_iter", "max_irls", "beta0_ratio",
            "mag_downsample", "n_clusters", "cluster_random_state",
            "target_vertices_per_contact",
            "depth_weight_beta_grav", "depth_weight_beta_mag", "depth_weight_z0",
            "borehole_properties_csv", "borehole_ref_model",
            "borehole_local_weight_radius",
            "local_origin_x", "local_origin_y", "z_datum",
            "output_dir", "plot", "bouguer_density",
            "borehole_uncertainty",
        ]:
            if hasattr(ns, field_name):
                val = getattr(ns, field_name)
                if val is not None:
                    kwargs[field_name] = val
        return cls(**kwargs)

    def auto_cross_grad_weight(self) -> float:
        """Normalized cross-gradient weight: cell_size^4 / 16.

        The reference SimPEG tutorial uses 30m cells → 30^4 = 810,000.
        The old run_joint_inversion.py used 60m → 60^4 = 12,960,000 (16x too large).
        This normalization makes the weight invariant to cell size choice.
        """
        if self.cross_grad_weight < 0:
            return self.base_cell_size ** 4 / 16.0
        return self.cross_grad_weight


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser for 'lithoseed invert' and 'lithoseed pipeline'."""

    p = argparse.ArgumentParser(
        description="LithoSeed — SimPEG starting-model pipeline for litho inversion",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Mode
    p.add_argument("--mode", default="joint",
                   choices=["gravity", "magnetic", "joint"],
                   help="Physics mode (default: joint)")

    # Mesh
    p.add_argument("--mesh-type", default="octree", choices=["octree", "tensor"])
    p.add_argument("--cell-size", type=float, default=30.0,
                   help="Base cell size in metres (default: 30)")
    p.add_argument("--depth-core", type=float, default=1000.0,
                   help="Core mesh depth in metres (default: 1000)")
    p.add_argument("--pad-distance", type=float, default=500.0,
                   help="XY padding beyond data extent (default: 500)")
    p.add_argument("--extraction-cell-size", type=float, default=30.0,
                   help="Tensor grid cell size for surface extraction (default: 30)")

    # Data
    p.add_argument("--gravity-csv", help="Path to observed gravity CSV")
    p.add_argument("--magnetic-csv", help="Path to observed magnetic CSV")
    p.add_argument("--dem-csv", help="Path to DEM CSV for topography")
    p.add_argument("--borehole-csv", help="Path to borehole lithology intervals CSV")

    # Bounds
    p.add_argument("--density-lower", type=float, default=-0.4)
    p.add_argument("--density-upper", type=float, default=0.6)
    p.add_argument("--susc-lower", type=float, default=0.0)
    p.add_argument("--susc-upper", type=float, default=0.15)

    # Uncertainties
    p.add_argument("--gravity-uncert", type=float, default=0.05,
                   help="Gravity floor uncertainty in mGal (default: 0.05)")
    p.add_argument("--magnetic-uncert", type=float, default=5.0,
                   help="Magnetic floor uncertainty in nT (default: 5.0)")

    # IGRF
    p.add_argument("--igrf-f", type=float, default=55000.0)
    p.add_argument("--igrf-i", type=float, default=-66.0)
    p.add_argument("--igrf-d", type=float, default=0.0)

    # Cross-gradient
    p.add_argument("--no-cross-grad", action="store_true",
                   help="Disable cross-gradient coupling")
    p.add_argument("--cross-grad-weight", type=float, default=-1.0,
                   help="Cross-gradient weight (default: auto = cell_size^4/16)")

    # Optimization
    p.add_argument("--max-iter", type=int, default=20)
    p.add_argument("--max-irls", type=int, default=10)

    # Clustering
    p.add_argument("--n-clusters", type=int, default=4)
    p.add_argument("--cluster-seed", type=int, default=42)

    # Decimation
    p.add_argument("--target-vertices", type=int, default=500,
                   help="Target vertices per contact surface after decimation")

    # Depth weighting
    p.add_argument("--depth-weight-beta-grav", type=float, default=2.0,
                   help="Depth weighting exponent for gravity (default: 2.0)")
    p.add_argument("--depth-weight-beta-mag", type=float, default=3.0,
                   help="Depth weighting exponent for magnetics (default: 3.0)")
    p.add_argument("--depth-weight-z0", type=float, default=-1.0,
                   help="Depth weighting reference depth (default: auto)")

    # Borehole geophysics
    p.add_argument("--borehole-properties-csv",
                   help="Path to downhole property measurements CSV")
    p.add_argument("--no-borehole-ref-model", action="store_true",
                   help="Disable building reference model from borehole properties")
    p.add_argument("--borehole-local-weight-radius", type=float, default=-1.0,
                   help="Radius to reduce regularization near boreholes (default: auto)")

    # Coordinates
    p.add_argument("--origin-x", type=float, default=0.0)
    p.add_argument("--origin-y", type=float, default=0.0)
    p.add_argument("--z-datum", type=float, default=0.0)

    # Output
    p.add_argument("--output-dir", default="lithoseed_output")
    p.add_argument("--no-plot", action="store_true", help="Disable diagnostic plots")

    # Borehole
    p.add_argument("--borehole-uncertainty", type=float, default=25.0,
                   help="Borehole intersection tolerance in metres")

    return p


def build_extract_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser for 'lithoseed extract' (no re-inversion)."""

    p = argparse.ArgumentParser(
        description="LithoSeed extract — extract surfaces from existing inversion",
    )
    p.add_argument("--mesh", required=True, help="Path to UBC .msh mesh file")
    p.add_argument("--density", required=True, help="Path to density .mod or .xyz file")
    p.add_argument("--susceptibility", default=None,
                   help="Path to susceptibility .mod or .xyz file")
    p.add_argument("--n-clusters", type=int, default=4)
    p.add_argument("--cluster-seed", type=int, default=42)
    p.add_argument("--target-vertices", type=int, default=500)
    p.add_argument("--origin-x", type=float, default=0.0)
    p.add_argument("--origin-y", type=float, default=0.0)
    p.add_argument("--z-datum", type=float, default=0.0)
    p.add_argument("--output-dir", default="lithoseed_output")
    return p
