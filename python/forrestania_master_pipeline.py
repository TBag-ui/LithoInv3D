#!/usr/bin/env python3
"""
Forrestania Master Pipeline -- Presentation-Ready Outputs
=========================================================

Complete end-to-end pipeline from raw survey data to presentation-ready
SimPEG inversion meshes, clustered lithology volumes, and C++ INI configs.

Produces EVERYTHING needed for the Forrestania presentation:
  - SimPEG density + susceptibility models (GA-compatible .mod + .msh)
  - Absolute-density version for visualization (background 2.67 g/cm³)
  - Model-space AND data-coordinate (UTM) copies of all meshes
  - GMM intersection clustering (12 groups: 4 density x 3 susceptibility)
  - Closed volume extraction (voxel-face, zero-gap, shared vertices)
  - Fragment cleanup + disconnected-component splitting
  - Single-property clustering (density-only + susceptibility-only)
  - INI configs for C++ litho inversion (baseline + perturbed properties)
  - Scatter plots for QC

Usage
-----
    cd python
    python forrestania_master_pipeline.py

Requirements
------------
    SimPEG, discretize, numpy, scipy, scikit-learn, scikit-image, pandas, matplotlib
    Local packages: cluster_api, lithoseed (auto-added to path)

Output Structure
----------------
    ../apps/forrestania/presentation/
        simpeg/
            simpeg_mesh.msh                    # UBC-format tensor mesh
            simpeg_density.mod                  # Density contrast (GA format)
            simpeg_susceptibility.mod           # Susceptibility (GA format)
            simpeg_density_absolute.mod         # Absolute density = contrast + 2.67
            data_coords/                        # UTM-coordinate copies
                simpeg_mesh.msh
                simpeg_density_absolute.mod
                simpeg_susceptibility.mod
        volumes/
            combination_clusters/               # GMM intersection (density x susc)
            density_clusters/                   # Density-only GMM volumes
            susc_clusters/                      # Susceptibility-only GMM volumes
            user_clusters/                      # Manual lithology-gated volumes
        configs/
            baseline/
                resolved_config.ini             # C++ inversion config (baseline properties)
            perturbed/
                resolved_config.ini             # C++ inversion config (perturbed properties)
        figures/
            cluster_scatter_intersection.png
            cluster_scatter_density.png
            cluster_scatter_susc.png
            data_fit.png
        cluster_properties.csv                  # Master cluster table
        observed_gravity.csv                    # Observed data (local coords)
        observed_magnetic.csv

Author: Generated 2026-06-05
"""

import csv
import os
import shutil
import sys
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import griddata, LinearNDInterpolator

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT / "python" / "cluster_api" / "src"))
sys.path.insert(0, str(PROJ_ROOT / "python" / "lithoseed" / "src"))

DATASETS = PROJ_ROOT / "datasets" / "Forrestania"
FORRESTANIA_APP = PROJ_ROOT / "apps" / "forrestania"
PRESENTATION = FORRESTANIA_APP / "presentation"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
IGRF_F = 58874.0       # nT -- IGRF total field at Forrestania (epoch 2026.0)
IGRF_I = -66.2         # deg -- inclination (negative = upward, southern hemisphere)
IGRF_D = -0.1          # deg -- declination
BACKGROUND_DENSITY = 2.67  # g/cm³ -- background density for absolute conversion
CELL_SIZE = 200.0      # m -- tensor mesh cell size
DEPTH_CORE = 1000.0    # m -- depth below minimum topography
MAX_ITER = 100         # SimPEG inversion max iterations
GRAV_UNCERTAINTY = 0.005  # mGal (tight — forces more iterations)
MAG_UNCERTAINTY = 5.0     # nT (tight — forces more iterations)
N_DENSITY_CLUSTERS = 4   # density bins for intersection clustering
N_SUSC_CLUSTERS = 3      # susceptibility bins for intersection clustering
SENTINEL = -9999.0       # Geoscience Analyst "no data" sentinel
CLEAN_MIN_CELLS = 25     # minimum cells per fragment before reassignment
RNG = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Utility: GA-compatible .mod writer
# ---------------------------------------------------------------------------

def write_mod_ga(path, values_flat, nx, ny, nz, air_mask=None, order='F'):
    """
    Write a UBC-GIF .mod file compatible with Geoscience Analyst.

    UBC-GIF format specification:
      - NO header line (GA reads geometry from the .msh file)
      - One value per line
      - Cell ordering: z fastest (top-to-bottom), then x (easting),
        then y (northing) slowest.
      - Cell [1,1,1] = top, south-west corner of the model.
      - Air cells set to SENTINEL (-9999) for transparent rendering.

    Parameters
    ----------
    path : str or Path
        Output file path.
    values_flat : np.ndarray
        1D array of shape (nx*ny*nz,) in SimPEG Fortran order
        (x fastest, then y, then z).
    nx, ny, nz : int
        Mesh dimensions (easting, northing, vertical).
    air_mask : np.ndarray or None
        (nx, ny, nz) boolean mask in same order as values_flat.
    order : str
        'F' = Fortran (SimPEG native: x fastest, then y, then z).
    """
    # Reshape to 3D in SimPEG Fortran order: data[i_x, j_y, k_z]
    # SimPEG stores z BOTTOM-TO-TOP (k=0 is deepest).
    data_3d = values_flat.reshape(nx, ny, nz, order=order)
    if air_mask is not None:
        data_3d[air_mask] = SENTINEL
    # UBC-GIF order: z TOP-TO-BOTTOM fastest, then x, then y.
    # Flip z-axis (::−1) so k=0 (SimPEG bottom) writes LAST,
    # and k=nz−1 (SimPEG top) writes FIRST.
    data_ubc = np.transpose(data_3d[:, :, ::-1], (1, 0, 2))  # (ny, nx, nz)
    data_flat = data_ubc.ravel(order='C')  # z top-down fastest, then x, then y
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, 'w') as f:
        for v in data_flat:
            f.write(f"{v:.6e}\n")

    n_masked = int(air_mask.sum()) if air_mask is not None else 0
    print(f"  Wrote {path} ({len(data_flat)} cells, {n_masked} air-masked)")


# ---------------------------------------------------------------------------
# Utility: UBC mesh writer (with optional origin override for data coords)
# ---------------------------------------------------------------------------

def write_mesh_ubc(path, x0, y0, z0, dx, dy, dz, nx, ny, nz):
    """
    Write a UBC-format tensor mesh file (.msh).

    Format::

        nx ny nz
        x0 y0 z0
        hx[0] hx[1] ... hx[nx-1]
        hy[0] hy[1] ... hy[ny-1]
        hz[0] hz[1] ... hz[nz-1]

    Origin (x0, y0, z0) is the southwest top corner.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, 'w') as f:
        f.write(f"{nx} {ny} {nz}\n")
        f.write(f"{x0:.6f} {y0:.6f} {z0:.6f}\n")
        f.write(" ".join([f"{dx:.6f}"] * nx) + "\n")
        f.write(" ".join([f"{dy:.6f}"] * ny) + "\n")
        f.write(" ".join([f"{dz:.6f}"] * nz) + "\n")
    print(f"  Wrote {path}")


# ---------------------------------------------------------------------------
# Utility: convergence CSV writer
# ---------------------------------------------------------------------------

def _write_conv_csv(history, path):
    """Write a SimPEG convergence CSV."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, 'w') as f:
        f.write("iter,beta,phi_d,phi_m,f\n")
        for h in history:
            f.write(f"{h['iter']},{h['beta']:.6e},{h['phi_d']:.6e},"
                    f"{h['phi_m']:.6e},{h['f']:.6e}\n")


def _make_convergence_recorder():
    """Return a ConvergenceRecorder class (must call after SimPEG import)."""
    from simpeg.directives import InversionDirective

    class _Recorder(InversionDirective):
        def __init__(self):
            super().__init__()
            self.history = []

        def initialize(self):
            pass

        def endIter(self):
            self.history.append({
                "iter": int(self.opt.iter),
                "beta": float(self.invProb.beta),
                "phi_d": float(self.invProb.phi_d),
                "phi_m": float(self.invProb.phi_m),
                "f": float(self.invProb.phi_d + self.invProb.beta * self.invProb.phi_m),
            })

        def write_csv(self, path):
            _write_conv_csv(self.history, path)

    return _Recorder


# ===================================================================
# STEP 1 -- Data Preparation
# ===================================================================

def step1_prepare_data():
    """
    Load and preprocess gravity and magnetic survey data.

    Returns
    -------
    grav_xyz : np.ndarray (N, 3)
        Gravity station locations in UTM coordinates (x=easting, y=northing, z=elevation).
    g_obs : np.ndarray (N,)
        Observed gravity (mGal), sign-flipped for SimPEG convention.
    mag_xyz : np.ndarray (M, 3)
        Magnetic station locations (UTM).
    t_obs : np.ndarray (M,)
        Detrended RMI (nT).
    topo : np.ndarray (P, 3)
        DEM points built from station elevations.
    local_origin : tuple (ox, oy)
        Survey center for local-coordinate transform.
    z_datum : float
        Mean survey elevation.
    """
    print("=" * 70)
    print("STEP 1 -- DATA PREPARATION")
    print("=" * 70)

    # -- Gravity --------------------------------------------------------
    grav_csv = DATASETS / "Forrestania_Gravity_Station_trim_.csv"
    if not grav_csv.is_file():
        raise FileNotFoundError(f"Gravity data not found: {grav_csv}")

    df_grav = pd.read_csv(grav_csv)
    grav_xyz = df_grav[["X", "Y", "Z"]].values.astype(np.float64)
    # SimPEG convention: positive gz is downward. Flip sign.
    g_obs = -df_grav["FGrav_mgal"].values.astype(np.float64)
    print(f"  Gravity: {len(df_grav)} stations, "
          f"g_obs range [{g_obs.min():.4f}, {g_obs.max():.4f}] mGal")

    # -- Magnetics ------------------------------------------------------
    mag_csv = DATASETS / "60472_AOI4_regional_mag_Zfixed.csv"
    if not mag_csv.is_file():
        raise FileNotFoundError(f"Magnetic data not found: {mag_csv}")

    df_mag = pd.read_csv(mag_csv)
    df_mag["RMI"] = df_mag["MAGCOMP"] - df_mag["IGRF"]

    # Downsample to ~500 stations for computational efficiency
    n_target = min(500, len(df_mag))
    idx = RNG.choice(len(df_mag), n_target, replace=False)
    sub = df_mag.iloc[idx].copy()
    mag_xyz = sub[["X", "Y", "Z"]].values.astype(np.float64)
    rmi = sub["RMI"].values.astype(np.float64)

    # Planar detrend
    A = np.c_[mag_xyz[:, 0], mag_xyz[:, 1], np.ones(len(mag_xyz))]
    coeff, _, _, _ = np.linalg.lstsq(A, rmi, rcond=None)
    t_obs = rmi - A @ coeff
    print(f"  Magnetics: {len(df_mag)} raw -> {len(sub)} downsampled")
    print(f"  Detrend: slope_x={coeff[0]:.6e}, slope_y={coeff[1]:.6e}, "
          f"intercept={coeff[2]:.1f} nT")
    print(f"  RMI_detrended range [{t_obs.min():.1f}, {t_obs.max():.1f}] nT")

    # -- DEM (from station elevations) ----------------------------------
    spacing = 120.0
    pad = 500.0
    all_pts = grav_xyz  # use gravity station elevations
    x_min, x_max = all_pts[:, 0].min() - pad, all_pts[:, 0].max() + pad
    y_min, y_max = all_pts[:, 1].min() - pad, all_pts[:, 1].max() + pad
    dem_x = np.arange(x_min, x_max + spacing, spacing)
    dem_y = np.arange(y_min, y_max + spacing, spacing)
    dem_xx, dem_yy = np.meshgrid(dem_x, dem_y)
    dem_zz = griddata(all_pts[:, :2], all_pts[:, 2],
                       (dem_xx, dem_yy), method="linear")
    mask_nan = np.isnan(dem_zz)
    if mask_nan.any():
        dem_zz_nn = griddata(all_pts[:, :2], all_pts[:, 2],
                              (dem_xx, dem_yy), method="nearest")
        dem_zz[mask_nan] = dem_zz_nn[mask_nan]
    topo = np.c_[dem_xx.ravel(), dem_yy.ravel(), dem_zz.ravel()]
    print(f"  DEM: {len(topo)} points, {spacing}m spacing, "
          f"Z=[{topo[:, 2].min():.1f}, {topo[:, 2].max():.1f}] m")

    # -- Coordinate system parameters -----------------------------------
    ox = (grav_xyz[:, 0].min() + grav_xyz[:, 0].max()) / 2.0
    oy = (grav_xyz[:, 1].min() + grav_xyz[:, 1].max()) / 2.0
    z_datum = float(grav_xyz[:, 2].mean())
    local_origin = (ox, oy)
    print(f"  Local origin: ({ox:.1f}, {oy:.1f}), z_datum={z_datum:.1f} m")

    return grav_xyz, g_obs, mag_xyz, t_obs, topo, local_origin, z_datum


# ===================================================================
# STEP 2 -- Build Tensor Mesh (FIXED: no overhead padding)
# ===================================================================

def step2_build_mesh(topo, cell_size=CELL_SIZE, depth_core=DEPTH_CORE):
    """
    Build a 3D tensor mesh for SimPEG inversions.

    **FIX applied 2026-06-05**: The original ``build_tensor_mesh`` added
    ``+ 2*cell_size`` to the top elevation, which placed 2 full cell layers
    entirely above the ground surface.  We now snap the top to the nearest
    cell boundary AT OR ABOVE the maximum topography -- eliminating the
    overhead padding while keeping the top cell at or above ground.

    Returns
    -------
    mesh : discretize.TensorMesh
    active : np.ndarray (bool)
        Active-cells mask (cells below topography).
    n_act : int
        Number of active cells.
    nx, ny, nz : int
        Mesh dimensions.
    dx, dy, dz : float
        Cell sizes.
    x0, y0, z0 : float
        Mesh origin (southwest top corner, UTM).
    """
    from discretize import TensorMesh
    from discretize.utils import active_from_xyz

    print("\n" + "=" * 70)
    print("STEP 2 -- BUILD TENSOR MESH (FIXED: no overhead padding)")
    print("=" * 70)

    xyz_all = topo
    x_min, x_max = xyz_all[:, 0].min(), xyz_all[:, 0].max()
    y_min, y_max = xyz_all[:, 1].min(), xyz_all[:, 1].max()
    z_min = xyz_all[:, 2].min() - depth_core
    z_max = xyz_all[:, 2].max()

    # Pad horizontally by 2 cells on each side
    x0 = np.floor(x_min / cell_size) * cell_size - 2 * cell_size
    y0 = np.floor(y_min / cell_size) * cell_size - 2 * cell_size

    # FIX: Snap z_top to the cell boundary at or just above max topography.
    # Old code: z_top = ceil(z_max/cell_size)*cell_size + 2*cell_size
    #           -> 800m with z_max=380m, cell_size=200m (2 layers of air!)
    # New code: just snap to the next cell boundary above ground.
    z_top = np.ceil(z_max / cell_size) * cell_size

    # z_bottom: at least depth_core below min topography, snapped to cell grid
    z_bot = np.floor(z_min / cell_size) * cell_size - 2 * cell_size

    nx = int(np.ceil((x_max + 2 * cell_size - x0) / cell_size))
    ny = int(np.ceil((y_max + 2 * cell_size - y0) / cell_size))
    nz = int(np.ceil((z_top - z_bot) / cell_size))

    hx = np.ones(nx) * cell_size
    hy = np.ones(ny) * cell_size
    hz = np.ones(nz) * cell_size
    origin = np.array([x0, y0, z_bot])

    mesh = TensorMesh([hx, hy, hz], origin=origin)

    # Active cells: those below topography
    active = active_from_xyz(mesh, topo)
    n_act = int(active.sum())

    # Calculate cell centers for diagnostics
    z_centers = mesh.cell_centers_z[::nx*ny]  # one per z-layer
    topo_max = topo[:, 2].max()

    print(f"  Mesh: {nx}x{ny}x{nz} = {mesh.n_cells} cells @ {cell_size}m")
    print(f"  Origin (UTM): ({x0:.1f}, {y0:.1f}, {z_bot:.1f})")
    print(f"  Z range: [{z_bot:.1f}, {z_top:.1f}] m")
    print(f"  Topo max: {topo_max:.1f} m")
    print(f"  Z cell centers: {np.array2string(z_centers, precision=0)}")
    print(f"  Air cells: {int((~active).sum())} / {mesh.n_cells}")
    print(f"  Active cells: {n_act} / {mesh.n_cells}")
    print(f"  FIX VERIFIED: top cell center at {z_centers[0]:.0f}m, "
          f"ground at ~{topo_max:.0f}m -> "
          f"{'OK' if z_centers[0] <= topo_max + cell_size else 'STILL TOO HIGH!'}")

    return mesh, active, n_act, nx, ny, nz, cell_size, cell_size, cell_size, x0, y0, z_bot


# ===================================================================
# STEP 3 -- SimPEG Gravity Inversion
# ===================================================================

def step3_gravity_inversion(mesh, active, n_act, nx, ny, nz,
                             grav_xyz, g_obs, output_dir):
    """
    Run SimPEG gravity-only inversion recovering density contrast.

    Matches the proven configuration from run_independent_inversions.py:
      - IdentityMap on full mesh (all cells)
      - Sparse regularization (L2 smallness + L1 smoothness) with IRLS
      - Blocky, focused features via norms=[2,1,1,1]
      - NO active_cells — regularizer sees full mesh

    Returns
    -------
    rec_density : np.ndarray (nc,)
        Recovered density contrast on ALL mesh cells.
    """
    from simpeg import maps, data_misfit, regularization, optimization
    from simpeg import inverse_problem, inversion, directives, data
    from simpeg.potential_fields import gravity

    ConvergenceRecorder = _make_convergence_recorder()
    nc = mesh.n_cells

    print("\n" + "=" * 70)
    print("STEP 3 -- SimPEG GRAVITY INVERSION")
    print("=" * 70)

    grav_rx = gravity.receivers.Point(grav_xyz, components=["gz"])
    grav_src = gravity.sources.SourceField(receiver_list=[grav_rx])
    grav_survey = gravity.survey.Survey(grav_src)
    grav_sim = gravity.simulation.Simulation3DIntegral(
        mesh, survey=grav_survey, rhoMap=maps.IdentityMap(nP=nc),
        store_sensitivities="ram",
    )

    grav_data = data.Data(
        grav_survey, dobs=g_obs,
        standard_deviation=np.ones(len(g_obs)) * GRAV_UNCERTAINTY,
    )
    dmis = data_misfit.L2DataMisfit(data=grav_data, simulation=grav_sim)

    # Sparse with L2 smallness + L1 smoothness (blocky edges)
    # alpha_z=5.0 penalizes vertical oscillations — prevents checkerboard
    # at depth where sensitivity is very low
    reg = regularization.Sparse(
        mesh=mesh, reference_model=np.zeros(nc),
        norms=[2, 1, 1, 1],
        alpha_s=1.0, alpha_x=1.0, alpha_y=1.0, alpha_z=5.0)

    m0 = np.ones(nc) * 1e-4

    opt = optimization.ProjectedGNCG(
        maxIter=MAX_ITER,
        lower=np.ones(nc) * -0.8, upper=np.ones(nc) * 0.8,
        maxIterLS=20)

    inv_prob = inverse_problem.BaseInvProblem(dmis, reg, opt)
    conv_recorder = ConvergenceRecorder()
    dir_list = [
        directives.Update_IRLS(max_irls_iterations=10, f_min_change=1e-3),
        directives.BetaEstimate_ByEig(beta0_ratio=1.0),
        conv_recorder,
    ]
    inv = inversion.BaseInversion(inv_prob, directiveList=dir_list)

    print(f"  Running gravity inversion (max_iter={MAX_ITER}, sparse+IRLS, alpha_z=5)...")
    t0 = time.time()
    rec_density = inv.run(m0)
    elapsed = time.time() - t0
    pred = np.asarray(inv_prob.dpred).ravel()
    rms = np.sqrt(np.mean((g_obs - pred)**2))
    print(f"  Completed in {elapsed:.1f}s  RMS={rms:.4f}")
    print(f"  Density contrast range: [{rec_density.min():.4f}, "
          f"{rec_density.max():.4f}] g/cm³")

    # Save convergence
    os.makedirs(output_dir, exist_ok=True)
    conv_recorder.write_csv(os.path.join(output_dir, "simpeg_gravity_convergence.csv"))

    return rec_density


# ===================================================================
# STEP 4 -- SimPEG Magnetic Inversion
# ===================================================================

def step4_magnetic_inversion(mesh, active, n_act, nx, ny, nz,
                              mag_xyz, t_obs, output_dir):
    """
    Run SimPEG magnetic-only inversion recovering susceptibility (SI).

    Matches the proven configuration from run_independent_inversions.py:
      - IdentityMap on full mesh
      - Sparse regularization (L2 smallness + L1 smoothness) with IRLS
      - Blocky, focused features via norms=[2,1,1,1]

    Returns
    -------
    rec_susc : np.ndarray (nc,)
        Recovered susceptibility on ALL mesh cells.
    """
    from simpeg import maps, data_misfit, regularization, optimization
    from simpeg import inverse_problem, inversion, directives, data
    from simpeg.potential_fields import magnetics

    ConvergenceRecorder = _make_convergence_recorder()
    nc = mesh.n_cells

    print("\n" + "=" * 70)
    print("STEP 4 -- SimPEG MAGNETIC INVERSION")
    print("=" * 70)

    mag_rx = magnetics.receivers.Point(mag_xyz, components=["tmi"])
    mag_src = magnetics.sources.UniformBackgroundField(
        receiver_list=[mag_rx],
        amplitude=IGRF_F,
        inclination=IGRF_I,
        declination=IGRF_D,
    )
    mag_survey = magnetics.survey.Survey(mag_src)
    mag_sim = magnetics.simulation.Simulation3DIntegral(
        mesh, survey=mag_survey, chiMap=maps.IdentityMap(nP=nc),
        store_sensitivities="ram",
    )

    mag_data = data.Data(
        mag_survey, dobs=t_obs,
        standard_deviation=np.ones(len(t_obs)) * MAG_UNCERTAINTY,
    )
    dmis = data_misfit.L2DataMisfit(data=mag_data, simulation=mag_sim)

    # Sparse with L2 smallness + L1 smoothness (blocky edges)
    # alpha_z=5.0 prevents checkerboard at depth
    reg = regularization.Sparse(
        mesh=mesh, reference_model=np.zeros(nc),
        norms=[2, 1, 1, 1],
        alpha_s=1.0, alpha_x=1.0, alpha_y=1.0, alpha_z=5.0)

    m0 = np.ones(nc) * 1e-4

    opt = optimization.ProjectedGNCG(
        maxIter=MAX_ITER,
        lower=np.ones(nc) * 0.0, upper=np.ones(nc) * 0.5,
        maxIterLS=20)

    inv_prob = inverse_problem.BaseInvProblem(dmis, reg, opt)
    conv_recorder = ConvergenceRecorder()
    dir_list = [
        directives.Update_IRLS(max_irls_iterations=10, f_min_change=1e-3),
        directives.BetaEstimate_ByEig(beta0_ratio=1.0),
        conv_recorder,
    ]
    inv = inversion.BaseInversion(inv_prob, directiveList=dir_list)

    print(f"  Running magnetic inversion (max_iter={MAX_ITER}, sparse+IRLS)...")
    t0 = time.time()
    rec_susc = inv.run(m0)
    elapsed = time.time() - t0
    pred = np.asarray(inv_prob.dpred).ravel()
    rms = np.sqrt(np.mean((t_obs - pred)**2))
    print(f"  Completed in {elapsed:.1f}s  RMS={rms:.4f}")
    print(f"  Susceptibility range: [{rec_susc.min():.6f}, "
          f"{rec_susc.max():.6f}] SI")

    # Save convergence
    os.makedirs(output_dir, exist_ok=True)
    conv_recorder.write_csv(os.path.join(output_dir, "simpeg_magnetic_convergence.csv"))

    return rec_susc


# ===================================================================
# STEP 5 -- Convert to Absolute Density
# ===================================================================

def step5_absolute_density(rec_density, background=BACKGROUND_DENSITY):
    """
    Convert density contrast to absolute density.

    Absolute density = recovered contrast + background density (2.67 g/cm³).

    This is for PRESENTATION -- the GA mesh should display physical density
    values (e.g., 2.3--3.2 g/cm³) rather than contrasts (e.g., -0.4--+0.5).
    The clustering pipeline still uses the contrast values.

    Returns
    -------
    rec_density_abs : np.ndarray (n_act,)
        Absolute density (g/cm³).
    """
    print("\n" + "=" * 70)
    print("STEP 5 -- CONVERT TO ABSOLUTE DENSITY")
    print("=" * 70)

    rec_density_abs = rec_density + background
    print(f"  Background density: {background:.2f} g/cm³")
    print(f"  Contrast range:    [{rec_density.min():.4f}, {rec_density.max():.4f}]")
    print(f"  Absolute range:    [{rec_density_abs.min():.4f}, {rec_density_abs.max():.4f}]")

    return rec_density_abs


# ===================================================================
# STEP 6 -- Build Air Mask + Export GA-Compatible Files (model space)
# ===================================================================

def step6_export_ga_model_space(mesh, active, rec_density, rec_density_abs,
                                 rec_susc, nx, ny, nz, x0, y0, z_bot,
                                 dx, dy, dz, output_dir):
    """
    Export SimPEG outputs in Geoscience Analyst-compatible format.

    Produces:
      - simpeg_mesh.msh            -- UBC-format tensor mesh
      - simpeg_density.mod         -- Density contrast (GA format: no header)
      - simpeg_density_absolute.mod -- Absolute density (GA format)
      - simpeg_susceptibility.mod  -- Susceptibility (GA format)

    All .mod files are written WITHOUT the ``nx ny nz`` header that UBC
    format normally prepends -- GA reads geometry from the .msh file.
    Air cells are set to -9999 (GA's "no data" sentinel) for transparent
    rendering above topography.
    """
    from discretize import TensorMesh

    print("\n" + "=" * 70)
    print("STEP 6 -- GA-COMPATIBLE EXPORT (model space)")
    print("=" * 70)

    os.makedirs(output_dir, exist_ok=True)

    # Build air mask: True where cell center is above topography
    xc = mesh.cell_centers_x[::nx*ny] if hasattr(mesh, 'cell_centers_x') else \
         np.linspace(x0 + dx/2, x0 + (nx - 0.5)*dx, nx)
    # Actually, use the full cell centers
    air_mask_3d = (~active).reshape(nx, ny, nz, order='F')

    n_cells_total = nx * ny * nz
    n_air = int((~active).sum())

    # -- Mesh file (always .msh for GA) ----------------------------------
    msh_path = os.path.join(output_dir, "simpeg_mesh.msh")
    write_mesh_ubc(msh_path, x0, y0, z_bot + nz * dz, dx, dy, dz, nx, ny, nz)

    # -- Density contrast (GA format) ------------------------------------
    # Expand active-only values to full mesh (inactive = 0, will be masked)
    full_density = np.zeros(n_cells_total)
    full_density[active] = rec_density

    mod_path = os.path.join(output_dir, "simpeg_density.mod")
    write_mod_ga(mod_path, full_density, nx, ny, nz,
                 air_mask=air_mask_3d, order='F')

    # -- Absolute density (GA format) -- for presentation -----------------
    full_density_abs = np.zeros(n_cells_total)
    full_density_abs[active] = rec_density_abs

    mod_abs_path = os.path.join(output_dir, "simpeg_density_absolute.mod")
    write_mod_ga(mod_abs_path, full_density_abs, nx, ny, nz,
                 air_mask=air_mask_3d.copy(), order='F')

    # -- Susceptibility (GA format) --------------------------------------
    full_susc = np.zeros(n_cells_total)
    full_susc[active] = rec_susc

    mod_susc_path = os.path.join(output_dir, "simpeg_susceptibility.mod")
    write_mod_ga(mod_susc_path, full_susc, nx, ny, nz,
                 air_mask=air_mask_3d.copy(), order='F')

    print(f"\n  GA-compatible files written to {output_dir}/")
    print(f"  Load in Geoscience Analyst:")
    print(f"    1. Open mesh:  simpeg_mesh.msh")
    print(f"    2. Load model: simpeg_density_absolute.mod (absolute density)")
    print(f"    3. Load model: simpeg_susceptibility.mod (susceptibility)")

    return air_mask_3d


# ===================================================================
# STEP 7 -- Data-Coordinate (UTM) Copies
# ===================================================================

def step7_export_data_coords(mesh, rec_density_abs, rec_susc, active,
                              nx, ny, nz, x0, y0, z_bot, dx, dy, dz,
                              output_dir):
    """
    Export copies of all meshes in the original data coordinate system (UTM).

    The main exports in step 6 are already in UTM coordinates -- the SimPEG
    mesh origin IS the UTM origin.  This step creates an explicit
    ``data_coords/`` subdirectory with copies for presentation clarity,
    and also writes a local-coordinate version (model space) by subtracting
    the local origin.

    Model space = UTM coordinates translated so survey center is (0, 0, z_datum).
    Data space = original UTM MGA Zone 50 coordinates.
    """
    print("\n" + "=" * 70)
    print("STEP 7 -- COORDINATE-SYSTEM COPIES")
    print("=" * 70)

    data_coords_dir = os.path.join(output_dir, "data_coords")
    model_space_dir = os.path.join(output_dir, "model_space")
    os.makedirs(data_coords_dir, exist_ok=True)
    os.makedirs(model_space_dir, exist_ok=True)

    # The main outputs are already in UTM (data coordinates).
    # Just copy them to the explicit directory.
    z_top = z_bot + nz * dz

    # -- Data coordinates (UTM) ------------------------------------------
    write_mesh_ubc(os.path.join(data_coords_dir, "simpeg_mesh.msh"),
                   x0, y0, z_top, dx, dy, dz, nx, ny, nz)

    air_mask_3d = (~active).reshape(nx, ny, nz, order='F')

    full_dens_abs = np.zeros(nx * ny * nz)
    full_dens_abs[active] = rec_density_abs
    write_mod_ga(os.path.join(data_coords_dir, "simpeg_density_absolute.mod"),
                 full_dens_abs, nx, ny, nz, air_mask=air_mask_3d.copy(), order='F')

    full_susc = np.zeros(nx * ny * nz)
    full_susc[active] = rec_susc
    write_mod_ga(os.path.join(data_coords_dir, "simpeg_susceptibility.mod"),
                 full_susc, nx, ny, nz, air_mask=air_mask_3d.copy(), order='F')

    print(f"  Data-coordinate copies: {data_coords_dir}/")
    print(f"  These use UTM MGA Zone 50 coordinates directly.")
    print(f"  Model-space copies: {model_space_dir}/")
    print(f"  (local origin subtracted -- produced by clustering pipeline)")


# ===================================================================
# STEP 8 -- MeshData Conversion
# ===================================================================

def step8_convert_to_meshdata(mesh, active, rec_density, n_act):
    """
    Convert SimPEG TensorMesh + active cells to cluster_api MeshData.

    MeshData is the bridge between SimPEG's discretize.TensorMesh and the
    lithoseed clustering + volume extraction pipeline.
    """
    from cluster_api._io import MeshData

    print("\n" + "=" * 70)
    print("STEP 8 -- CONVERT TO MeshData")
    print("=" * 70)

    dx = float(mesh.h[0][0])
    dy = float(mesh.h[1][0])
    dz = float(mesh.h[2][0])
    nx = len(mesh.h[0])
    ny = len(mesh.h[1])
    nz = len(mesh.h[2])

    origin = mesh.origin
    x0 = float(origin[0])
    y0 = float(origin[1])
    # z0 is top of mesh (origin[2] is bottom)
    z0 = float(origin[2]) + nz * dz

    cc = mesh.cell_centers[active]
    x_center = cc[:, 0]
    y_center = cc[:, 1]
    z_center = cc[:, 2]

    # 1-based indices (cluster_api convention)
    ix = np.round((x_center - x0) / dx - 0.5).astype(int) + 1
    iy = np.round((y_center - y0) / dy - 0.5).astype(int) + 1
    iz_from_bottom = np.round((z_center - float(origin[2])) / dz - 0.5).astype(int)
    iz = nz - iz_from_bottom  # convert to top-down indexing

    mesh_data = MeshData(
        nx=nx, ny=ny, nz=nz,
        x0=x0, y0=y0, z0=z0,
        dx=dx, dy=dy, dz=dz,
        n_active=n_act,
        ix=ix, iy=iy, iz=iz,
        x_center=x_center, y_center=y_center, z_center=z_center,
    )
    print(f"  MeshData: {mesh_data.nx}x{mesh_data.ny}x{mesh_data.nz}, "
          f"{mesh_data.n_active} active cells")
    return mesh_data


# ===================================================================
# STEP 9 -- Intersection Clustering + Contact Surfaces
# ===================================================================

def step9_clustering(mesh_data, rec_density, rec_susc, output_dir,
                      local_origin, z_datum, grav_xyz, g_obs, mag_xyz, t_obs):
    """
    Run GMM intersection clustering (1D density x 1D susceptibility).

    Configuration: N_DENSITY_CLUSTERS density bins x
                   N_SUSC_CLUSTERS susceptibility bins
                 -> up to N_DENSITY_CLUSTERS * N_SUSC_CLUSTERS intersection groups.

    This produces the 12-group starting model used for the C++ litho inversion.
    Contact surfaces between adjacent groups are extracted via marching cubes.
    """

    from cluster_api._cluster import cluster_intersection
    from lithoseed._pipeline import run_extract_from_labels

    print("\n" + "=" * 70)
    print("STEP 9 -- INTERSECTION CLUSTERING + CONTACT SURFACES")
    print("=" * 70)

    labels, summary = cluster_intersection(
        rec_density, rec_susc,
        n_density=N_DENSITY_CLUSTERS,
        n_susc=N_SUSC_CLUSTERS,
        random_state=42,
    )
    print(f"  Intersection clustering: {N_DENSITY_CLUSTERS} density x "
          f"{N_SUSC_CLUSTERS} susceptibility -> {summary.n_groups} non-empty groups")

    # Print group properties
    print(f"\n  {'Group':>6s}  {'Density mean':>14s}  {'Density std':>14s}  "
          f"{'Susc mean':>14s}  {'Susc std':>14s}  {'Count':>8s}")
    print(f"  {'-'*6}  {'-'*14}  {'-'*14}  {'-'*14}  {'-'*14}  {'-'*8}")
    for i in range(summary.n_groups):
        gid = summary.labels[i]
        d_mean = summary.density_mean[i]
        d_std = summary.density_std[i]
        s_mean = summary.susc_mean[i]
        s_std = summary.susc_std[i]
        count = summary.counts[i]
        print(f"  {gid:6d}  {d_mean:14.6f}  {d_std:14.6f}  "
              f"{s_mean:14.6e}  {s_std:14.6e}  {count:8d}")

    # Contact surface extraction (for INI config contacts)
    # Run extract_from_labels to get contact surfaces + INI
    result = run_extract_from_labels(
        mesh_data,
        rec_density,
        rec_susc,
        labels,
        summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=output_dir,
        export_formats=("ts", "obj"),
        gravity_xyz=grav_xyz,
        gravity_obs=g_obs,
        gravity_std=None,
        magnetic_xyz=mag_xyz,
        magnetic_obs=t_obs,
        magnetic_std=None,
        igrf_f=IGRF_F,
        igrf_i=IGRF_I,
        igrf_d=IGRF_D,
    )

    print(f"\n  Extracted {result.n_contacts} contact surface(s)")
    for cs in result.contacts:
        print(f"    Contact {cs.group_above}->{cs.group_below}: "
              f"{len(cs.vertices)} verts, depth={cs.median_depth:.0f}m")

    return labels, summary, result


# ===================================================================
# STEP 10 -- Volume Extraction (voxel-face, zero-gap)
# ===================================================================

def step10_volume_extraction(mesh_data, labels, summary, output_dir,
                              local_origin, z_datum):
    """
    Extract closed volumes via voxel-face boundary extraction.

    Every face in the 3D label grid is iterated once.  Shared quads between
    adjacent groups are assigned to both with opposite winding.  This
    guarantees zero-gap, shared-vertex, perfectly squared-off closed surfaces.

    Output: combination_clusters/ (cleaned, min {CLEAN_MIN_CELLS} cells, then split)
    """
    from lithoseed._pipeline import run_extract_group_volumes

    print("\n" + "=" * 70)
    print("STEP 10 -- VOLUME EXTRACTION (voxel-face, zero-gap)")
    print("=" * 70)

    # Cleaned: reassign small fragments
    vol_clean_dir = os.path.join(output_dir, "combination_clusters")
    print(f"\n  Cleaned ({CLEAN_MIN_CELLS} cell min) -> combination_clusters/")
    clean_result = run_extract_group_volumes(
        mesh_data, labels, summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=vol_clean_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=CLEAN_MIN_CELLS,
    )

    return clean_result


# ===================================================================
# STEP 11 -- All QC Plots (convergence, scatter, spatial maps)
# ===================================================================

def step11_all_plots(rec_density, rec_susc, labels, summary,
                      mesh, active, nx, ny, nz,
                      grav_xyz, g_obs, mag_xyz, t_obs,
                      simpeg_dir, output_dir):
    """
    Generate all QC plots for the presentation:

    1. Convergence plots -- iteration vs misfit, beta, model norm
       for both gravity and magnetic SimPEG inversions.
    2. Intersection scatter -- density vs susceptibility colored by cluster.
    3. XY spatial maps -- map-view slices of density and susceptibility
       at representative depths, with cluster boundaries.
    4. Single-property scatters -- histogram + cell-index plots.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    print("\n" + "=" * 70)
    print("STEP 11 -- QC PLOTS (convergence, scatter, spatial maps)")
    print("=" * 70)

    os.makedirs(output_dir, exist_ok=True)

    valid = labels >= 0
    d_val = rec_density[valid]
    s_val = rec_susc[valid]
    l_val = labels[valid]

    unique_labels = sorted(set(l_val))
    n_groups = len(unique_labels)
    cmap = plt.cm.tab10
    colors = [cmap(i % 10) for i in range(n_groups)]

    # -- 11a. Convergence plots ------------------------------------------
    _plot_convergence(simpeg_dir, output_dir)

    # -- 11b. Intersection scatter (density vs susceptibility) ------------
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, label in enumerate(unique_labels):
        mask = l_val == label
        ax.scatter(d_val[mask], s_val[mask], c=[colors[i]], alpha=0.5,
                   s=1, label=f"Group {label}")
    if hasattr(summary, 'density_mean') and hasattr(summary, 'susc_mean'):
        for i in range(len(summary.density_mean)):
            ax.scatter(summary.density_mean[i], summary.susc_mean[i],
                       c="black", marker="x", s=100, linewidths=2, zorder=5)
    ax.set_xlabel("Density contrast (g/cm³)")
    ax.set_ylabel("Magnetic Susceptibility (SI)")
    ax.set_title(f"GMM Intersection: {N_DENSITY_CLUSTERS}x{N_SUSC_CLUSTERS} -> {n_groups} groups")
    ax.legend(markerscale=5, fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = os.path.join(output_dir, "cluster_scatter_intersection.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Wrote {path}")

    # -- 11c. XY spatial maps --------------------------------------------
    _plot_spatial_maps(mesh, active, rec_density, rec_susc, labels,
                        nx, ny, nz, summary, output_dir)

    # -- 11d. Single-property scatters ------------------------------------
    def _single_scatter(vals, lbls, path, xlabel):
        v = vals[valid]
        l = lbls[valid]
        ul = sorted(set(l))
        cls = [cmap(i % 10) for i in range(len(ul))]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        for i, label in enumerate(ul):
            mask = l == label
            ax1.hist(v[mask], bins=50, alpha=0.6, color=cls[i], label=f"Group {label}")
        ax1.set_xlabel(xlabel); ax1.set_ylabel("Count")
        ax1.set_title("Per-Cluster Histogram"); ax1.legend(fontsize=7); ax1.grid(True, alpha=0.3)
        indices = np.arange(len(v))
        for i, label in enumerate(ul):
            mask = l == label
            ax2.scatter(indices[mask], v[mask], c=[cls[i]], alpha=0.3, s=1, label=f"Group {label}")
        ax2.set_xlabel("Cell Index"); ax2.set_ylabel(xlabel)
        ax2.set_title("Per-Cluster Scatter (by cell index)")
        ax2.legend(markerscale=5, fontsize=7); ax2.grid(True, alpha=0.3)
        fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
        print(f"  Wrote {path}")

    from cluster_api._cluster import cluster_lithology
    dens_labels, _ = cluster_lithology(
        rec_density, np.zeros_like(rec_density),
        n_clusters=N_DENSITY_CLUSTERS, random_state=42)
    _single_scatter(rec_density, dens_labels,
                     os.path.join(output_dir, "cluster_scatter_density.png"),
                     xlabel="Density contrast (g/cm³)")

    susc_labels, _ = cluster_lithology(
        np.zeros_like(rec_susc), rec_susc,
        n_clusters=N_SUSC_CLUSTERS, random_state=42)
    _single_scatter(rec_susc, susc_labels,
                     os.path.join(output_dir, "cluster_scatter_susc.png"),
                     xlabel="Magnetic Susceptibility (SI)")


def _plot_convergence(simpeg_dir, output_dir):
    """Plot convergence curves for gravity and magnetic inversions."""
    import matplotlib.pyplot as plt

    for name, csv_name in [("gravity", "simpeg_gravity_convergence.csv"),
                            ("magnetic", "simpeg_magnetic_convergence.csv")]:
        csv_path = os.path.join(simpeg_dir, csv_name)
        if not os.path.isfile(csv_path):
            print(f"  SKIP convergence plot: {csv_name} not found")
            continue

        try:
            data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[0] < 2:
                print(f"  SKIP convergence plot: {csv_name} has only {data.shape[0]} iterations")
                continue
        except Exception as e:
            print(f"  WARNING: Could not read {csv_name}: {e}")
            continue

        iters = data[:, 0].astype(int)
        beta  = data[:, 1]
        phi_d = data[:, 2]
        phi_m = data[:, 3]
        f     = data[:, 4]

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(f"SimPEG {name.capitalize()} Inversion Convergence", fontsize=14)

        # Data misfit
        axes[0, 0].semilogy(iters, phi_d, 'b.-', markersize=4)
        axes[0, 0].set_xlabel("Iteration"); axes[0, 0].set_ylabel("phi_d (data misfit)")
        axes[0, 0].grid(True, alpha=0.3)
        axes[0, 0].set_title(f"Final phi_d = {phi_d[-1]:.4e}")

        # Model norm
        axes[0, 1].plot(iters, phi_m, 'r.-', markersize=4)
        axes[0, 1].set_xlabel("Iteration"); axes[0, 1].set_ylabel("phi_m (model norm)")
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title(f"Final phi_m = {phi_m[-1]:.4e}")

        # Total objective
        axes[1, 0].semilogy(iters, f, 'g.-', markersize=4)
        axes[1, 0].set_xlabel("Iteration"); axes[1, 0].set_ylabel("Phi = phi_d + beta·phi_m")
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].set_title(f"Final Phi = {f[-1]:.4e}")

        # Beta schedule
        axes[1, 1].semilogy(iters, beta, 'k.-', markersize=4)
        axes[1, 1].set_xlabel("Iteration"); axes[1, 1].set_ylabel("beta (trade-off)")
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_title(f"Final beta = {beta[-1]:.4e}")

        fig.tight_layout()
        path = os.path.join(output_dir, f"convergence_{name}.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  Wrote {path}")


def _plot_spatial_maps(mesh, active, rec_density, rec_susc, labels,
                        nx, ny, nz, summary, output_dir):
    """
    XY map-view slices of recovered density and susceptibility at key depths.

    Produces:
      - density_map.png       -- depth slice through peak density contrast
      - susceptibility_map.png -- depth slice through peak susceptibility
    """
    import matplotlib.pyplot as plt

    # Build 3D grids (nx, ny, nz) from active-cell values
    d_grid = np.full((nx, ny, nz), np.nan)
    s_grid = np.full((nx, ny, nz), np.nan)
    l_grid = np.full((nx, ny, nz), -1)

    active_3d = active.reshape(nx, ny, nz, order='F')
    d_3d = np.zeros((nx, ny, nz))
    s_3d = np.zeros((nx, ny, nz))
    l_3d = np.full((nx, ny, nz), -1)

    d_flat = np.zeros(nx * ny * nz)
    s_flat = np.zeros(nx * ny * nz)
    l_flat = np.full(nx * ny * nz, -1)
    d_flat[active] = rec_density
    s_flat[active] = rec_susc
    l_flat[active] = labels
    d_3d = d_flat.reshape(nx, ny, nz, order='F')
    s_3d = s_flat.reshape(nx, ny, nz, order='F')
    l_3d = l_flat.reshape(nx, ny, nz, order='F')

    # Get cell center coordinates
    dx = float(mesh.h[0][0])
    dy = float(mesh.h[1][0])
    dz = float(mesh.h[2][0])
    x0 = float(mesh.origin[0])
    y0 = float(mesh.origin[1])
    z0 = float(mesh.origin[2])

    xc = np.linspace(x0 + dx/2, x0 + (nx - 0.5)*dx, nx)
    yc = np.linspace(y0 + dy/2, y0 + (ny - 0.5)*dy, ny)
    zc = np.linspace(z0 + dz/2, z0 + (nz - 0.5)*dz, nz)

    # Pick depth slice with most active cells (typically near-surface)
    active_per_layer = active_3d.sum(axis=(0, 1))  # sum over x,y for each z
    best_k = int(np.argmax(active_per_layer))
    print(f"  Spatial maps: using depth slice k={best_k} (z={zc[best_k]:.0f}m, "
          f"{active_per_layer[best_k]} active cells)")

    # Also pick a deeper slice
    deep_k = max(0, best_k + max(1, nz // 4))
    if deep_k >= nz:
        deep_k = nz - 1

    unique_labels = sorted(set(labels[labels >= 0]))
    n_groups = len(unique_labels)
    cmap = plt.cm.tab10
    cluster_colors = [cmap(i % 10) for i in range(n_groups)]

    for prop_name, grid, prop_label in [
        ("density", d_3d, "Density contrast (g/cm³)"),
        ("susceptibility", s_3d, "Magnetic Susceptibility (SI)"),
    ]:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle(f"{prop_name.capitalize()} -- XY Map Slices", fontsize=14)

        for ax_idx, (k, depth_label) in enumerate([
            (best_k, f"Near-surface (z={zc[best_k]:.0f}m)"),
            (deep_k, f"Deep (z={zc[deep_k]:.0f}m)"),
        ]):
            ax = axes[ax_idx]
            slice_data = grid[:, :, k].T  # transpose: x=columns, y=rows

            # Mask inactive cells
            slice_active = active_3d[:, :, k].T
            slice_masked = np.where(slice_active, slice_data, np.nan)

            im = ax.pcolormesh(xc, yc, slice_masked,
                                cmap='RdBu_r', shading='auto',
                                vmin=np.nanpercentile(slice_masked, 5),
                                vmax=np.nanpercentile(slice_masked, 95))
            ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
            ax.set_title(depth_label)
            ax.set_aspect('equal')
            plt.colorbar(im, ax=ax, label=prop_label, shrink=0.8)

        fig.tight_layout()
        path = os.path.join(output_dir, f"{prop_name}_map.png")
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"  Wrote {path}")

    # -- Cluster label map at best depth ----------------------------------
    fig, ax = plt.subplots(figsize=(10, 8))
    slice_l = l_3d[:, :, best_k].T.astype(float)
    slice_active_l = active_3d[:, :, best_k].T
    slice_l[~slice_active_l] = np.nan

    # Use discrete colormap for clusters
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cluster_cmap = ListedColormap(
        [cmap(i % 10) for i in range(max(unique_labels) + 1)])
    bounds = np.arange(-0.5, max(unique_labels) + 1.5, 1)
    norm = BoundaryNorm(bounds, cluster_cmap.N)

    im = ax.pcolormesh(xc, yc, slice_l, cmap=cluster_cmap, norm=norm,
                        shading='auto')
    ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
    ax.set_title(f"Cluster Labels at z={zc[best_k]:.0f}m")
    ax.set_aspect('equal')
    cbar = plt.colorbar(im, ax=ax, label="Cluster ID", ticks=unique_labels)
    fig.tight_layout()
    path = os.path.join(output_dir, "cluster_labels_map.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Wrote {path}")


# ===================================================================
# STEP 12 -- Split Disconnected Volumes
# ===================================================================

def step12_split_volumes(vol_clean_dir, output_dir):
    """
    Detect and split disconnected components in cleaned volume .ts meshes.

    A single cluster may contain multiple spatially disconnected volumes
    (e.g., two separate mafic lenses at the same stratigraphic level).
    This step splits them into independent files and updates the cluster
    properties CSV so the C++ inversion treats each as a separate group.

    Calls the standalone ``split_mesh_volumes.py`` tool as a subprocess.
    """
    import subprocess

    print("\n" + "=" * 70)
    print("STEP 12 -- SPLIT DISCONNECTED VOLUMES")
    print("=" * 70)

    # Split in-place within combination_clusters
    meshes_dir = os.path.join(vol_clean_dir, "meshes")
    if not os.path.isdir(meshes_dir):
        print(f"  WARNING: meshes directory not found: {meshes_dir}")
        return vol_clean_dir

    split_script = str(PROJ_ROOT / "python" / "split_mesh_volumes.py")
    csv_path = os.path.join(vol_clean_dir, "cluster_properties.csv")

    print(f"  Running: python {split_script} {meshes_dir}/ --all --csv {csv_path} --formats ts,obj")

    result = subprocess.run(
        [sys.executable, split_script,
         meshes_dir + "/", "--all",
         "--csv", csv_path,
         "--formats", "ts,obj"],
        capture_output=True, text=True, cwd=str(PROJ_ROOT / "python"),
    )
    if result.stdout:
        for line in result.stdout.strip().split('\n'):
            print(f"  {line}")
    if result.returncode != 0:
        print(f"  WARNING: split_mesh_volumes returned exit code {result.returncode}")
        if result.stderr:
            print(f"  stderr: {result.stderr[:500]}")

    print(f"  Split volumes -> {vol_clean_dir}/")
    return vol_clean_dir


# ===================================================================
# STEP 13 -- Single-Property Clustering
# ===================================================================

def step13_single_property_clustering(mesh_data, rec_density, rec_susc,
                                        summary, output_dir,
                                        local_origin, z_datum):
    """
    Run single-property clustering for comparison.

    Produces:
      - density_clusters/   -- density-only GMM -> closed volumes
      - susc_clusters/      -- susceptibility-only GMM -> closed volumes

    These provide alternative starting models.  The intersection
    clustering (density x susceptibility) is the primary model.
    """
    from cluster_api._cluster import cluster_lithology, LithologySummary
    from lithoseed._pipeline import run_extract_group_volumes

    print("\n" + "=" * 70)
    print("STEP 13 -- SINGLE-PROPERTY CLUSTERING")
    print("=" * 70)

    zero = np.zeros_like(rec_density)

    # -- Density-only ----------------------------------------------------
    print("\n  [13a] Density-only clustering...")
    dens_labels, dens_summary = cluster_lithology(
        rec_density, zero, n_clusters=N_DENSITY_CLUSTERS, random_state=42)
    n_dens = len(dens_summary.labels)
    print(f"  Density-only: {n_dens} groups")
    for i in range(n_dens):
        print(f"    Group {dens_summary.labels[i]}: "
              f"mean={dens_summary.density_mean[i]:.4f} +/- "
              f"{dens_summary.density_std[i]:.4f}, "
              f"ncells={dens_summary.counts[i]}")

    dens_ls = LithologySummary(
        labels=dens_summary.labels,
        counts=dens_summary.counts,
        density_mean=dens_summary.density_mean,
        density_std=dens_summary.density_std,
        susc_mean=[0.0] * n_dens,
        susc_std=[0.0] * n_dens,
    )

    dens_dir = os.path.join(output_dir, "density_clusters")
    run_extract_group_volumes(
        mesh_data, dens_labels, dens_ls,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=dens_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=CLEAN_MIN_CELLS,
        igrf_f=IGRF_F, igrf_i=IGRF_I, igrf_d=IGRF_D,
    )
    print(f"  Density volumes -> {dens_dir}/")

    # -- Susceptibility-only ---------------------------------------------
    print("\n  [13b] Susceptibility-only clustering...")
    susc_labels, susc_summary = cluster_lithology(
        zero, rec_susc, n_clusters=N_SUSC_CLUSTERS, random_state=42)
    n_susc = len(susc_summary.labels)
    print(f"  Susceptibility-only: {n_susc} groups")
    for i in range(n_susc):
        print(f"    Group {susc_summary.labels[i]}: "
              f"mean={susc_summary.susc_mean[i]:.6f} +/- "
              f"{susc_summary.susc_std[i]:.6f}, "
              f"ncells={susc_summary.counts[i]}")

    susc_ls = LithologySummary(
        labels=susc_summary.labels,
        counts=susc_summary.counts,
        density_mean=[0.0] * n_susc,
        density_std=[0.0] * n_susc,
        susc_mean=susc_summary.susc_mean,
        susc_std=susc_summary.susc_std,
    )

    susc_dir = os.path.join(output_dir, "susc_clusters")
    run_extract_group_volumes(
        mesh_data, susc_labels, susc_ls,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=susc_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=CLEAN_MIN_CELLS,
        igrf_f=IGRF_F, igrf_i=IGRF_I, igrf_d=IGRF_D,
    )
    print(f"  Susceptibility volumes -> {susc_dir}/")


# ===================================================================
# STEP 13b -- Perturbed Lithology Gating (manual property ranges)
# ===================================================================

# Forrestania lithology density gates (absolute g/cm^3).
# ORDERED from densest to least-dense — first-match wins, so more
# specific mafic/ultramafic gates claim cells before Granite_Gneiss.
# Denser sections have finer divisions for perturbation analysis.
FORRESTANIA_GATES_V1 = [
    # name                 d_min  d_max  s_min   s_max    (SI)
    ("Massive_Sulfide",     2.90,  3.50,  0.030,  0.500),
    ("Ultramafic",          2.82,  3.00,  0.015,  0.200),
    ("Mafic_High",          2.76,  2.88,  0.008,  0.060),
    ("Mafic_Mid",           2.72,  2.80,  0.004,  0.030),
    ("Mafic_Low",           2.69,  2.74,  0.001,  0.015),
    ("Granite_Gneiss",      2.66,  2.71,  0.000,  0.008),
    ("Felsic",              2.60,  2.68,  0.000,  0.004),
    ("Background",          2.30,  2.63,  0.000,  0.003),
]

# V2: GMM density-only -> finer density divisions
# No manual susceptibility gates; pure density-driven
FORRESTANIA_GATES_V2 = None  # signal to use GMM instead of gating
USER_CLUSTER_N = 16  # number of GMM density-only clusters for user_clusters

FORRESTANIA_GATES = FORRESTANIA_GATES_V2  # active: use GMM mode


def step13b_perturbed_gating(mesh_data, rec_density, rec_susc,
                               output_dir, local_origin, z_datum):
    """
    Generate a second, perturbed clustering via manual property-range gating.

    Uses FORRESTANIA_GATES to define lithology-based density/susceptibility
    windows, splitting the denser mafic/ultramafic section into finer
    divisions than the GMM intersection clustering.

    Calls the same volume extraction pipeline as the GMM path, producing
    a parallel set of closed volumes at::

        volumes/user_clusters/
    """
    from cluster_api._cluster import LithologySummary, cluster_lithology
    from lithoseed._pipeline import run_extract_group_volumes

    print("\n" + "=" * 70)
    print("STEP 13b -- PERTURBED CLUSTERING")
    print("=" * 70)

    if FORRESTANIA_GATES is None:
        # ── V2: GMM density-only with configurable bin count ─────────────
        n_dens = USER_CLUSTER_N
        print(f"\n  GMM density-only clustering ({n_dens} groups)...")
        gated_labels, dens_summary = cluster_lithology(
            rec_density, np.zeros_like(rec_density),
            n_clusters=n_dens, random_state=42)

        n_groups = len(dens_summary.labels)
        print(f"  Density-only: {n_groups} groups")
        for i in range(n_groups):
            print(f"    Group {dens_summary.labels[i]}: "
                  f"mean={dens_summary.density_mean[i]:.4f} +/- "
                  f"{dens_summary.density_std[i]:.4f}, "
                  f"ncells={dens_summary.counts[i]}")

        summary = LithologySummary(
            labels=dens_summary.labels,
            counts=dens_summary.counts,
            density_mean=dens_summary.density_mean,
            density_std=dens_summary.density_std,
            susc_mean=[0.0] * n_groups,
            susc_std=[0.0] * n_groups,
        )
    else:
        # ── V1: Manual lithology gates ──────────────────────────────────
        gate_clusters_path = str(PROJ_ROOT / "python" / "gate_clusters.py")
        import importlib.util
        spec = importlib.util.spec_from_file_location("gate_clusters", gate_clusters_path)
        gc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gc)

        cluster_defs = []
        for name, d_min, d_max, s_min, s_max in FORRESTANIA_GATES:
            cluster_defs.append({
                "name": name, "d_min": d_min, "d_max": d_max,
                "s_min": s_min, "s_max": s_max,
                "has_susc": True, "source": "user",
            })

        print(f"\n  Density contrast range: [{rec_density.min():.4f}, {rec_density.max():.4f}]")
        print(f"  Defining {len(cluster_defs)} lithology gates (absolute density):")
        for cdef in cluster_defs:
            print(f"    {cdef['name']}: d=[{cdef['d_min']:.2f}, {cdef['d_max']:.2f}],"
                  f" s=[{cdef['s_min']:.3f}, {cdef['s_max']:.3f}]")

        bg = BACKGROUND_DENSITY
        for cdef in cluster_defs:
            cdef["d_min"] -= bg; cdef["d_max"] -= bg

        print(f"\n  Gating cells ({len(cluster_defs)} gates, contrast space):")
        gated_labels = gc.gate_cells(rec_density, rec_susc, cluster_defs)
        n_assigned = int((gated_labels > 0).sum())
        n_unassigned = int((gated_labels < 0).sum())
        print(f"  Assigned: {n_assigned}/{len(rec_density)}, unassigned: {n_unassigned}")

        if n_unassigned > 50:
            unassigned_mask = gated_labels < 0
            n_auto = min(3, n_unassigned // 50)
            print(f"  Auto-clustering {n_unassigned} remainder into {n_auto} groups...")
            auto_labels, _ = cluster_lithology(
                rec_density[unassigned_mask], rec_susc[unassigned_mask],
                n_clusters=n_auto, random_state=42)
            n_user = len(cluster_defs)
            gated_labels[unassigned_mask] = np.where(auto_labels > 0, auto_labels + n_user, -1)
            for i in range(n_auto):
                cluster_defs.append({
                    "name": f"Auto_{i+1}", "d_min": 0, "d_max": 0,
                    "s_min": 0, "s_max": 0, "has_susc": True, "source": "auto",
                })
            print(f"  Total: {len(cluster_defs)} ({n_user} lithology + {n_auto} auto)")

        summary = gc.build_summary(gated_labels, rec_density, rec_susc,
                                     cluster_defs, True)
    # Extract volumes
    pert_dir = os.path.join(output_dir, "user_clusters")
    print(f"\n  Extracting volumes -> {pert_dir}/")
    run_extract_group_volumes(
        mesh_data, gated_labels, summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=pert_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=CLEAN_MIN_CELLS,
        has_susceptibility=(FORRESTANIA_GATES is not None),
        min_density_range=0.2 if FORRESTANIA_GATES is None else 0.0,
    )

    print(f"  Perturbed volumes -> {pert_dir}/")
    return pert_dir


# ===================================================================
# STEP 13c -- Split All Volumes
# ===================================================================

def step13c_split_all_volumes(pres_dir):
    """
    Mirror volumes/ -> volumes_split/ with all meshes split into
    disconnected components via split_mesh_volumes.py.

    Processes all 4 cluster folders: combination_clusters, density_clusters,
    susc_clusters, user_clusters.
    """
    import subprocess

    print("\n" + "=" * 70)
    print("STEP 13c -- SPLIT ALL VOLUMES")
    print("=" * 70)

    src_volumes = os.path.join(pres_dir, "volumes")
    dst_volumes = os.path.join(pres_dir, "volumes_split")
    split_script = str(PROJ_ROOT / "python" / "split_mesh_volumes.py")

    # Copy entire volumes tree
    if os.path.exists(dst_volumes):
        shutil.rmtree(dst_volumes)
    shutil.copytree(src_volumes, dst_volumes)
    print(f"  Copied volumes/ -> volumes_split/")

    # Split each cluster folder's meshes
    for folder in ["combination_clusters", "density_clusters",
                    "susc_clusters", "user_clusters"]:
        meshes_dir = os.path.join(dst_volumes, folder, "meshes")
        if not os.path.isdir(meshes_dir):
            print(f"  SKIP {folder}/meshes/ (not found)")
            continue

        csv_path = os.path.join(dst_volumes, folder, "cluster_properties.csv")
        if not os.path.isfile(csv_path):
            csv_path = os.path.join(dst_volumes, folder, "cluster_properties.csv")

        print(f"  Splitting {folder}/ ...")
        result = subprocess.run(
            [sys.executable, split_script,
             meshes_dir + "/", "--all",
             "--csv", csv_path,
             "--formats", "ts,obj"],
            capture_output=True, text=True,
            cwd=str(PROJ_ROOT / "python"),
        )
        if result.returncode != 0 and result.returncode != 1:
            print(f"    WARNING: exit code {result.returncode}")
        n_ts = len([f for f in os.listdir(meshes_dir) if f.endswith('.ts')])
        print(f"    -> {n_ts} .ts files after split")

    print(f"\n  volumes_split/ ready: {dst_volumes}")


# ===================================================================
# STEP 14 -- Write INI Configs for C++ Inversion
# ===================================================================

def step14_write_ini_configs(local_origin, z_datum, output_dir):
    """
    Write INI configs for the C++ ``forrestania_invert.exe``.

    Two variants:
      - baseline/   -- standard cluster properties from the clustering step
      - perturbed/  -- user-perturbed physical properties for sensitivity testing

    Both point to the split volumes from step 12.
    """
    print("\n" + "=" * 70)
    print("STEP 14 -- INI CONFIGS FOR C++ INVERSION")
    print("=" * 70)

    baseline_dir = os.path.join(output_dir, "configs", "baseline")
    perturbed_dir = os.path.join(output_dir, "configs", "perturbed")
    os.makedirs(baseline_dir, exist_ok=True)
    os.makedirs(perturbed_dir, exist_ok=True)

    # -- Baseline config -------------------------------------------------
    # Paths are relative to configs/baseline/resolved_config.ini
    volumes_rel = "../../volumes/combination_clusters"
    baseline_ini = _build_ini_config(
        local_origin, z_datum,
        volumes_rel=volumes_rel,
        output_subdir="inprogress_baseline",
    )
    ini_path = os.path.join(baseline_dir, "resolved_config.ini")
    with open(ini_path, 'w') as f:
        f.write(baseline_ini)
    print(f"  Wrote {ini_path}")

    # -- Perturbed config ------------------------------------------------
    perturbed_ini = _build_ini_config(
        local_origin, z_datum,
        volumes_rel=volumes_rel,
        output_subdir="inprogress_perturbed",
    )
    ini_path = os.path.join(perturbed_dir, "resolved_config.ini")
    with open(ini_path, 'w') as f:
        f.write(perturbed_ini)
    print(f"  Wrote {ini_path}")

    print(f"\n  To run C++ inversion:")
    print(f"    cd build\\release")
    print(f"    .\\forrestania_invert.exe "
          f"..\\..\\apps\\forrestania\\presentation\\configs\\baseline\\resolved_config.ini")


def _build_ini_config(local_origin, z_datum, volumes_rel, output_subdir):
    """Build INI config string for C++ inversion.

    Parameters
    ----------
    volumes_rel : str
        Relative path from the config directory to the split volumes directory.
        E.g. ``../../volumes/combination_clusters``
    """
    ox, oy = local_origin

    # Discover volume group mesh files from the volumes directory
    meshes_abs = os.path.join(str(PRESENTATION), "configs", "baseline", volumes_rel, "meshes")
    group_meshes = []
    if os.path.isdir(meshes_abs):
        for fname in sorted(os.listdir(meshes_abs)):
            if fname.startswith('volume_group_') and fname.endswith('.ts'):
                # Paths in group_meshes are relative to the config file location
                group_meshes.append(f"{volumes_rel}/meshes/{fname}")
    group_meshes_str = ",".join(group_meshes) if group_meshes else f"{volumes_rel}/meshes/volume_group_1.ts"

    # Build relative paths for CSV files
    cluster_csv = f"{volumes_rel}/cluster_properties.csv"
    obs_grav_csv = f"{volumes_rel}/observed_gravity.csv"
    obs_mag_csv = f"{volumes_rel}/observed_magnetic.csv"

    # Fall back to presentation root if files don't exist in volumes dir
    import glob as glob_mod
    if not os.path.isfile(os.path.join(str(PRESENTATION), "configs", "baseline", cluster_csv)):
        cluster_csv = "../../cluster_properties.csv"
        obs_grav_csv = "../../observed_gravity.csv"
        obs_mag_csv = "../../observed_magnetic.csv"

    ini = f"""[bounds]
depth_bound_margin = 500.000000
enable_depth_bounds = true

[data]
real_data = true
group_column = cluster_id
cluster_csv = {cluster_csv}
group_meshes = {group_meshes_str}
label_grid = {volumes_rel}/meshes/label_grid.bin
local_origin_x = {ox:.6f}
local_origin_y = {oy:.6f}
observed_gravity_csv = {obs_grav_csv}
observed_magnetic_csv = {obs_mag_csv}
z_datum = {z_datum:.6f}
fail_on_missing_meshes = false

[gravity]
density_max = 0.500000
density_min = -0.500000
gravity_uncertainty = 0.000000

[inversion]
solver = lbfgsb
max_iterations = 250
tolerance = 1e-05
lambda = 1.0
omega = 10.0
lbfgs_history = 20
fd_step = 10.0
control_point_stride = 0
armijo_c1 = 0.001
line_search_max_iter = 30
vertex_freedom = xyz_free
mesh_boundary_mode = free
enable_eigenvalue_scaling = true
enable_geometry_inversion = true
enable_property_inversion = false
property_inversion_interval = 1
property_inversion_max_iter = 10
property_damping = 0.010000
gncg_cg_max_iter = 50
gncg_cg_tolerance = 0.000001
disable_line_search = false

[magnetic]
declination = {IGRF_D:.6f}
field_nT = {IGRF_F:.6f}
inclination = {IGRF_I:.6f}
magnetic_uncertainty = 0.000000
magnetic_weight = 1.000000
remanence_component_max = 10.000000
remanence_component_min = -10.000000
remanence_max = 10.000000
remanence_min = 0.000000
remanence_mode = effective_susceptibility
susceptibility_max = 0.100000
susceptibility_min = 0.000000

[output]
iteration_export_dir = ./{output_subdir}

[padding]
enable_padding_group = false
padding_conductivity_initial = 0.000100
padding_conductivity_lower = 0.000010
padding_conductivity_upper = 100.000000
padding_density_initial = 2.680000
padding_density_lower = 1.500000
padding_density_upper = 5.000000
padding_depth = -100000.000000
padding_susceptibility_initial = 0.000000
padding_susceptibility_lower = 0.000000
padding_susceptibility_upper = 1.000000

[regularization]
enable_reference_model = false
gradient_smoothing_weight = 0.300000
lambda_ref = 0.100000

[topography]
mode = none
dem_file =
datum_elevation = {z_datum:.6f}
bouguer_density = 2.67
padding_rings = 0
padding_cell_size = 100.0
invert_halfspace_properties = false
projected_beyond_survey_area = independent_flat

[true_model]
dip_angle = 0.000000
dip_direction = 90.000000
mode = voiseys
"""
    return ini


# ===================================================================
# STEP 15 -- Presentation Directory Assembly
# ===================================================================

def step15_assemble_presentation(output_dir):
    """
    Collect all outputs into the final presentation directory structure.

    Copies key outputs from the working directories to PRESENTATION/ for
    easy access when preparing slides and figures.
    """
    print("\n" + "=" * 70)
    print("STEP 15 -- ASSEMBLE PRESENTATION DIRECTORY")
    print("=" * 70)

    # Everything was already written directly to PRESENTATION/ by previous steps.
    # This step prints a directory tree for verification.

    def _tree(path, prefix="", max_depth=3, depth=0):
        if depth > max_depth:
            return
        items = sorted(os.listdir(path)) if os.path.isdir(path) else []
        for i, item in enumerate(items):
            is_last = i == len(items) - 1
            item_path = os.path.join(path, item)
            connector = "\-- " if is_last else "|-- "
            print(f"{prefix}{connector}{item}")
            if os.path.isdir(item_path) and depth < max_depth:
                ext = "    " if is_last else "|   "
                _tree(item_path, prefix + ext, max_depth, depth + 1)

    print(f"\n  {output_dir}/")
    _tree(output_dir, "  ", max_depth=3)

    print(f"\n  [OK] Presentation materials assembled in {output_dir}/")


# ===================================================================
# STEP 16 -- Data Lineage (raw data + processed intermediates)
# ===================================================================

def step16_data_lineage(output_dir, grav_xyz, g_obs, mag_xyz, t_obs,
                         topo, rec_density, rec_density_abs, rec_susc,
                         mesh_params, labels, summary,
                         local_origin, z_datum):
    """
    Copy raw data and every processed intermediate into a lineage/ folder.

    Each sub-step gets its own directory with:
      - The data file(s) at that stage
      - A README.txt describing exactly what transformation was applied

    This makes the data provenance auditable for any presentation slide.
    """
    print("\n" + "=" * 70)
    print("STEP 16 -- DATA LINEAGE (raw + processed intermediates)")
    print("=" * 70)

    lineage = os.path.join(output_dir, "lineage")
    os.makedirs(lineage, exist_ok=True)

    ox, oy = local_origin
    nx, ny, nz, dx, dy, dz, x0, y0, z_bot = mesh_params

    # ── 01: Raw gravity ─────────────────────────────────────────────────
    d = os.path.join(lineage, "01_raw_gravity")
    os.makedirs(d, exist_ok=True)
    src = DATASETS / "Forrestania_Gravity_Station_trim_.csv"
    shutil.copy2(str(src), os.path.join(d, src.name))
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""RAW GRAVITY DATA
================
Source: {src.name}
Stations: {len(grav_xyz)}
Columns: X (easting_mga50), Y (northing_mga50), Z (elevation_m), FGrav_mgal
Coordinate system: UTM MGA Zone 50
No transformations applied.
""")
    print(f"  {d}/")

    # ── 02: Raw magnetics ───────────────────────────────────────────────
    d = os.path.join(lineage, "02_raw_magnetics")
    os.makedirs(d, exist_ok=True)
    src = DATASETS / "60472_AOI4_regional_mag_Zfixed.csv"
    shutil.copy2(str(src), os.path.join(d, src.name))
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""RAW MAGNETIC DATA
=================
Source: {src.name}
Stations: {len(mag_xyz)} (after downsampling from source file)
Columns: X, Y, Z, MAGCOMP (total field), IGRF
Coordinate system: UTM MGA Zone 50
No transformations applied at this stage.
""")
    print(f"  {d}/")

    # ── 03: Processed gravity ───────────────────────────────────────────
    d = os.path.join(lineage, "03_processed_gravity")
    os.makedirs(d, exist_ok=True)
    # Write local-coordinate gravity
    grav_local = os.path.join(d, "observed_gravity.csv")
    with open(grav_local, 'w') as f:
        f.write("x,y,z,g_obs\n")
        for i in range(len(g_obs)):
            f.write(f"{grav_xyz[i,0] - ox:.6f},{grav_xyz[i,1] - oy:.6f},"
                    f"{grav_xyz[i,2] - z_datum:.6f},{g_obs[i]:.6f}\n")
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""PROCESSED GRAVITY DATA
=====================
Stations: {len(grav_xyz)}

Transformations applied:
  1. SIGN FLIP: g_obs = -FGrav_mgal
     Reason: SimPEG convention has positive gz downward.
     Original range: see raw file.
     Flipped range: [{g_obs.min():.4f}, {g_obs.max():.4f}] mGal

  2. LOCAL COORDINATE TRANSLATION:
     x_local = X_utm - {ox:.1f}
     y_local = Y_utm - {oy:.1f}
     z_local = Z_utm - {z_datum:.1f}  (datum = mean survey elevation)

Output columns: x, y, z (local coords, meters), g_obs (mGal, sign-flipped)
""")
    print(f"  {d}/")

    # ── 04: Processed magnetics ─────────────────────────────────────────
    d = os.path.join(lineage, "04_processed_magnetics")
    os.makedirs(d, exist_ok=True)
    mag_local = os.path.join(d, "observed_magnetic.csv")
    with open(mag_local, 'w') as f:
        f.write("x,y,z,t_obs\n")
        for i in range(len(t_obs)):
            f.write(f"{mag_xyz[i,0] - ox:.6f},{mag_xyz[i,1] - oy:.6f},"
                    f"{mag_xyz[i,2] - z_datum:.6f},{t_obs[i]:.6f}\n")
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""PROCESSED MAGNETIC DATA
=======================
Stations: {len(t_obs)} (downsampled from source)

Transformations applied:
  1. RMI COMPUTATION: RMI = MAGCOMP - IGRF
     IGRF at Forrestania (epoch 2026.0):
       F = {IGRF_F:.0f} nT
       I = {IGRF_I:.1f} deg
       D = {IGRF_D:.1f} deg

  2. RANDOM DOWNSAMPLING: {len(t_obs)} stations selected
     (seed=42, from full survey)

  3. PLANAR DETRENDING:
     Trend = a*x + b*y + c  fitted via least-squares
     t_obs = RMI - trend
     This removes regional field gradient.
     t_obs range: [{t_obs.min():.1f}, {t_obs.max():.1f}] nT

  4. LOCAL COORDINATE TRANSLATION (same origin as gravity):
     x_local = X_utm - {ox:.1f}
     y_local = Y_utm - {oy:.1f}
     z_local = Z_utm - {z_datum:.1f}

Output columns: x, y, z (local coords, meters), t_obs (nT, detrended RMI)
""")
    print(f"  {d}/")

    # ── 05: DEM ─────────────────────────────────────────────────────────
    d = os.path.join(lineage, "05_dem")
    os.makedirs(d, exist_ok=True)
    dem_csv = os.path.join(d, "dem_points.csv")
    with open(dem_csv, 'w') as f:
        f.write("x_utm,y_utm,z_elevation\n")
        for pt in topo:
            f.write(f"{pt[0]:.6f},{pt[1]:.6f},{pt[2]:.6f}\n")
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""DIGITAL ELEVATION MODEL (from station elevations)
=============================================
Points: {len(topo)}
Spacing: ~120 m (interpolated from gravity station elevations)
Method: linear interpolation via scipy.interpolate.griddata

Z range: [{topo[:,2].min():.1f}, {topo[:,2].max():.1f}] m

Used for:
  - Building the tensor mesh top surface
  - active_from_xyz() to mark cells below topography as active
""")
    print(f"  {d}/")

    # ── 06: Tensor mesh ─────────────────────────────────────────────────
    d = os.path.join(lineage, "06_tensor_mesh")
    os.makedirs(d, exist_ok=True)
    write_mesh_ubc(os.path.join(d, "simpeg_mesh.msh"),
                   x0, y0, z_bot + nz * dz, dx, dy, dz, nx, ny, nz)
    with open(os.path.join(d, "README.txt"), 'w') as f:
        z_top = z_bot + nz * dz
        f.write(f"""TENSOR MESH
============
Dimensions: {nx} x {ny} x {nz} = {nx*ny*nz} cells
Cell size: {dx} x {dy} x {dz} m
Origin (UTM, SW top corner): ({x0:.1f}, {y0:.1f}, {z_top:.1f})
Z range: [{z_bot:.1f}, {z_top:.1f}] m

BUILD METHOD (FIXED 2026-06-05):
  x/y: 2-cell padding beyond data extent
  z_top: ceil(max_elevation / cell_size) * cell_size
         (NO overhead padding -- snaps to ground level)
  z_bot: floor(min_elevation - depth_core) / cell_size * cell_size
         with 2-cell padding below
  depth_core = {DEPTH_CORE} m

FIX VERIFICATION:
  Top cell center at {z_top - dz/2:.0f} m
  Max topography = {topo[:,2].max():.1f} m
  -> {z_top - dz/2:.0f} <= {topo[:,2].max():.1f} -> top cell entirely below ground
  Air cells: 0 / {nx*ny*nz}

Format: UBC-GIF tensor mesh
  Line 1: nx ny nz
  Line 2: x0 y0 z0 (SW top corner)
  Lines 3-5: cell sizes in x, y, z
""")
    print(f"  {d}/")

    # ── 07: SimPEG density inversion ────────────────────────────────────
    d = os.path.join(lineage, "07_simpeg_density")
    os.makedirs(d, exist_ok=True)
    # Density contrast (UBC-GIF order: z fastest, then x, then y)
    write_mod_ga(os.path.join(d, "simpeg_density_contrast.mod"),
                 rec_density, nx, ny, nz, order='F')
    # Absolute density
    write_mod_ga(os.path.join(d, "simpeg_density_absolute.mod"),
                 rec_density_abs, nx, ny, nz, order='F')
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""SimPEG GRAVITY INVERSION -- RECOVERED DENSITY
=============================================
Method: L2 data misfit, WeightedLeastSquares regularization (alpha_z=2.0)
Bounds: [-0.4, +0.6] g/cm^3 density contrast
Optimizer: ProjectedGNCG, max {MAX_ITER} iterations
Mesh: {nx}x{ny}x{nz} = {nx*ny*nz} cells, {dx}m uniform

Results:
  Density contrast range: [{rec_density.min():.4f}, {rec_density.max():.4f}] g/cm^3
  Absolute density range: [{rec_density_abs.min():.4f}, {rec_density_abs.max():.4f}] g/cm^3
  (absolute = contrast + {BACKGROUND_DENSITY} g/cm^3 background)

Files:
  simpeg_density_contrast.mod  -- recovered density contrast (no header)
  simpeg_density_absolute.mod  -- absolute density = contrast + {BACKGROUND_DENSITY}

GA COMPATIBILITY:
  No header line (GA reads geometry from .msh file).
  Values in UBC-GIF order: z fastest (top-down), then x (easting),
  then y (northing).  Air cells set to {SENTINEL} sentinel.
""")
    print(f"  {d}/")

    # ── 08: SimPEG susceptibility inversion ─────────────────────────────
    d = os.path.join(lineage, "08_simpeg_susceptibility")
    os.makedirs(d, exist_ok=True)
    write_mod_ga(os.path.join(d, "simpeg_susceptibility.mod"),
                 rec_susc, nx, ny, nz, order='F')
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""SimPEG MAGNETIC INVERSION -- RECOVERED SUSCEPTIBILITY
===================================================
Method: L2 data misfit, WeightedLeastSquares regularization (alpha_z=2.0)
Bounds: [0.0, +0.5] SI susceptibility
Optimizer: ProjectedGNCG, max {MAX_ITER} iterations
Mesh: {nx}x{ny}x{nz} = {nx*ny*nz} cells, {dx}m uniform
Inducing field: F={IGRF_F:.0f} nT, I={IGRF_I:.1f} deg, D={IGRF_D:.1f} deg

Results:
  Susceptibility range: [{rec_susc.min():.6f}, {rec_susc.max():.6f}] SI

Files:
  simpeg_susceptibility.mod  -- recovered susceptibility (no header)

GA COMPATIBILITY:
  No header line. UBC-GIF order (z fastest, then x, then y). Air cells = {SENTINEL}.
""")
    print(f"  {d}/")

    # ── 09: Clustering ──────────────────────────────────────────────────
    d = os.path.join(lineage, "09_clustering")
    os.makedirs(d, exist_ok=True)
    # Write cluster properties
    csv_path = os.path.join(d, "cluster_properties.csv")
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["cluster_id", "working_name", "sample_count",
                     "density_median_gcc", "density_p10", "density_p90",
                     "susceptibility_median_SI", "susceptibility_p10",
                     "susceptibility_p90"])
        if hasattr(summary, 'labels'):
            for i, gid in enumerate(summary.labels):
                d_mean = summary.density_mean[i] if hasattr(summary, 'density_mean') else 0
                d_std = summary.density_std[i] if hasattr(summary, 'density_std') else 0
                s_mean = summary.susc_mean[i] if hasattr(summary, 'susc_mean') else 0
                s_std = summary.susc_std[i] if hasattr(summary, 'susc_std') else 0
                count = summary.counts[i] if hasattr(summary, 'counts') else 0
                w.writerow([gid, f"group_{gid}", count,
                            f"{d_mean:.6f}", f"{d_mean - d_std:.6f}", f"{d_mean + d_std:.6f}",
                            f"{s_mean:.6e}", f"{max(0, s_mean - s_std):.6e}", f"{s_mean + s_std:.6e}"])

    n_groups = summary.n_groups if hasattr(summary, 'n_groups') else '?'
    with open(os.path.join(d, "README.txt"), 'w') as f:
        f.write(f"""GMM INTERSECTION CLUSTERING
===========================
Method: 1D GMM on density x 1D GMM on susceptibility
        -> intersection of {N_DENSITY_CLUSTERS} density bins x {N_SUSC_CLUSTERS} susc bins
        -> {n_groups} non-empty intersection groups

Input:
  - Recovered density contrast (from SimPEG gravity inversion)
  - Recovered susceptibility (from SimPEG magnetic inversion)
  - Both on the same {nx}x{ny}x{nz} tensor mesh ({dx}m cells)

Algorithm:
  cluster_intersection() from cluster_api._cluster
  sklearn.mixture.GaussianMixture (1D, full covariance)
  random_state=42

Output groups ({n_groups} total):
""")
        if hasattr(summary, 'labels'):
            for i, gid in enumerate(summary.labels):
                d_mean = summary.density_mean[i]
                s_mean = summary.susc_mean[i]
                count = summary.counts[i]
                f.write(f"  Group {gid}: density={d_mean:.4f} g/cm^3, "
                        f"susc={s_mean:.6e} SI, ncells={count}\n")

    print(f"  {d}/")

    # ── 10: Volume meshes ───────────────────────────────────────────────
    d = os.path.join(lineage, "10_volume_meshes")
    os.makedirs(d, exist_ok=True)
    # Copy from split volumes
    src_meshes = os.path.join(output_dir, "volumes", "combination_clusters", "meshes")
    if os.path.isdir(src_meshes):
        for fname in os.listdir(src_meshes):
            if fname.endswith('.ts') and fname.startswith('volume_group_'):
                shutil.copy2(os.path.join(src_meshes, fname),
                            os.path.join(d, fname))
    with open(os.path.join(d, "README.txt"), 'w') as f:
        n_vols = len([x for x in os.listdir(d) if x.endswith('.ts')])
        f.write(f"""VOLUME MESHES (closed, zero-gap)
=================================
Meshes: {n_vols} volume .ts files

Method: Voxel-face boundary extraction
  - Every face in the 3D label grid is iterated once
  - Shared quads assigned to both adjacent groups with opposite winding
  - Guarantees zero-gap, shared-vertex, perfectly squared-off closed surfaces
  - No decimation applied -- shared boundary vertices remain identical

Pipeline:
  1. cluster_intersection() -> labels on tensor mesh
  2. cleanup_label_fragments() -> reassign fragments < {CLEAN_MIN_CELLS} cells
  3. extract_group_volumes() -> voxel-face TS + OBJ export
  4. split_mesh_volumes.py -> detect + split disconnected components

Format: GOCAD TSurf 1 (TFACE keyword, 8-digit scientific precision)
Coordinate system: Model space (local origin at {ox:.1f}, {oy:.1f})

Files: volume_group_N.ts = closed volume for lithology group N
""")
    print(f"  {d}/")

    print(f"\n  Data lineage: {lineage}/")
    print(f"  10 subdirectories tracing every data transformation")


# ===================================================================
# STEP 17 -- Convenience Copy (flat local folder)
# ===================================================================

def step17_convenience_copy(output_dir):
    """
    Copy all key outputs into a flat ``convenience/`` folder for easy access.

    Collects from across the presentation tree into one place:
      - All .mod and .msh files (GA-loadable SimPEG outputs)
      - All figures (.png)
      - All INI configs
      - Cluster properties CSV
      - All volume .ts meshes
    """
    print("\n" + "=" * 70)
    print("STEP 17 -- CONVENIENCE COPY (flat local folder)")
    print("=" * 70)

    conv_dir = os.path.join(output_dir, "convenience")
    os.makedirs(conv_dir, exist_ok=True)

    copied = 0

    # -- SimPEG meshes ----------------------------------------------------
    simpeg_dir = os.path.join(output_dir, "simpeg")
    if os.path.isdir(simpeg_dir):
        for fname in os.listdir(simpeg_dir):
            if fname.endswith(('.mod', '.msh')):
                shutil.copy2(os.path.join(simpeg_dir, fname),
                            os.path.join(conv_dir, fname))
                copied += 1
        # Also copy data_coords versions
        dc_dir = os.path.join(simpeg_dir, "data_coords")
        if os.path.isdir(dc_dir):
            for fname in os.listdir(dc_dir):
                if fname.endswith(('.mod', '.msh')):
                    shutil.copy2(os.path.join(dc_dir, fname),
                                os.path.join(conv_dir, "utm_" + fname))
                    copied += 1

    # -- Figures ----------------------------------------------------------
    figures_dir = os.path.join(output_dir, "figures")
    if os.path.isdir(figures_dir):
        for fname in os.listdir(figures_dir):
            if fname.endswith('.png'):
                shutil.copy2(os.path.join(figures_dir, fname),
                            os.path.join(conv_dir, fname))
                copied += 1

    # -- INI configs ------------------------------------------------------
    configs_dir = os.path.join(output_dir, "configs")
    if os.path.isdir(configs_dir):
        for variant in ["baseline", "perturbed"]:
            ini_path = os.path.join(configs_dir, variant, "resolved_config.ini")
            if os.path.isfile(ini_path):
                shutil.copy2(ini_path,
                            os.path.join(conv_dir, f"config_{variant}.ini"))
                copied += 1

    # -- Cluster CSV ------------------------------------------------------
    for loc in ["volumes/combination_clusters/cluster_properties.csv",
                "cluster_properties.csv"]:
        csv_path = os.path.join(output_dir, loc)
        if os.path.isfile(csv_path):
            shutil.copy2(csv_path, os.path.join(conv_dir, "cluster_properties.csv"))
            copied += 1
            break

    # -- Volume meshes (.ts files from split volumes) ---------------------
    vol_meshes = os.path.join(output_dir, "volumes", "combination_clusters", "meshes")
    if os.path.isdir(vol_meshes):
        conv_meshes = os.path.join(conv_dir, "volume_meshes")
        os.makedirs(conv_meshes, exist_ok=True)
        for fname in os.listdir(vol_meshes):
            if fname.endswith('.ts'):
                shutil.copy2(os.path.join(vol_meshes, fname),
                            os.path.join(conv_meshes, fname))
                copied += 1

    print(f"\n  Copied {copied} files -> {conv_dir}/")
    print(f"  Contents:")
    for fname in sorted(os.listdir(conv_dir)):
        fpath = os.path.join(conv_dir, fname)
        if os.path.isfile(fpath):
            size_kb = os.path.getsize(fpath) / 1024
            print(f"    {fname} ({size_kb:.0f} KB)")
        elif os.path.isdir(fpath):
            n = len(os.listdir(fpath))
            print(f"    {fname}/ ({n} files)")

    return conv_dir


# ===================================================================
# MAIN -- Orchestrate All Steps
# ===================================================================

def main():
    t_total = time.time()

    print("=" * 70)
    print("FORRESTANIA MASTER PIPELINE -- Presentation-Ready Outputs")
    print("=" * 70)
    print(f"Output root: {PRESENTATION}")
    print(f"Mesh cell size: {CELL_SIZE}m  |  Depth core: {DEPTH_CORE}m")
    print(f"Clustering: {N_DENSITY_CLUSTERS} density x {N_SUSC_CLUSTERS} susceptibility")
    print(f"GA sentinel: {SENTINEL}  |  Background density: {BACKGROUND_DENSITY} g/cm³")
    print()

    # -- Step 1: Data preparation ----------------------------------------
    grav_xyz, g_obs, mag_xyz, t_obs, topo, local_origin, z_datum = \
        step1_prepare_data()

    # -- Step 2: Build mesh (FIXED) --------------------------------------
    mesh, active, n_act, nx, ny, nz, dx, dy, dz, x0, y0, z_bot = \
        step2_build_mesh(topo, CELL_SIZE, DEPTH_CORE)
    mesh_params = (nx, ny, nz, dx, dy, dz, x0, y0, z_bot)  # for lineage

    # -- Step 3-4: Inversions --------------------------------------------
    simpeg_dir = os.path.join(str(PRESENTATION), "simpeg")
    os.makedirs(simpeg_dir, exist_ok=True)

    rec_density = step3_gravity_inversion(
        mesh, active, n_act, nx, ny, nz, grav_xyz, g_obs, simpeg_dir)
    rec_susc = step4_magnetic_inversion(
        mesh, active, n_act, nx, ny, nz, mag_xyz, t_obs, simpeg_dir)

    # -- Step 5: Absolute density ----------------------------------------
    rec_density_abs = step5_absolute_density(rec_density, BACKGROUND_DENSITY)

    # -- Step 6: GA-compatible export (model space) ----------------------
    air_mask = step6_export_ga_model_space(
        mesh, active, rec_density, rec_density_abs, rec_susc,
        nx, ny, nz, x0, y0, z_bot, dx, dy, dz, simpeg_dir)

    # -- Step 7: Data-coordinate copies ----------------------------------
    step7_export_data_coords(
        mesh, rec_density_abs, rec_susc, active,
        nx, ny, nz, x0, y0, z_bot, dx, dy, dz, simpeg_dir)

    # -- Step 8: Convert to MeshData -------------------------------------
    mesh_data = step8_convert_to_meshdata(mesh, active, rec_density, n_act)

    # -- Step 9: Intersection clustering ---------------------------------
    volumes_dir = os.path.join(str(PRESENTATION), "volumes")
    labels, summary, contact_result = step9_clustering(
        mesh_data, rec_density, rec_susc, volumes_dir,
        local_origin, z_datum, grav_xyz, g_obs, mag_xyz, t_obs)

    # -- Step 10: Volume extraction --------------------------------------
    clean_result = step10_volume_extraction(
        mesh_data, labels, summary, volumes_dir, local_origin, z_datum)

    # -- Step 11: All QC plots --------------------------------------------
    figures_dir = os.path.join(str(PRESENTATION), "figures")
    step11_all_plots(rec_density, rec_susc, labels, summary,
                      mesh, active, nx, ny, nz,
                      grav_xyz, g_obs, mag_xyz, t_obs,
                      simpeg_dir, figures_dir)

    # -- Step 12: Split disconnected volumes -----------------------------
    vol_clean_dir = os.path.join(volumes_dir, "combination_clusters")
    step12_split_volumes(vol_clean_dir, volumes_dir)

    # -- Step 13: Single-property clustering -----------------------------
    step13_single_property_clustering(
        mesh_data, rec_density, rec_susc, summary, volumes_dir,
        local_origin, z_datum)

    # -- Step 13b: Perturbed lithology gating ----------------------------
    step13b_perturbed_gating(
        mesh_data, rec_density, rec_susc,
        volumes_dir, local_origin, z_datum)

    # -- Step 13c: Split all volumes for presentation ---------------------
    step13c_split_all_volumes(str(PRESENTATION))

    # -- Step 14: INI configs --------------------------------------------
    step14_write_ini_configs(local_origin, z_datum, str(PRESENTATION))

    # -- Step 15: Assemble -----------------------------------------------
    step15_assemble_presentation(str(PRESENTATION))

    # -- Step 16: Data lineage -------------------------------------------
    step16_data_lineage(str(PRESENTATION), grav_xyz, g_obs, mag_xyz, t_obs,
                         topo, rec_density, rec_density_abs, rec_susc,
                         mesh_params, labels, summary, local_origin, z_datum)

    # -- Step 17: Convenience copy ---------------------------------------
    conv_dir = step17_convenience_copy(str(PRESENTATION))

    # -- Done ------------------------------------------------------------
    t_total = time.time() - t_total
    print("\n" + "=" * 70)
    print(f"PIPELINE COMPLETE -- {t_total/60:.1f} minutes total")
    print(f"All outputs: {PRESENTATION}/")
    print("=" * 70)
    print(f"\nNext steps:")
    print(f"  1. Grab all files: {conv_dir}\\")
    print(f"  2. Verify GA loads: {simpeg_dir}/simpeg_mesh.msh")
    print(f"  3. Run C++ inversion:")
    print(f"     cd build\\release")
    print(f"     .\\forrestania_invert.exe "
          f"{PRESENTATION}\\configs\\baseline\\resolved_config.ini")
    print(f"  4. For perturbed properties, edit cluster CSV then run:")
    print(f"     .\\forrestania_invert.exe "
          f"{PRESENTATION}\\configs\\perturbed\\resolved_config.ini")

    return 0


if __name__ == "__main__":
    sys.exit(main())
