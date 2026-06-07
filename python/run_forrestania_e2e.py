#!/usr/bin/env python3
"""
Forrestania end-to-end pipeline (Steps 1-2).

Generates starting-model contact surfaces + INI config for the C++ litho
inversion (Step 3).

Steps:
  1-3. Prep gravity, magnetic data + DEM
  4. Build shared tensor mesh → writes simpeg_mesh.msh
  5. Run SimPEG gravity inversion → writes simpeg_density.mod + convergence CSV
  6. Run SimPEG magnetic inversion → writes simpeg_susceptibility.mod + convergence CSV
  7. Intersection clustering (1D GMM density × 1D GMM susceptibility)
     → contact surface extraction + group closed volumes via LithoSeed
  8. Single-property clustering (density-only + susceptibility-only)
     → group closed volumes + scatter plots (histogram + cell-index scatter)
  9. Write local-coordinate CSVs + intersection scatter plot + validation
"""
import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

# Add local packages to path
PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))

DATASETS = os.path.join(PROJ_ROOT, "..", "datasets", "Forrestania")
OUTPUT_DIR = os.path.join(PROJ_ROOT, "..", "apps", "forrestania", "inversion_output")
FORRESTANIA_APP = os.path.join(PROJ_ROOT, "..", "apps", "forrestania")

IGRF_F = 58874.0
IGRF_I = -66.2
IGRF_D = -0.1
GRAV_UNCERTAINTY = 0.05
MAG_UNCERTAINTY = 50.0
CELL_SIZE = 200.0
DEPTH_CORE = 1000.0
MAX_ITER = 60
N_CLUSTERS = 4

RNG = np.random.default_rng(42)


# ── Shared utilities ────────────────────────────────────────────────────────

def _write_conv_csv(history: list, path: str) -> None:
    """Write a SimPEG convergence CSV from `ConvergenceRecorder.history`."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write("iter,beta,phi_d,phi_m,f\n")
        for h in history:
            f.write(f"{h['iter']},{h['beta']:.6e},{h['phi_d']:.6e},"
                    f"{h['phi_m']:.6e},{h['f']:.6e}\n")


def _make_convergence_recorder():
    """Return a ``ConvergenceRecorder`` class inheriting from SimPEG's
    ``InversionDirective``.  Must be called with SimPEG already imported.
    """
    from simpeg.directives import InversionDirective

    class _Recorder(InversionDirective):
        def __init__(self) -> None:
            super().__init__()
            self.history: list[dict[str, float]] = []

        def initialize(self) -> None:
            pass

        def endIter(self) -> None:
            self.history.append({
                "iter": int(self.opt.iter),
                "beta": float(self.invProb.beta),
                "phi_d": float(self.invProb.phi_d),
                "phi_m": float(self.invProb.phi_m),
                "f": float(self.invProb.phi_d + self.invProb.beta * self.invProb.phi_m),
            })

        def write_csv(self, path: str) -> None:
            _write_conv_csv(self.history, path)

    return _Recorder


def _save_model_ubc(path: str, model: np.ndarray, active: np.ndarray,
                    nx: int, ny: int, nz: int) -> None:
    """Write a UBC-GIF .mod file (model on full mesh, inactive cells = 0)."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    full = np.zeros(len(active))
    full[active] = model
    with open(path, "w") as f:
        f.write(f"{nx} {ny} {nz}\n")
        for val in full:
            f.write(f"{val:.6e}\n")


# ── Data preparation ────────────────────────────────────────────────────────


def prep_gravity(grav_csv):
    """Load and prep gravity data: sign flip for SimPEG convention."""
    df = pd.read_csv(grav_csv)
    xyz = df[["X", "Y", "Z"]].values
    g_obs = -df["FGrav_mgal"].values
    print(f"  Gravity: {len(df)} stations, g_obs range [{g_obs.min():.4f}, {g_obs.max():.4f}] mGal")
    return xyz, g_obs


def prep_magnetics(mag_csv, n_target=500):
    """Load and prep magnetic data: RMI, downsample, planar detrend."""
    df = pd.read_csv(mag_csv)
    df["RMI"] = df["MAGCOMP"] - df["IGRF"]

    n_target = min(n_target, len(df))
    idx = RNG.choice(len(df), n_target, replace=False)
    sub = df.iloc[idx].copy()

    xyz = sub[["X", "Y", "Z"]].values
    rmi = sub["RMI"].values
    A = np.c_[xyz[:, 0], xyz[:, 1], np.ones(len(xyz))]
    coeff, _, _, _ = np.linalg.lstsq(A, rmi, rcond=None)
    trend = A @ coeff
    t_obs = rmi - trend

    print(f"  Magnetics: {len(df)} raw -> {len(sub)} downsampled")
    print(f"  Trend coefficients: {coeff[0]:.6e}, {coeff[1]:.6e}, intercept={coeff[2]:.2f} nT")
    print(f"  RMI_detrended range [{t_obs.min():.1f}, {t_obs.max():.1f}] nT")
    return xyz, t_obs


def build_dem(grav_xyz, mag_xyz, mag_ground_z=None, spacing=120.0, pad=500.0):
    """Build DEM grid from ground elevations."""
    if mag_ground_z is not None:
        all_pts = np.vstack([
            grav_xyz,
            np.c_[mag_xyz[:, 0], mag_xyz[:, 1], mag_ground_z],
        ])
    else:
        all_pts = grav_xyz
    x_min, x_max = all_pts[:, 0].min() - pad, all_pts[:, 0].max() + pad
    y_min, y_max = all_pts[:, 1].min() - pad, all_pts[:, 1].max() + pad
    dem_x = np.arange(x_min, x_max + spacing, spacing)
    dem_y = np.arange(y_min, y_max + spacing, spacing)
    dem_xx, dem_yy = np.meshgrid(dem_x, dem_y)
    dem_zz = griddata(all_pts[:, :2], all_pts[:, 2], (dem_xx, dem_yy), method="linear")
    mask_nan = np.isnan(dem_zz)
    if mask_nan.any():
        dem_zz_nn = griddata(all_pts[:, :2], all_pts[:, 2], (dem_xx, dem_yy), method="nearest")
        dem_zz[mask_nan] = dem_zz_nn[mask_nan]
    topo = np.c_[dem_xx.ravel(), dem_yy.ravel(), dem_zz.ravel()]
    print(f"  DEM: {len(topo)} points, {spacing}m spacing, Z=[{topo[:, 2].min():.1f}, {topo[:, 2].max():.1f}]")
    return topo


def build_tensor_mesh(xyz, topo, cell_size=CELL_SIZE, depth_core=DEPTH_CORE):
    """Build a 3D tensor mesh covering the data extent + topographic padding.

    Returns (mesh, active, n_active, nx, ny, nz).
    """
    from discretize import TensorMesh
    from discretize.utils import active_from_xyz

    x_min, x_max = xyz[:, 0].min(), xyz[:, 0].max()
    y_min, y_max = xyz[:, 1].min(), xyz[:, 1].max()
    z_min = topo[:, 2].min() - depth_core
    z_max = topo[:, 2].max()

    x0 = np.floor(x_min / cell_size) * cell_size - 2 * cell_size
    y0 = np.floor(y_min / cell_size) * cell_size - 2 * cell_size
    z_bot = np.floor(z_min / cell_size) * cell_size - 2 * cell_size
    z_top = np.ceil(z_max / cell_size) * cell_size + 2 * cell_size

    nx = int(np.ceil((x_max + 2 * cell_size - x0) / cell_size))
    ny = int(np.ceil((y_max + 2 * cell_size - y0) / cell_size))
    nz = int(np.ceil((z_top - z_bot) / cell_size))

    hx = np.ones(nx) * cell_size
    hy = np.ones(ny) * cell_size
    hz = np.ones(nz) * cell_size
    origin = np.array([x0, y0, z_bot])

    mesh = TensorMesh([hx, hy, hz], origin=origin)

    active = active_from_xyz(mesh, topo)
    n_act = int(active.sum())

    print(f"  Tensor mesh: {nx}x{ny}x{nz} = {mesh.n_cells} cells @ {cell_size}m")
    print(f"  Active cells: {n_act} / {mesh.n_cells}")

    return mesh, active, n_act, nx, ny, nz


# ── SimPEG inversions ──────────────────────────────────────────────────────


def run_simpeg_gravity_inversion(mesh, active, n_act, nx, ny, nz,
                                 grav_xyz, g_obs, output_dir):
    """Run SimPEG gravity-only inversion on a pre-built tensor mesh.

    Returns recovered density model (active cells only).
    """
    from simpeg import maps, data_misfit, regularization, optimization
    from simpeg import inverse_problem, inversion, directives, data
    from simpeg.potential_fields import gravity

    ConvergenceRecorder = _make_convergence_recorder()

    act_map = maps.InjectActiveCells(mesh, active, 0.0)

    grav_rx = gravity.receivers.Point(grav_xyz, components=["gz"])
    grav_src = gravity.sources.SourceField(receiver_list=[grav_rx])
    grav_survey = gravity.survey.Survey(grav_src)
    grav_sim = gravity.simulation.Simulation3DIntegral(
        mesh, survey=grav_survey, rhoMap=act_map,
        store_sensitivities="ram",
    )

    grav_data = data.Data(
        grav_survey, dobs=g_obs,
        standard_deviation=np.ones(len(g_obs)) * GRAV_UNCERTAINTY,
    )
    dmis = data_misfit.L2DataMisfit(data=grav_data, simulation=grav_sim)

    reg = regularization.WeightedLeastSquares(
        mesh, active_cells=active, alpha_z=2.0)

    m0 = np.ones(n_act) * 1e-4

    opt = optimization.ProjectedGNCG(
        maxIter=MAX_ITER, lower=-0.4, upper=0.6,
        tolF=1e-8, tolX=1e-8)

    inv_prob = inverse_problem.BaseInvProblem(dmis, reg, opt)

    conv_recorder = ConvergenceRecorder()
    dir_list = [
        directives.UpdateSensitivityWeights(every_iteration=False),
        directives.BetaEstimate_ByEig(beta0_ratio=10.0),
        directives.BetaSchedule(coolingFactor=2.0, coolingRate=5),
        directives.TargetMisfit(),
        conv_recorder,
    ]
    inv = inversion.BaseInversion(inv_prob, directiveList=dir_list)

    print(f"  Running gravity inversion (max_iter={MAX_ITER})...")
    t0 = time.time()
    rec_model = inv.run(m0)
    elapsed = time.time() - t0
    print(f"  Inversion completed in {elapsed:.1f}s")
    print(f"  Density range: [{rec_model.min():.4f}, {rec_model.max():.4f}] g/cc")

    # Save mesh
    os.makedirs(output_dir, exist_ok=True)
    msh_path = os.path.join(output_dir, "simpeg_mesh.msh")
    mesh.write_UBC(msh_path)
    print(f"  Wrote mesh to {msh_path}")

    # Save density model
    mod_path = os.path.join(output_dir, "simpeg_density.mod")
    _save_model_ubc(mod_path, rec_model, active, nx, ny, nz)
    print(f"  Wrote density model to {mod_path}")

    # Save convergence
    conv_path = os.path.join(output_dir, "simpeg_gravity_convergence.csv")
    conv_recorder.write_csv(conv_path)
    print(f"  Wrote convergence CSV to {conv_path}")

    return rec_model


def run_simpeg_magnetic_inversion(mesh, active, n_act, nx, ny, nz,
                                  mag_xyz, t_obs, output_dir):
    """Run SimPEG magnetic-only inversion on a pre-built tensor mesh.

    Uses TMI data with IGRF-directed inducing field. Returns recovered
    susceptibility model (active cells only, SI units).
    """
    from simpeg import maps, data_misfit, regularization, optimization
    from simpeg import inverse_problem, inversion, directives, data
    from simpeg.potential_fields import magnetics

    ConvergenceRecorder = _make_convergence_recorder()

    act_map = maps.InjectActiveCells(mesh, active, 0.0)

    mag_rx = magnetics.receivers.Point(mag_xyz, components=["tmi"])
    mag_src = magnetics.sources.UniformBackgroundField(
        receiver_list=[mag_rx],
        amplitude=IGRF_F,
        inclination=IGRF_I,
        declination=IGRF_D,
    )
    mag_survey = magnetics.survey.Survey(mag_src)
    mag_sim = magnetics.simulation.Simulation3DIntegral(
        mesh, survey=mag_survey, chiMap=act_map,
        store_sensitivities="ram",
    )

    mag_data = data.Data(
        mag_survey, dobs=t_obs,
        standard_deviation=np.ones(len(t_obs)) * MAG_UNCERTAINTY,
    )
    dmis = data_misfit.L2DataMisfit(data=mag_data, simulation=mag_sim)

    reg = regularization.WeightedLeastSquares(
        mesh, active_cells=active, alpha_z=2.0)

    m0 = np.ones(n_act) * 1e-4

    opt = optimization.ProjectedGNCG(
        maxIter=MAX_ITER, lower=0.0, upper=0.5,
        tolF=1e-8, tolX=1e-8)

    inv_prob = inverse_problem.BaseInvProblem(dmis, reg, opt)

    conv_recorder = ConvergenceRecorder()
    dir_list = [
        directives.UpdateSensitivityWeights(every_iteration=False),
        directives.BetaEstimate_ByEig(beta0_ratio=10.0),
        directives.BetaSchedule(coolingFactor=2.0, coolingRate=5),
        directives.TargetMisfit(),
        conv_recorder,
    ]
    inv = inversion.BaseInversion(inv_prob, directiveList=dir_list)

    print(f"  Running magnetic inversion (max_iter={MAX_ITER})...")
    t0 = time.time()
    rec_susc = inv.run(m0)
    elapsed = time.time() - t0
    print(f"  Inversion completed in {elapsed:.1f}s")
    print(f"  Susceptibility range: [{rec_susc.min():.6f}, {rec_susc.max():.6f}] SI")

    # Save susceptibility model
    mod_path = os.path.join(output_dir, "simpeg_susceptibility.mod")
    _save_model_ubc(mod_path, rec_susc, active, nx, ny, nz)
    print(f"  Wrote susceptibility model to {mod_path}")

    # Save convergence
    conv_path = os.path.join(output_dir, "simpeg_magnetic_convergence.csv")
    conv_recorder.write_csv(conv_path)
    print(f"  Wrote convergence CSV to {conv_path}")

    return rec_susc


def simpeg_to_meshdata(simpeg_mesh, active, density, n_act):
    """Convert SimPEG TensorMesh + active cells to MeshData."""
    from cluster_api._io import MeshData

    dx = float(simpeg_mesh.h[0][0])
    dy = float(simpeg_mesh.h[1][0])
    dz = float(simpeg_mesh.h[2][0])
    nx = len(simpeg_mesh.h[0])
    ny = len(simpeg_mesh.h[1])
    nz = len(simpeg_mesh.h[2])

    origin = simpeg_mesh.origin
    x0 = float(origin[0])
    y0 = float(origin[1])
    z0 = float(origin[2]) + nz * dz

    cc = simpeg_mesh.cell_centers[active]
    x_center = cc[:, 0]
    y_center = cc[:, 1]
    z_center = cc[:, 2]

    ix = np.round((x_center - x0) / dx - 0.5).astype(int) + 1
    iy = np.round((y_center - y0) / dy - 0.5).astype(int) + 1
    iz_from_bottom = np.round((z_center - float(origin[2])) / dz - 0.5).astype(int)
    iz = nz - iz_from_bottom

    mesh = MeshData(
        nx=nx, ny=ny, nz=nz,
        x0=x0, y0=y0, z0=z0,
        dx=dx, dy=dy, dz=dz,
        n_active=n_act,
        ix=ix, iy=iy, iz=iz,
        x_center=x_center, y_center=y_center, z_center=z_center,
    )
    return mesh


def validate_outputs(result):
    """Check that all expected output files exist and are non-empty."""
    print("\n=== Validation ===")
    checks_passed = 0
    checks_total = 0

    checks_total += 1
    if result.n_contacts > 0:
        print(f"  [PASS] {result.n_contacts} contact surface(s) extracted")
        checks_passed += 1
    else:
        print(f"  [FAIL] No contact surfaces extracted")

    checks_total += 1
    if os.path.isfile(result.ini_path) and os.path.getsize(result.ini_path) > 0:
        print(f"  [PASS] INI file: {result.ini_path}")
        checks_passed += 1
    else:
        print(f"  [FAIL] INI file missing or empty")

    for cs in result.contacts:
        name = f"contact_{cs.group_above}_{cs.group_below}"
        for ext in ("ts", "obj"):
            path = os.path.join(OUTPUT_DIR, "meshes", f"{name}.{ext}")
            checks_total += 1
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                print(f"  [PASS] {ext.upper()}: {name}")
                checks_passed += 1
            else:
                print(f"  [FAIL] Missing: {path}")

    for fname, label in [("cluster_properties.csv", "Cluster CSV"),
                          ("observed_gravity.csv", "Gravity CSV"),
                          ("observed_magnetic.csv", "Magnetic CSV"),
                          ("simpeg_mesh.msh", "SimPEG mesh"),
                          ("simpeg_density.mod", "SimPEG density model"),
                          ("simpeg_susceptibility.mod", "SimPEG susceptibility model"),
                          ("simpeg_gravity_convergence.csv", "Gravity convergence"),
                          ("simpeg_magnetic_convergence.csv", "Magnetic convergence"),
                          ("cluster_scatter_density.png", "Density scatter plot"),
                          ("cluster_scatter_susc.png", "Susc scatter plot"),
                          ("cluster_scatter_intersection.png", "Intersection scatter plot"),
                          ("density_clusters", "Density-only cluster volumes"),
                          ("susc_clusters", "Susc-only cluster volumes")]:
        path = os.path.join(OUTPUT_DIR, fname)
        checks_total += 1
        if os.path.exists(path):
            print(f"  [PASS] {label}: {path}")
            checks_passed += 1
        else:
            print(f"  [FAIL] Missing {label}")

    print(f"\n  Result: {checks_passed}/{checks_total} checks passed")
    return checks_passed == checks_total


def compute_local_origin(grav_xyz):
    """Compute local origin from data center."""
    ox = (grav_xyz[:, 0].min() + grav_xyz[:, 0].max()) / 2.0
    oy = (grav_xyz[:, 1].min() + grav_xyz[:, 1].max()) / 2.0
    return ox, oy


def write_local_csvs(grav_xyz, g_obs, mag_xyz, t_obs, local_origin, z_datum, output_dir):
    """Write gravity and magnetic CSVs in local coordinates for C++ step."""
    ox, oy = local_origin
    os.makedirs(output_dir, exist_ok=True)

    # Gravity
    grav_path = os.path.join(output_dir, "gravity_local.csv")
    with open(grav_path, "w") as f:
        f.write("x,y,z,g_obs\n")
        for i in range(len(g_obs)):
            x = grav_xyz[i, 0] - ox
            y = grav_xyz[i, 1] - oy
            z = grav_xyz[i, 2] - z_datum
            f.write(f"{x:.6f},{y:.6f},{z:.6f},{g_obs[i]:.6f}\n")
    print(f"  Wrote {grav_path} ({len(g_obs)} points)")

    # Magnetic
    mag_path = os.path.join(output_dir, "magnetics_local.csv")
    with open(mag_path, "w") as f:
        f.write("x,y,z,t_obs\n")
        for i in range(len(t_obs)):
            x = mag_xyz[i, 0] - ox
            y = mag_xyz[i, 1] - oy
            z = mag_xyz[i, 2] - z_datum
            f.write(f"{x:.6f},{y:.6f},{z:.6f},{t_obs[i]:.6f}\n")
    print(f"  Wrote {mag_path} ({len(t_obs)} points)")

    return grav_path, mag_path


def plot_cluster_scatter(density, susceptibility, labels, summary, path):
    """Generate density vs susceptibility scatter plot colored by cluster."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = labels >= 0
    d_val = density[valid]
    s_val = susceptibility[valid]
    l_val = labels[valid]

    if len(d_val) == 0:
        print("  No valid cells for scatter plot")
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    n_clusters = summary["n_groups"] if hasattr(summary, '__contains__') and "n_groups" in summary else len(set(l_val))
    unique_labels = sorted(set(l_val))
    cmap = plt.cm.tab10
    colors = [cmap(i % 10) for i in range(len(unique_labels))]

    for i, label in enumerate(unique_labels):
        mask = l_val == label
        ax.scatter(d_val[mask], s_val[mask], c=[colors[i]], alpha=0.5,
                   s=1, label=f"Group {label}")

    # Plot cluster centers
    if "density_mean" in summary and "susc_mean" in summary:
        for i in range(len(summary["density_mean"])):
            gid = summary["group_id"][i]
            ax.scatter(summary["density_mean"][i], summary["susc_mean"][i],
                       c="black", marker="x", s=100, linewidths=2,
                       zorder=5)

    ax.set_xlabel("Density (g/cm$^3$)")
    ax.set_ylabel("Magnetic Susceptibility (SI)")
    ax.set_title("GMM Clustering: Density vs Magnetic Susceptibility")
    ax.legend(markerscale=5, fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Wrote cluster scatter plot to {path}")


def plot_single_property_scatter(prop_vals, labels, summary, path, xlabel):
    """Generate single-property scatter (histogram + cell-index scatter)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = labels >= 0
    vals = prop_vals[valid]
    l_val = labels[valid]

    if len(vals) == 0:
        print(f"  No valid cells for {path}")
        return

    unique_labels = sorted(set(l_val))
    cmap = plt.cm.tab10
    colors = [cmap(i % 10) for i in range(len(unique_labels))]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram per cluster
    for i, label in enumerate(unique_labels):
        mask = l_val == label
        ax1.hist(vals[mask], bins=50, alpha=0.6, color=colors[i],
                 label=f"Group {label}")
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("Count")
    ax1.set_title("Per-Cluster Histogram")
    ax1.legend(fontsize=7)
    ax1.grid(True, alpha=0.3)

    # Scatter: cell index vs value
    indices = np.arange(len(vals))
    for i, label in enumerate(unique_labels):
        mask = l_val == label
        ax2.scatter(indices[mask], vals[mask], c=[colors[i]], alpha=0.3,
                    s=1, label=f"Group {label}")
    ax2.set_xlabel("Cell Index")
    ax2.set_ylabel(xlabel)
    ax2.set_title("Per-Cluster Scatter (by cell index)")
    ax2.legend(markerscale=5, fontsize=7)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Wrote {path}")


def main():
    print("=" * 70)
    print("FORRESTANIA E2E PIPELINE (SimPEG Gravity + Magnetic + Intersection Clustering)")
    print("=" * 70)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    grav_csv = os.path.join(DATASETS, "Forrestania_Gravity_Station_trim_.csv")
    mag_csv = os.path.join(DATASETS, "60472_AOI4_regional_mag_Zfixed.csv")

    if not os.path.isfile(grav_csv):
        print(f"ERROR: Gravity data not found: {grav_csv}")
        return 1
    if not os.path.isfile(mag_csv):
        print(f"ERROR: Magnetic data not found: {mag_csv}")
        return 1

    # ── Data preparation ────────────────────────────────────────────────
    print("\n[1/9] Preparing gravity data...")
    grav_xyz, g_obs = prep_gravity(grav_csv)

    print("\n[2/9] Preparing magnetic data...")
    mag_xyz, t_obs = prep_magnetics(mag_csv, n_target=500)

    print("\n[3/9] Building DEM...")
    topo = build_dem(grav_xyz, mag_xyz)

    local_origin = compute_local_origin(grav_xyz)
    z_datum = float(grav_xyz[:, 2].mean())
    print(f"  Local origin: ({local_origin[0]:.1f}, {local_origin[1]:.1f}), z_datum={z_datum:.1f}")

    # ── Tensor mesh (shared by both inversions) ─────────────────────────
    print("\n[4/9] Building tensor mesh...")
    simpeg_mesh, active, n_act, nx, ny, nz = build_tensor_mesh(
        grav_xyz, topo, CELL_SIZE, DEPTH_CORE,
    )

    # ── Separate SimPEG inversions ──────────────────────────────────────
    print("\n[5/9] Running SimPEG gravity inversion...")
    rec_dens = run_simpeg_gravity_inversion(
        simpeg_mesh, active, n_act, nx, ny, nz,
        grav_xyz, g_obs, OUTPUT_DIR,
    )

    print("\n[6/9] Running SimPEG magnetic inversion...")
    rec_susc = run_simpeg_magnetic_inversion(
        simpeg_mesh, active, n_act, nx, ny, nz,
        mag_xyz, t_obs, OUTPUT_DIR,
    )

    # ── Intersection clustering + contact extraction ────────────────────
    print("\n[7/9] Intersection clustering + contact surface extraction...")
    from cluster_api._cluster import cluster_intersection, cluster_lithology, LithologySummary
    from cluster_api._io import MeshData

    labels, summary = cluster_intersection(
        rec_dens, rec_susc,
        n_density=N_CLUSTERS, n_susc=3,
        random_state=42,
    )
    print(f"  Intersection: {summary.n_groups} groups"
          f" ({N_CLUSTERS} density x 3 susceptibility ->"
          f" {summary.n_groups} non-empty)")

    mesh_data = simpeg_to_meshdata(simpeg_mesh, active, rec_dens, n_act)
    print(f"  MeshData: {mesh_data.nx}x{mesh_data.ny}x{mesh_data.nz}, {mesh_data.n_active} active")

    # Use run_extract_from_labels to skip internal clustering
    from lithoseed._pipeline import run_extract_from_labels

    result = run_extract_from_labels(
        mesh_data,
        rec_dens,
        rec_susc,
        labels,
        summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=OUTPUT_DIR,
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
        print(f"    Contact {cs.group_above}-{cs.group_below}: {len(cs.vertices)} verts, depth={cs.median_depth:.0f}m")

    # ── Group closed volumes (voxel-face extraction: zero-gap, shared edges) ─
    from lithoseed._pipeline import run_extract_group_volumes

    # Raw: keep all fragments for debugging
    vol_raw_dir = os.path.join(OUTPUT_DIR, "volumes_voxel")
    vol_raw_result = run_extract_group_volumes(
        mesh_data, labels, summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=vol_raw_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=0,
    )

    # Cleaned: reassign fragments < 25 cells to neighbouring groups
    vol_clean_dir = os.path.join(OUTPUT_DIR, "volumes_voxel_cleaned")
    vol_clean_result = run_extract_group_volumes(
        mesh_data, labels, summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=vol_clean_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=25,
    )

    # ── Single-property cluster volumes (density-only + susceptibility-only) ─
    print("\n[8/9] Single-property clustering + volume extraction + scatter plots...")

    zero_susc = np.zeros_like(rec_dens)
    zero_dens = np.zeros_like(rec_susc)

    # Density-only clustering
    dens_labels, dens_summary = cluster_lithology(
        rec_dens, zero_susc, n_clusters=N_CLUSTERS, random_state=42)
    n_dens_groups = len(dens_summary.labels)
    print(f"  Density-only: {n_dens_groups} groups")
    for i in range(n_dens_groups):
        print(f"    Group {dens_summary.labels[i]}:"
              f" mean={dens_summary.density_mean[i]:.4f} +/- {dens_summary.density_std[i]:.4f},"
              f" ncells={dens_summary.counts[i]}")

    dens_ls = LithologySummary(
        labels=dens_summary.labels,
        counts=dens_summary.counts,
        density_mean=dens_summary.density_mean,
        density_std=dens_summary.density_std,
        susc_mean=[0.0] * n_dens_groups,
        susc_std=[0.0] * n_dens_groups,
    )

    dens_vol_dir = os.path.join(OUTPUT_DIR, "density_clusters")
    run_extract_group_volumes(
        mesh_data, dens_labels, dens_ls,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=dens_vol_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=25,
        igrf_f=IGRF_F, igrf_i=IGRF_I, igrf_d=IGRF_D,
    )
    print(f"  Density volumes: {dens_vol_dir}")

    # Density-only scatter (histogram + cell-index scatter)
    dens_scatter_path = os.path.join(OUTPUT_DIR, "cluster_scatter_density.png")
    plot_single_property_scatter(rec_dens, dens_labels, dens_summary,
                                 dens_scatter_path, xlabel="Density (g/cm³)")

    # Susceptibility-only clustering
    susc_labels, susc_summary = cluster_lithology(
        zero_dens, rec_susc, n_clusters=N_CLUSTERS, random_state=42)
    n_susc_groups = len(susc_summary.labels)
    print(f"  Susceptibility-only: {n_susc_groups} groups")
    for i in range(n_susc_groups):
        print(f"    Group {susc_summary.labels[i]}:"
              f" mean={susc_summary.susc_mean[i]:.6f} +/- {susc_summary.susc_std[i]:.6f},"
              f" ncells={susc_summary.counts[i]}")

    susc_ls = LithologySummary(
        labels=susc_summary.labels,
        counts=susc_summary.counts,
        density_mean=[0.0] * n_susc_groups,
        density_std=[0.0] * n_susc_groups,
        susc_mean=susc_summary.susc_mean,
        susc_std=susc_summary.susc_std,
    )

    susc_vol_dir = os.path.join(OUTPUT_DIR, "susc_clusters")
    run_extract_group_volumes(
        mesh_data, susc_labels, susc_ls,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=susc_vol_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=25,
        igrf_f=IGRF_F, igrf_i=IGRF_I, igrf_d=IGRF_D,
    )
    print(f"  Susceptibility volumes: {susc_vol_dir}")

    # Susceptibility-only scatter (histogram + cell-index scatter)
    susc_scatter_path = os.path.join(OUTPUT_DIR, "cluster_scatter_susc.png")
    plot_single_property_scatter(rec_susc, susc_labels, susc_summary,
                                 susc_scatter_path, xlabel="Magnetic Susceptibility (SI)")

    # ── Local-coordinate CSVs for C++ step ──────────────────────────────
    print("\n[9/9] Writing local-coordinate CSVs + intersection scatter...")
    grav_local_path, mag_local_path = write_local_csvs(
        grav_xyz, g_obs, mag_xyz, t_obs, local_origin, z_datum, OUTPUT_DIR,
    )

    # Cluster properties table
    if result.cluster_summary:
        print("\n  Cluster properties:")
        for i, gid in enumerate(result.cluster_summary.get("group_id", [])):
            d_mean = result.cluster_summary["density_mean"][i]
            d_std = result.cluster_summary["density_std"][i]
            s_mean = result.cluster_summary.get("susc_mean", [0.0] * len(result.cluster_summary["group_id"]))[i]
            s_std = result.cluster_summary.get("susc_std", [0.0] * len(result.cluster_summary["group_id"]))[i]
            print(f"    Group {gid}: density={d_mean:.4f} +/- {d_std:.4f}, susc={s_mean:.6f} +/- {s_std:.6f}")

        # Intersection scatter plot (density vs susceptibility)
        scatter_inter_path = os.path.join(OUTPUT_DIR, "cluster_scatter_intersection.png")
        plot_cluster_scatter(rec_dens, rec_susc, result.labels,
                             result.cluster_summary, scatter_inter_path)

    success = validate_outputs(result)

    print("\n" + "=" * 70)
    if success:
        print("FORRESTANIA E2E: ALL CHECKS PASSED")
        print(f"\nNext: Run C++ litho inversion using:")
        print(f"  config = {result.ini_path}")
        print(f"  gravity_local_csv = {grav_local_path}")
        print(f"  magnetics_local_csv = {mag_local_path}")
    else:
        print("FORRESTANIA E2E: SOME CHECKS FAILED")
    print("=" * 70)

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())

