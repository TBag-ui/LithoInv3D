"""Re-run volumetric extraction using existing SimPEG results.
Generates group_meshes INI + label_grid.bin for C++ volumetric loading path.
No SimPEG dependency — parses UBC mesh directly.
"""
import os
import sys
import numpy as np

PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))

from cluster_api._cluster import cluster_intersection
from cluster_api._io import MeshData
from lithoseed._pipeline import run_extract_group_volumes

OUTPUT_DIR = os.path.join(PROJ_ROOT, "..", "apps", "forrestania", "inversion_output")
DATASETS = os.path.join(PROJ_ROOT, "..", "datasets", "Forrestania")

IGRF_F = 58874.0
IGRF_I = -66.2
IGRF_D = -0.1
N_CLUSTERS = 4


def read_ubc_mesh(path):
    """Parse a UBC-format tensor mesh file (as written by SimPEG write_UBC).

    Returns (nx, ny, nz, x0, y0, z0, dx, dy, dz).
    """
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    nx, ny, nz = map(int, lines[0].split())
    x0, y0, z0 = map(float, lines[1].split())

    dx_vals = list(map(float, lines[2].split()))
    dy_vals = list(map(float, lines[3].split()))
    dz_vals = list(map(float, lines[4].split()))

    dx = dx_vals[0] if len(set(dx_vals)) == 1 else None
    dy = dy_vals[0] if len(set(dy_vals)) == 1 else None
    dz = dz_vals[0] if len(set(dz_vals)) == 1 else None

    if dx is None or dy is None or dz is None:
        raise ValueError("Non-uniform cell sizes not supported")

    return nx, ny, nz, x0, y0, z0, dx, dy, dz


def read_ubc_model(path):
    """Read a UBC-format model file (header + all nx*ny*nz values).

    Handles both one-value-per-line and np.savetxt multi-value-per-line formats.
    """
    with open(path) as f:
        text = f.read()

    # Split all whitespace-separated tokens
    tokens = text.split()
    nx, ny, nz = int(tokens[0]), int(tokens[1]), int(tokens[2])
    values = np.array([float(t) for t in tokens[3:]])
    expected = nx * ny * nz
    if len(values) != expected:
        raise ValueError(f"Model has {len(values)} values, expected {expected} ({nx}x{ny}x{nz})")
    return nx, ny, nz, values


def main():
    print("=" * 70)
    print("REGENERATING VOLUMETRIC INI + LABEL GRID")
    print("=" * 70)

    # Load existing SimPEG results
    print("\n[1/3] Loading existing SimPEG results...")
    msh_path = os.path.join(OUTPUT_DIR, "simpeg_mesh.msh")
    dens_path = os.path.join(OUTPUT_DIR, "simpeg_density.mod")
    susc_path = os.path.join(OUTPUT_DIR, "simpeg_susceptibility.mod")

    if not os.path.isfile(msh_path):
        print(f"ERROR: Mesh not found: {msh_path}")
        print("Run run_forrestania_e2e.py first to generate SimPEG results.")
        return 1

    nx, ny, nz, x0, y0, z0, dx, dy, dz = read_ubc_mesh(msh_path)
    print(f"  Mesh: {nx}x{ny}x{nz} cells @ {dx}m, origin=({x0:.1f}, {y0:.1f}, {z0:.1f})")

    _, _, _, rec_dens_full = read_ubc_model(dens_path)
    _, _, _, rec_susc_full = read_ubc_model(susc_path)
    print(f"  Density range: [{rec_dens_full.min():.4f}, {rec_dens_full.max():.4f}]")
    print(f"  Susceptibility range: [{rec_susc_full.min():.6f}, {rec_susc_full.max():.6f}]")

    # Determine active cells (model value != 0, since inactive = 0)
    active_mask = np.abs(rec_dens_full) > 1e-10
    n_act = int(active_mask.sum())
    print(f"  Active cells: {n_act} / {len(rec_dens_full)}")

    # Get active cell values
    rec_dens = rec_dens_full[active_mask]
    rec_susc = rec_susc_full[active_mask]

    # Compute 1-based 3D indices for active cells (UBC order: i fastest, j, k slowest)
    active_indices_flat = np.where(active_mask)[0]
    ix = np.zeros(n_act, dtype=int)
    iy = np.zeros(n_act, dtype=int)
    iz = np.zeros(n_act, dtype=int)
    for idx in range(n_act):
        flat_idx = active_indices_flat[idx]
        iz[idx] = flat_idx // (nx * ny) + 1      # 1-based
        iy[idx] = (flat_idx % (nx * ny)) // nx + 1
        ix[idx] = flat_idx % nx + 1

    # Build MeshData
    mesh_data = MeshData(
        nx=nx, ny=ny, nz=nz,
        x0=x0, y0=y0, z0=z0,
        dx=dx, dy=dy, dz=dz,
        n_active=n_act,
        ix=ix, iy=iy, iz=iz,
        x_center=x0 + (ix - 0.5) * dx,
        y_center=y0 + (iy - 0.5) * dy,
        z_center=z0 - (iz - 0.5) * dz,
    )
    print(f"  MeshData: {nx}x{ny}x{nz}, {n_act} active")

    # Intersection clustering (deterministic, same as e2e script)
    print("\n[2/3] Running intersection clustering...")
    labels, summary = cluster_intersection(
        rec_dens, rec_susc,
        n_density=N_CLUSTERS, n_susc=3,
        random_state=42,
    )
    print(f"  Groups: {summary.n_groups} ({N_CLUSTERS} density x 3 susceptibility)")

    # Load gravity + magnetic data for CSVs
    grav_csv = os.path.join(DATASETS, "Forrestania_Gravity_Station_trim_.csv")
    mag_csv = os.path.join(DATASETS, "60472_AOI4_regional_mag_Zfixed.csv")

    import pandas as pd
    grav_xyz = None
    g_obs = None
    mag_xyz = None
    t_obs = None
    z_datum = 0.0
    local_origin = (0.0, 0.0, 0.0)

    rng = np.random.default_rng(42)

    if os.path.isfile(grav_csv):
        gdf = pd.read_csv(grav_csv)
        grav_xyz = gdf[["X", "Y", "Z"]].values.astype(float)
        g_obs = -gdf["FGrav_mgal"].values.astype(float)  # sign flip per e2e prep_gravity
        z_datum = float(grav_xyz[:, 2].mean())
        local_origin = (float(grav_xyz[:, 0].min()), float(grav_xyz[:, 1].min()), 0.0)
        print(f"  Loaded {len(g_obs)} gravity points, g_obs range [{g_obs.min():.2f}, {g_obs.max():.2f}]")

    if os.path.isfile(mag_csv):
        mdf = pd.read_csv(mag_csv)
        mdf["RMI"] = mdf["MAGCOMP"] - mdf["IGRF"]
        n_target = min(500, len(mdf))
        idx = rng.choice(len(mdf), n_target, replace=False)
        sub = mdf.iloc[idx].copy()
        mag_xyz = sub[["X", "Y", "Z"]].values.astype(float)
        rmi = sub["RMI"].values.astype(float)
        # Plane detrend
        A = np.c_[mag_xyz[:, 0], mag_xyz[:, 1], np.ones(len(mag_xyz))]
        coeff, _, _, _ = np.linalg.lstsq(A, rmi, rcond=None)
        trend = A @ coeff
        t_obs = rmi - trend
        print(f"  Loaded {len(t_obs)} magnetic points, t_obs range [{t_obs.min():.1f}, {t_obs.max():.1f}]")

    # Volumetric extraction with updated code
    print("\n[3/3] Running volumetric extraction with updated pipeline...")

    vol_clean_dir = os.path.join(OUTPUT_DIR, "volumes_voxel_cleaned")
    result = run_extract_group_volumes(
        mesh_data, labels, summary,
        local_origin=(local_origin[0], local_origin[1], 0.0),
        z_datum=z_datum,
        output_dir=vol_clean_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=25,
        gravity_xyz=grav_xyz,
        gravity_obs=g_obs,
        magnetic_xyz=mag_xyz,
        magnetic_obs=t_obs,
        igrf_f=IGRF_F,
        igrf_i=IGRF_I,
        igrf_d=IGRF_D,
        has_susceptibility=True,
    )

    print(f"\n  Generated INI: {result.ini_path}")
    print(f"  Group meshes: {result.n_contacts}")
    print(f"  Groups: {result.group_order}")

    # Verify
    ini_path = result.ini_path
    label_grid_path = os.path.join(vol_clean_dir, "meshes", "label_grid.bin")

    if os.path.isfile(ini_path):
        with open(ini_path) as f:
            content = f.read()
        has_group = "group_meshes" in content
        has_label = "label_grid" in content
        print(f"  INI: group_meshes={'YES' if has_group else 'NO'}, label_grid={'YES' if has_label else 'NO'}")

    if os.path.isfile(label_grid_path):
        size_mb = os.path.getsize(label_grid_path) / (1024 * 1024)
        print(f"  Label grid: {size_mb:.2f} MB")

    print("\nDone. Ready for C++ volumetric inversion with group_meshes INI key.")
    return 0


if __name__ == "__main__":
    exit(main())

