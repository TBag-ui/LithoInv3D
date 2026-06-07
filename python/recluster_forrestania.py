#!/usr/bin/env python3
"""
Re-cluster existing SimPEG models and regenerate volume output with fixed
density conversion (contrast → absolute + scaling).

Skips the SimPEG inversion steps (which take hours). Uses the existing
simpeg_density.mod, simpeg_susceptibility.mod, and simpeg_mesh.msh files.
"""
import os, sys
import numpy as np

PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))

OUTPUT_DIR = os.path.join(PROJ_ROOT, "..", "apps", "forrestania", "inversion_output")
FORRESTANIA_APP = os.path.join(PROJ_ROOT, "..", "apps", "forrestania")

IGRF_F = 58874.0
IGRF_I = -66.2
IGRF_D = -0.1
N_CLUSTERS = 4


def load_ubc_mod(path):
    """Load a UBC-format .mod file (full grid, nx*ny*nz values)."""
    with open(path) as f:
        nx, ny, nz = map(int, f.readline().split())
    data = np.loadtxt(path, skiprows=1)
    return data, nx, ny, nz


def load_ubc_msh(path):
    """Load a UBC-format .msh file, returning (nx, ny, nz, x0, y0, z0, dx, dy, dz)."""
    with open(path) as f:
        nx, ny, nz = map(int, f.readline().split())
        x0, y0, z0 = map(float, f.readline().split())
        hx = np.fromstring(f.readline(), sep=" ")
        hy = np.fromstring(f.readline(), sep=" ")
        hz = np.fromstring(f.readline(), sep=" ")
    dx = float(hx[0]) if len(hx) > 0 else 200.0
    dy = float(hy[0]) if len(hy) > 0 else 200.0
    dz = float(hz[0]) if len(hz) > 0 else 200.0
    return nx, ny, nz, x0, y0, z0, dx, dy, dz


def build_meshdata(nx, ny, nz, x0, y0, z0, dx, dy, dz, active):
    """Build MeshData from active-cell mask on a full grid."""
    from cluster_api._io import MeshData

    n_active = int(active.sum())
    n_total = nx * ny * nz

    full_indices = np.arange(n_total)[active]
    ix = (full_indices % nx) + 1
    iy = ((full_indices // nx) % ny) + 1
    iz = (full_indices // (nx * ny)) + 1

    x_center = x0 + (ix - 0.5) * dx
    y_center = y0 + (iy - 0.5) * dy
    z_center = z0 - (iz - 0.5) * dz

    return MeshData(
        nx=nx, ny=ny, nz=nz,
        x0=x0, y0=y0, z0=z0,
        dx=dx, dy=dy, dz=dz,
        n_active=n_active,
        ix=ix.astype(np.int32),
        iy=iy.astype(np.int32),
        iz=iz.astype(np.int32),
        x_center=x_center,
        y_center=y_center,
        z_center=z_center,
    )


def main():
    print("=" * 70)
    print("FORRESTANIA RE-CLUSTER (skip SimPEG, use existing .mod files)")
    print("=" * 70)

    # --- Load existing models ---
    density_path = os.path.join(OUTPUT_DIR, "simpeg_density.mod")
    susc_path = os.path.join(OUTPUT_DIR, "simpeg_susceptibility.mod")
    msh_path = os.path.join(OUTPUT_DIR, "simpeg_mesh.msh")

    print(f"\nLoading density model: {density_path}")
    full_density, nx_d, ny_d, nz_d = load_ubc_mod(density_path)
    print(f"  Grid: {nx_d}x{ny_d}x{nz_d} = {len(full_density)} cells")
    print(f"  Range: [{full_density.min():.6f}, {full_density.max():.6f}]")

    print(f"\nLoading susceptibility model: {susc_path}")
    full_susc, nx_s, ny_s, nz_s = load_ubc_mod(susc_path)
    print(f"  Grid: {nx_s}x{ny_s}x{nz_s} = {len(full_susc)} cells")
    print(f"  Range: [{full_susc.min():.6f}, {full_susc.max():.6f}]")

    print(f"\nLoading mesh: {msh_path}")
    nx, ny, nz, x0, y0, z0, dx, dy, dz = load_ubc_msh(msh_path)
    print(f"  Mesh: {nx}x{ny}x{nz}, origin=({x0:.1f}, {y0:.1f}, {z0:.1f}),"
          f" cell=({dx}, {dy}, {dz})")

    assert nx == nx_d and ny == ny_d and nz == nz_d, \
        f"Mesh dimensions ({nx},{ny},{nz}) != density dimensions ({nx_d},{ny_d},{nz_d})"

    # --- Build active mask ---
    # In UBC format from _save_model_ubc: inactive cells = 0.0 exactly
    active = np.abs(full_density) > 1e-12
    print(f"\nActive cells: {active.sum()} / {len(full_density)} "
          f"({100 * active.sum() / len(full_density):.1f}%)")

    rec_dens = full_density[active].copy()
    rec_susc = full_susc[active].copy()

    # --- Build MeshData ---
    mesh_data = build_meshdata(nx, ny, nz, x0, y0, z0, dx, dy, dz, active)
    print(f"  MeshData: {mesh_data.nx}x{mesh_data.ny}x{mesh_data.nz}, "
          f"{mesh_data.n_active} active")

    # --- Intersection clustering ---
    print("\n--- Intersection clustering ---")
    from cluster_api._cluster import cluster_intersection
    labels, summary = cluster_intersection(
        rec_dens, rec_susc,
        n_density=N_CLUSTERS, n_susc=3,
        random_state=42,
    )
    print(f"  {summary.n_groups} groups ({N_CLUSTERS} density x 3 susceptibility -> "
          f"{summary.n_groups} non-empty)")

    # --- Volume extraction (the production path, with density fix) ---
    print("\n--- Volume extraction (volumes_voxel_cleaned) ---")
    from lithoseed._pipeline import run_extract_group_volumes

    # Compute local origin from data extents
    local_origin = (
        float(x0 + 0.5 * nx * dx),
        float(y0 + 0.5 * ny * dy),
        0.0,
    )
    z_datum = float(z0)
    print(f"  Local origin: ({local_origin[0]:.1f}, {local_origin[1]:.1f}), z_datum={z_datum:.1f}")

    # Cleaned: reassign fragments < 25 cells to neighbouring groups
    vol_clean_dir = os.path.join(OUTPUT_DIR, "volumes_voxel_cleaned")
    vol_clean_result = run_extract_group_volumes(
        mesh_data, labels, summary,
        local_origin=local_origin,
        z_datum=z_datum,
        output_dir=vol_clean_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=25,
        igrf_f=IGRF_F, igrf_i=IGRF_I, igrf_d=IGRF_D,
    )

    # --- Print summary ---
    csv_path = os.path.join(vol_clean_dir, "cluster_properties.csv")
    print(f"\n--- Cluster CSV: {csv_path} ---")
    with open(csv_path) as f:
        print(f.read())

    print("\n" + "=" * 70)
    print("RE-CLUSTER COMPLETE")
    print(f"Output: {vol_clean_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()

