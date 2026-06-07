#!/usr/bin/env python3
"""Generate separate gravity-only and magnetic-only cluster volumes + scatter plots.

Uses the existing SimPEG recovery models.  Does NOT re-run inversions.
Output dirs: inversion_output/density_clusters/  and  inversion_output/susceptibility_clusters/
"""
import os
import sys
import numpy as np

PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))

from cluster_api._cluster import cluster_lithology, LithologySummary
from cluster_api._io import MeshData
from lithoseed._pipeline import run_extract_group_volumes

OUTPUT_DIR = os.path.join(PROJ_ROOT, "..", "apps", "forrestania", "inversion_output")
DATASETS = os.path.join(PROJ_ROOT, "..", "datasets", "Forrestania")

N_CLUSTERS = 4


def read_ubc_mesh(path):
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    nx, ny, nz = map(int, lines[0].split())
    x0, y0, z0 = map(float, lines[1].split())
    dx_vals = list(map(float, lines[2].split()))
    dy_vals = list(map(float, lines[3].split()))
    dz_vals = list(map(float, lines[4].split()))
    dx = dx_vals[0]; dy = dy_vals[0]; dz = dz_vals[0]
    return nx, ny, nz, x0, y0, z0, dx, dy, dz


def read_ubc_model(path):
    with open(path) as f:
        text = f.read()
    tokens = text.split()
    nx, ny, nz = int(tokens[0]), int(tokens[1]), int(tokens[2])
    values = np.array([float(t) for t in tokens[3:]])
    return nx, ny, nz, values


def build_meshdata(nx, ny, nz, x0, y0, z0, dx, dy, dz, active_mask):
    n_act = int(active_mask.sum())
    active_indices_flat = np.where(active_mask)[0]
    ix = np.zeros(n_act, dtype=int)
    iy = np.zeros(n_act, dtype=int)
    iz = np.zeros(n_act, dtype=int)
    for idx in range(n_act):
        flat_idx = active_indices_flat[idx]
        iz[idx] = flat_idx // (nx * ny) + 1
        iy[idx] = (flat_idx % (nx * ny)) // nx + 1
        ix[idx] = flat_idx % nx + 1

    return MeshData(
        nx=nx, ny=ny, nz=nz,
        x0=x0, y0=y0, z0=z0,
        dx=dx, dy=dy, dz=dz,
        n_active=n_act,
        ix=ix, iy=iy, iz=iz,
        x_center=x0 + (ix - 0.5) * dx,
        y_center=y0 + (iy - 0.5) * dy,
        z_center=z0 - (iz - 0.5) * dz,
    )


def plot_scatter(property_vals, labels, summary, path, xlabel="Density (g/cm³)"):
    """Generate a scatter plot of property values colored by cluster."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = labels >= 0
    vals = property_vals[valid]
    l_val = labels[valid]

    if len(vals) == 0:
        print(f"  No valid cells for {path}")
        return

    unique_labels = sorted(set(l_val))
    n_groups = len(unique_labels)
    cmap = plt.cm.tab10
    colors = [cmap(i % 10) for i in range(n_groups)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram
    for i, label in enumerate(unique_labels):
        mask = l_val == label
        ax1.hist(vals[mask], bins=50, alpha=0.6, color=colors[i],
                 label=f"Group {label}")
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("Count")
    ax1.set_title("Per-Cluster Histogram")
    ax1.legend(fontsize=7)
    ax1.grid(True, alpha=0.3)

    # Index scatter (cell index vs value)
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


def run_property_extraction(property_name, prop_values, output_subdir):
    """Run clustering + volume extraction for a single property."""
    print(f"\n{'='*60}")
    print(f"  {property_name.upper()}-ONLY CLUSTERING ({N_CLUSTERS} clusters)")
    print(f"{'='*60}")

    # Cluster
    zero_susc = np.zeros_like(prop_values)
    labels, summary = cluster_lithology(
        prop_values, zero_susc,
        n_clusters=N_CLUSTERS,
        random_state=42,
    )

    print(f"  Clusters: {summary.n_groups}")
    for i in range(summary.n_groups):
        print(f"    Group {summary.labels[i]}:"
              f" mean={summary.density_mean[i]:.6f} +/- {summary.density_std[i]:.6f},"
              f" ncells={summary.counts[i]}")

    # Build a fake summary with the property in density_mean for volume extraction
    # (run_extract_group_volumes reads density_mean from summary)
    fake_summary = LithologySummary(
        labels=summary.labels,
        counts=summary.counts,
        density_mean=summary.density_mean,  # actual property values
        density_std=summary.density_std,
        susc_mean=[0.0] * summary.n_groups,
        susc_std=[0.0] * summary.n_groups,
    )

    output_dir = os.path.join(OUTPUT_DIR, output_subdir)
    result = run_extract_group_volumes(
        mesh_data, labels, fake_summary,
        local_origin=(local_origin_x, local_origin_y, 0.0),
        z_datum=z_datum,
        output_dir=output_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=25,
        igrf_f=58874.0, igrf_i=-66.2, igrf_d=-0.1,
    )

    print(f"  Output: {output_dir}")
    print(f"  Config: {result.ini_path}")
    return labels, summary


def main():
    global mesh_data, local_origin_x, local_origin_y, z_datum

    # ── Load SimPEG results ──────────────────────────────────────────
    msh_path = os.path.join(OUTPUT_DIR, "simpeg_mesh.msh")
    dens_path = os.path.join(OUTPUT_DIR, "simpeg_density.mod")
    susc_path = os.path.join(OUTPUT_DIR, "simpeg_susceptibility.mod")

    for p in [msh_path, dens_path, susc_path]:
        if not os.path.isfile(p):
            print(f"ERROR: {p} not found. Run run_forrestania_e2e.py first.")
            return 1

    nx, ny, nz, x0, y0, z0, dx, dy, dz = read_ubc_mesh(msh_path)
    _, _, _, rec_dens_full = read_ubc_model(dens_path)
    _, _, _, rec_susc_full = read_ubc_model(susc_path)

    active_mask = np.abs(rec_dens_full) > 1e-10
    n_act = int(active_mask.sum())
    print(f"Mesh: {nx}x{ny}x{nz}, {n_act} active cells")
    print(f"Density: [{rec_dens_full.min():.4f}, {rec_dens_full.max():.4f}]")
    print(f"Susceptibility: [{rec_susc_full.min():.6f}, {rec_susc_full.max():.6f}]")

    dens_active = rec_dens_full[active_mask]
    susc_active = rec_susc_full[active_mask]

    mesh_data = build_meshdata(nx, ny, nz, x0, y0, z0, dx, dy, dz, active_mask)

    # Local origin / datum from gravity data
    import pandas as pd
    global local_origin_x, local_origin_y, z_datum
    grav_csv = os.path.join(DATASETS, "Forrestania_Gravity_Station_trim_.csv")
    if os.path.isfile(grav_csv):
        gdf = pd.read_csv(grav_csv)
        gxyz = gdf[["X", "Y", "Z"]].values.astype(float)
        local_origin_x = (gxyz[:, 0].min() + gxyz[:, 0].max()) / 2.0
        local_origin_y = (gxyz[:, 1].min() + gxyz[:, 1].max()) / 2.0
        z_datum = float(gxyz[:, 2].mean())
    else:
        local_origin_x = 744682.28
        local_origin_y = 6415294.0
        z_datum = 381.232137

    # ── Density-only clusters + volumes ──────────────────────────────
    dens_labels, dens_summary = run_property_extraction(
        "density", dens_active, "density_clusters")

    # Density scatter plot
    scatter_path = os.path.join(OUTPUT_DIR, "cluster_scatter_density.png")
    plot_scatter(dens_active, dens_labels, dens_summary, scatter_path,
                 xlabel="Density (g/cm³)")

    # ── Susceptibility-only clusters + volumes ───────────────────────
    susc_labels, susc_summary = run_property_extraction(
        "susceptibility", susc_active, "susceptibility_clusters")

    # Susceptibility scatter plot
    scatter_path = os.path.join(OUTPUT_DIR, "cluster_scatter_susceptibility.png")
    plot_scatter(susc_active, susc_labels, susc_summary, scatter_path,
                 xlabel="Magnetic Susceptibility (SI)")

    print(f"\n{'='*60}")
    print("DONE. Outputs:")
    print(f"  Density clusters:  {os.path.join(OUTPUT_DIR, 'density_clusters')}")
    print(f"  Susc clusters:     {os.path.join(OUTPUT_DIR, 'susceptibility_clusters')}")
    print(f"  Scatter plots:     {os.path.join(OUTPUT_DIR, 'cluster_scatter_density.png')}")
    print(f"                     {os.path.join(OUTPUT_DIR, 'cluster_scatter_susceptibility.png')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

