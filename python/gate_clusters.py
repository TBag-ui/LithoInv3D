#!/usr/bin/env python3
"""
Extract closed volumes by manually gating physical-property ranges.

Replaces GMM clustering with user-defined property ranges.  Load a SimPEG
mesh + property files, define per-cluster density/susceptibility windows,
and extract closed volumes via the existing LithoSeed volume pipeline.

Usage:
    python gate_clusters.py \
        --mesh simpeg_mesh.msh \
        --density simpeg_density.mod \
        --susc simpeg_susceptibility.mod \
        --clusters my_gates.csv \
        --background-density 2.67 \
        --output-dir gated_output/
"""

import argparse, csv, os, sys
from typing import List
import numpy as np

PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))

from cluster_api._io import MeshData
from cluster_api._cluster import LithologySummary, cluster_lithology
from lithoseed._pipeline import run_extract_group_volumes


def _load_ubc_msh(path):
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


def _load_ubc_mod(path):
    """Load a UBC-format .mod file (full grid, nx*ny*nz values)."""
    with open(path) as f:
        first = f.readline().split()
    skip = 1 if len(first) == 3 else 0
    data = np.loadtxt(path, skiprows=skip)
    return data


def _build_meshdata(nx, ny, nz, x0, y0, z0, dx, dy, dz, active):
    """Build MeshData from an active-cell mask on a full grid."""
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


def read_cluster_gates(csv_path: str) -> List[dict]:
    """Read cluster definition CSV.

    Required columns: cluster_name, density_min, density_max.
    Optional columns: susc_min, susc_max (default: -inf, +inf).

    Returns list of dicts with keys: name, d_min, d_max, s_min, s_max.
    """
    rows = []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = [fn.strip() for fn in (reader.fieldnames or [])]

        required = {"cluster_name", "density_min", "density_max"}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(
                f"Cluster CSV missing required columns: {', '.join(sorted(missing))}"
            )

        has_susc = "susc_min" in fieldnames and "susc_max" in fieldnames

        for r in reader:
            name = r["cluster_name"].strip()
            d_min = float(r["density_min"])
            d_max = float(r["density_max"])
            s_min = float(r.get("susc_min", -np.inf))
            s_max = float(r.get("susc_max", np.inf))
            rows.append({
                "name": name,
                "d_min": d_min,
                "d_max": d_max,
                "s_min": s_min,
                "s_max": s_max,
                "has_susc": has_susc,
            })

    if not rows:
        raise ValueError(f"No cluster definitions found in {csv_path}")

    return rows


def gate_cells(
    density: np.ndarray,
    susceptibility: np.ndarray,
    cluster_defs: List[dict],
) -> np.ndarray:
    """Assign each cell to a cluster by property-range gating.

    Clusters are processed in definition order (first-match wins).
    Cells not matching any cluster get label -1.

    Returns (N,) int array of labels (-1 = unassigned, 1..N = cluster index).
    """
    labels = np.full(len(density), -1, dtype=np.int32)

    for ci, cdef in enumerate(cluster_defs):
        mask = (
            (density >= cdef["d_min"])
            & (density <= cdef["d_max"])
        )
        if cdef["has_susc"]:
            mask &= (
                (susceptibility >= cdef["s_min"])
                & (susceptibility <= cdef["s_max"])
            )
        # First-match: only assign cells not already claimed
        mask &= (labels < 0)
        labels[mask] = ci + 1  # 1-indexed for pipeline compatibility
        n_assigned = int(mask.sum())
        print(f"  [{ci + 1}] {cdef['name']}: density=[{cdef['d_min']:.3f}, {cdef['d_max']:.3f}]"
              f" -> {n_assigned} cells")

    return labels


def build_summary(
    labels: np.ndarray,
    density: np.ndarray,
    susceptibility: np.ndarray,
    cluster_defs: List[dict],
    do_susc: bool,
) -> LithologySummary:
    """Compute per-cluster property statistics for LithologySummary."""
    n_clusters = len(cluster_defs)
    summary_labels = list(range(1, n_clusters + 1))  # 1-indexed for pipeline
    counts = []
    d_mean_list = []
    d_std_list = []
    s_mean_list = []
    s_std_list = []

    for ci in range(n_clusters):
        mask = labels == (ci + 1)  # labels are 1-indexed
        n = int(mask.sum())
        counts.append(n)
        if n > 0:
            d_mean_list.append(float(np.mean(density[mask])))
            d_std_list.append(float(np.std(density[mask])) if n > 1 else 0.0)
            if do_susc:
                s_mean_list.append(float(np.mean(susceptibility[mask])))
                s_std_list.append(float(np.std(susceptibility[mask])) if n > 1 else 0.0)
            else:
                s_mean_list.append(0.0)
                s_std_list.append(0.0)
        else:
            d_mean_list.append(0.0)
            d_std_list.append(0.0)
            s_mean_list.append(0.0)
            s_std_list.append(0.0)

    return LithologySummary(
        labels=summary_labels,
        counts=counts,
        density_mean=d_mean_list,
        density_std=d_std_list,
        susc_mean=s_mean_list,
        susc_std=s_std_list,
    )


def _tag_cluster_source(csv_path: str, cluster_defs: List[dict]) -> None:
    """Post-process cluster CSV: prefix working_name with user_ or auto_."""
    if not os.path.exists(csv_path):
        return

    fieldnames, rows = [], []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    for r in rows:
        try:
            cid = int(r["cluster_id"])
        except (KeyError, ValueError):
            continue
        if cid < len(cluster_defs):
            source = cluster_defs[cid].get("source", "user")
            name = cluster_defs[cid]["name"]
            r["working_name"] = f"{source}_{name}"

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Extract closed volumes by manually gating property ranges")
    parser.add_argument(
        "--mesh", required=True,
        help="SimPEG tensor mesh (.msh) file")
    parser.add_argument(
        "--density", required=True,
        help="Density model (.mod) file")
    parser.add_argument(
        "--susc", default=None,
        help="Susceptibility model (.mod) file (optional)")
    parser.add_argument(
        "--clusters", required=True,
        help="CSV with cluster definitions: cluster_name,density_min,density_max,susc_min,susc_max")
    parser.add_argument(
        "--output-dir", default="gated_output",
        help="Output directory (default: gated_output/)")
    parser.add_argument(
        "--local-origin", default=None,
        help="x,y coordinate origin (auto-computed from mesh if omitted)")
    parser.add_argument(
        "--clean-min-cells", type=int, default=25,
        help="Min cells per fragment, smaller ones get reassigned (default: 25)")
    parser.add_argument(
        "--background-density", type=float, default=0.0,
        help="Background density to add to contrast values (SimPEG default: 2.67)"
    )
    parser.add_argument(
        "--auto-cluster-remainder", type=int, default=0,
        help="Auto-cluster unassigned cells into N GMM groups (0 = leave unassigned)"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("MANUAL PROPERTY GATING")
    print("=" * 70)

    # ---- Load data ----
    print(f"\nLoading mesh: {args.mesh}")
    print(f"  Density: {args.density}")
    print(f"  Susceptibility: {args.susc or '(none)'}")
    nx, ny, nz, x0, y0, z0, dx, dy, dz = _load_ubc_msh(args.mesh)
    print(f"  Mesh: {nx}x{ny}x{nz}, origin=({x0:.1f}, {y0:.1f}, {z0:.1f}),"
          f" cell=({dx}, {dy}, {dz})")

    full_density = _load_ubc_mod(args.density)
    print(f"  Density range: [{full_density.min():.4f}, {full_density.max():.4f}]")

    active = np.abs(full_density) > 1e-12
    print(f"  Active cells: {active.sum()} / {len(full_density)}")

    density = full_density[active].copy()

    if args.susc:
        full_susc = _load_ubc_mod(args.susc)
        susceptibility = full_susc[active].copy()
        has_susc = True
        print(f"  Susceptibility range: [{susceptibility.min():.6f}, {susceptibility.max():.6f}]")
    else:
        susceptibility = np.full(len(density), np.nan)
        has_susc = False

    if args.background_density != 0.0:
        d_label = "contrast"
    else:
        d_label = "absolute"
    print(f"  Density range ({d_label}): [{density.min():.4f}, {density.max():.4f}]")

    mesh_data = _build_meshdata(nx, ny, nz, x0, y0, z0, dx, dy, dz, active)

    # ---- Read cluster gates ----
    print(f"\nCluster definitions: {args.clusters}")
    cluster_defs = read_cluster_gates(args.clusters)
    print(f"  {len(cluster_defs)} cluster(s) defined")

    # Convert user's absolute density gate ranges to contrast space.
    # The volume pipeline internally adds BACKGROUND_DENSITY=2.67 to
    # density_mean values, so we must gate in contrast space.
    if args.background_density != 0.0:
        print(f"  Converting gate ranges: absolute -> contrast"
              f" (subtract {args.background_density:.2f})")
        for cdef in cluster_defs:
            cdef["d_min"] -= args.background_density
            cdef["d_max"] -= args.background_density

    # ---- Gate cells ----
    print("\nGating cells:")
    # Mark user-defined clusters
    for cdef in cluster_defs:
        cdef["source"] = "user"
    labels = gate_cells(density, susceptibility, cluster_defs)

    assigned = int((labels > 0).sum())  # 1-indexed: >0 means assigned
    if assigned == 0:
        print("\nERROR: No cells matched any cluster definition. Check your ranges.")
        sys.exit(1)
    print(f"  Total assigned: {assigned} / {len(density)} cells")

    # ---- Auto-cluster remainder ----
    if args.auto_cluster_remainder > 0:
        unassigned_mask = labels < 0
        n_unassigned = int(unassigned_mask.sum())
        n_auto = args.auto_cluster_remainder
        if n_unassigned >= n_auto * 10:  # need at least 10 cells per cluster
            print(f"\nAuto-clustering {n_unassigned} remainder cells into"
                  f" {n_auto} GMM groups:")
            auto_labels, auto_summary = cluster_lithology(
                density[unassigned_mask],
                susceptibility[unassigned_mask],
                n_clusters=n_auto,
                sentinel=-9999.0,
            )
            n_user = len(cluster_defs)
            auto_offset = np.where(auto_labels > 0, auto_labels + n_user, -1)
            labels[unassigned_mask] = auto_offset

            for i in range(n_auto):
                cluster_defs.append({
                    "name": f"Auto_{i + 1}",
                    "d_min": 0.0, "d_max": 0.0,
                    "s_min": 0.0, "s_max": 0.0,
                    "has_susc": has_susc,
                    "source": "auto",
                })
            print(f"  Total clusters: {len(cluster_defs)}"
                  f" ({n_user} user + {n_auto} auto)")
        else:
            print(f"\n  Too few remainder cells ({n_unassigned}) for"
                  f" {n_auto} auto-clusters — skipping")
    else:
        unassigned = int((labels < 0).sum())
        if unassigned:
            print(f"  {unassigned} cells unassigned"
                  f" (use --auto-cluster-remainder N to fill gaps)")

    # ---- Build summary ----
    summary = build_summary(labels, density, susceptibility, cluster_defs, has_susc)
    print(f"\nPer-cluster statistics:")
    for i, cdef in enumerate(cluster_defs):
        ci = summary.labels[i] - 1
        tag = "[user]" if cdef.get("source") == "user" else "[auto]"
        print(f"  {tag} [{ci + 1}] {cdef['name']}: "
              f"d={summary.density_mean[ci]:.3f}±{summary.density_std[ci]:.3f}, "
              f"s={summary.susc_mean[ci]:.6f}±{summary.susc_std[ci]:.6f}, "
              f"n={summary.counts[ci]}")

    # ---- Local origin ----
    if args.local_origin:
        parts = [float(x) for x in args.local_origin.split(",")]
        local_origin = (parts[0], parts[1], 0.0)
    else:
        local_origin = (
            float(mesh_data.x0 + 0.5 * mesh_data.nx * mesh_data.dx),
            float(mesh_data.y0 + 0.5 * mesh_data.ny * mesh_data.dy),
            0.0,
        )
    z_datum = float(mesh_data.z0)
    print(f"\nLocal origin: ({local_origin[0]:.1f}, {local_origin[1]:.1f}), "
          f"z_datum={z_datum:.1f}")

    # ---- Extract volumes ----
    print(f"\nExtracting volumes -> {args.output_dir}")
    result = run_extract_group_volumes(
        mesh_data, labels, summary,
        local_origin=local_origin,
        z_datum=z_datum,
        output_dir=args.output_dir,
        export_formats=("ts", "obj"),
        clean_min_cells=args.clean_min_cells,
        has_susceptibility=has_susc,
        min_density_range=0.0,  # user already specified exact ranges
    )

    print("\n" + "=" * 70)
    print("GATING COMPLETE")

    # ---- Post-process CSV: tag working names with user/auto source ----
    csv_path = os.path.join(args.output_dir, "cluster_properties.csv")
    _tag_cluster_source(csv_path, cluster_defs)
    print(f"Cluster CSV: {csv_path}")
    print("=" * 70)


if __name__ == "__main__":
    main()

