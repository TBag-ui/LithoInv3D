#!/usr/bin/env python3
"""
Split multi-component volume meshes into separate files.

Detects disconnected components in .ts mesh files and splits them into
independent volumes.  Updates cluster_properties.csv with copied rows
so the C++ inversion can treat split components as separate groups.

Usage:
    python split_mesh_volumes.py volume_group_3.ts --csv cluster_properties.csv
    python split_mesh_volumes.py meshes/ --all --csv cluster_properties.csv --formats ts,obj
"""

import argparse, csv, os, re, sys
from typing import Optional, Tuple as TupleType
import numpy as np

PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))

from lithoseed._surfaces import ContactSurface, separate_connected_components
from lithoseed._export import write_volume_ts, write_volume_obj


def parse_ts(path: str) -> TupleType[np.ndarray, np.ndarray, str]:
    """Parse a GOCAD TSurf file, returning (vertices, faces, name).

    vertices: (N, 3) float array
    faces:    (M, 3) int array, 0-based indices
    name:     value from the HEADER name: field, or filename stem
    """
    verts = []  # list of (x,y,z) float triples
    faces = []  # list of (v0,v1,v2) int triples
    name = os.path.splitext(os.path.basename(path))[0]
    in_header = False

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == "GOCAD TSurf 1":
                continue
            if line == "HEADER {":
                in_header = True
                continue
            if in_header:
                if line == "}":
                    in_header = False
                elif line.startswith("name:"):
                    name = line.split(":", 1)[1].strip()
                continue
            if line == "END":
                break
            if line == "TFACE":
                continue
            if line.startswith("VRTX"):
                parts = line.split()
                verts.append((float(parts[2]), float(parts[3]), float(parts[4])))
            elif line.startswith("TRGL"):
                parts = line.split()
                # Convert 1-based → 0-based
                faces.append((int(parts[1]) - 1, int(parts[2]) - 1, int(parts[3]) - 1))

    return np.array(verts, dtype=np.float64), np.array(faces, dtype=np.int32), name


def _cluster_id_from_stem(stem: str) -> Optional[int]:
    """Extract cluster/label ID from a volume_group_N or csv_start_cluster_id_N stem."""
    m = re.search(r"(\d+)$", stem)
    return int(m.group(1)) if m else None


def _natural_sort_key(name: str) -> TupleType[int, int]:
    """Key function for natural sort: volume_group_2 < volume_group_11."""
    m = re.search(r"(\d+)$", os.path.splitext(name)[0])
    if m:
        return (0, int(m.group(1)))
    return (1, 0)


def _build_positional_mapping(mesh_files: list, csv_rows: list) -> dict:
    """Map each mesh filename to its cluster_id by position in sorted order.

    Mesh files sorted naturally (1,2,...,9,11,12) map to CSV rows 0..N-1 in order.
    Returns {filename: cluster_id} dict.
    """
    sorted_files = sorted(mesh_files, key=_natural_sort_key)
    mapping = {}
    for i, fname in enumerate(sorted_files):
        if i < len(csv_rows):
            mapping[fname] = int(csv_rows[i]["cluster_id"])
    return mapping


def _read_cluster_csv(path: str) -> TupleType[list, list]:
    """Read cluster_properties.csv, returning (fieldnames, rows)."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return fieldnames, rows


def _write_cluster_csv(path: str, fieldnames: list, rows: list) -> None:
    """Write cluster_properties.csv with copied rows for split components."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def split_mesh(
    mesh_path,
    csv_path,
    output_dir,
    formats=("ts",),
    override_cluster_id=None,
):
    """Split a single .ts mesh into its disconnected components.

    Returns True if the mesh was split (had >1 component).
    """
    stem = os.path.splitext(os.path.basename(mesh_path))[0]
    cluster_id = override_cluster_id if override_cluster_id is not None else _cluster_id_from_stem(stem)

    print(f"Loading: {mesh_path}")
    verts, faces, name = parse_ts(mesh_path)
    print(f"  {len(verts)} vertices, {len(faces)} faces")

    surface = ContactSurface(
        group_above=cluster_id or 0,
        group_below=0,
        vertices=verts,
        faces=faces,
        median_depth=0.0,
    )

    components = separate_connected_components(surface, min_faces=1)
    n_comp = len(components)

    if n_comp <= 1:
        print(f"  Single component — no split needed")
        # Copy to output dir if different from source
        if os.path.dirname(mesh_path) != output_dir:
            for fmt in formats:
                src = mesh_path if fmt == "ts" else mesh_path.replace(".ts", f".{fmt}")
                dst = os.path.join(output_dir, os.path.basename(src))
                if os.path.exists(src):
                    import shutil
                    shutil.copy2(src, dst)
        return False

    print(f"  {n_comp} components — splitting")
    for i, comp in enumerate(components):
        vol_name = f"{stem}_vol{i}"
        print(f"    {vol_name}: {len(comp.vertices)} verts, {len(comp.faces)} faces")
        comp.group_above = cluster_id or 0
        if "ts" in formats:
            ts_path = os.path.join(output_dir, f"{vol_name}.ts")
            write_volume_ts(ts_path, comp, i)
        if "obj" in formats:
            obj_path = os.path.join(output_dir, f"{vol_name}.obj")
            write_volume_obj(obj_path, comp, i)

    # Update CSV — copy rows for extra components
    if cluster_id is not None and os.path.exists(csv_path):
        fieldnames, rows = _read_cluster_csv(csv_path)
        matching = [r for r in rows if r.get("cluster_id") == str(cluster_id)]
        if matching:
            template = matching[0]
            max_id = max(int(r.get("cluster_id", 0)) for r in rows)
            for i in range(1, n_comp):
                new_row = dict(template)
                max_id += 1
                new_row["cluster_id"] = str(max_id)
                new_row["working_name"] = f"group_{max_id}"
                rows.append(new_row)
                print(f"    CSV: copied row cluster_id={max_id} from {cluster_id}")
            _write_cluster_csv(csv_path, fieldnames, rows)

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Split multi-component volume meshes into separate files")
    parser.add_argument(
        "mesh_path", help="Path to a .ts file, or directory with --all")
    parser.add_argument(
        "--csv", required=True, help="Path to cluster_properties.csv")
    parser.add_argument(
        "--all", action="store_true", dest="all_files",
        help="Process all .ts files in the given directory")
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory (default: same directory as input)")
    parser.add_argument(
        "--formats", default="ts",
        help="Export formats, comma-separated: ts,obj (default: ts)")
    args = parser.parse_args()

    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())

    if args.all_files:
        mesh_dir = args.mesh_path
        output_dir = args.output_dir or os.path.join(mesh_dir, "split_mesh_volumes")
        os.makedirs(output_dir, exist_ok=True)
        ts_files = sorted(
            f for f in os.listdir(mesh_dir)
            if f.endswith(".ts") and ("volume_group_" in f or "csv_start_cluster_id_" in f)
        )
        print(f"Processing {len(ts_files)} .ts files in {mesh_dir}")
        print(f"Output: {output_dir}")

        # Work from a copy of the CSV in the output directory
        csv_output = os.path.join(output_dir, os.path.basename(args.csv))
        import shutil
        shutil.copy2(args.csv, csv_output)

        # Build positional mapping: naturally-sorted mesh files → CSV rows
        _, rows = _read_cluster_csv(csv_output)
        fname_to_cluster = _build_positional_mapping(ts_files, rows)
        for fname, cid in fname_to_cluster.items():
            print(f"  Mapping: {fname} -> cluster_id={cid}")

        any_split = False
        for ts_file in ts_files:
            full_path = os.path.join(mesh_dir, ts_file)
            cid = fname_to_cluster.get(ts_file)
            if split_mesh(full_path, csv_output, output_dir, formats,
                          override_cluster_id=cid):
                any_split = True
            print()
        if not any_split:
            print("All meshes are single-component — no splits performed.")
        else:
            print(f"Updated CSV: {csv_output}")
    else:
        mesh_dir = os.path.dirname(args.mesh_path) or "."
        output_dir = args.output_dir or os.path.join(mesh_dir, "split_mesh_volumes")
        os.makedirs(output_dir, exist_ok=True)

        # Work from a copy of the CSV in the output directory
        csv_output = os.path.join(output_dir, os.path.basename(args.csv))
        import shutil
        shutil.copy2(args.csv, csv_output)

        split_mesh(args.mesh_path, csv_output, output_dir, formats)


if __name__ == "__main__":
    main()

