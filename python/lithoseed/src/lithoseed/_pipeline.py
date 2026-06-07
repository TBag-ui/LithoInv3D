"""Pipeline orchestrator for LithoSeed."""

import os
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional

from ._surfaces import (
    ContactSurface, extract_contact_surfaces,
    extract_group_volumes, export_label_grid,
    cleanup_label_fragments,
    decimate_mesh, enforce_stratigraphic_ordering, separate_connected_components,
    merge_volume_components,
)
from ._export import (
    write_contact_ts, write_contact_obj,
    write_volume_ts, write_volume_obj,
    write_cluster_csv,
    write_observed_gravity_csv, write_observed_magnetic_csv,
    write_borehole_constraints_csv,
)
from ._ini import write_ini_config


@dataclass
class PipelineResult:
    """Output of the LithoSeed pipeline."""
    contacts: List[ContactSurface] = field(default_factory=list)
    group_order: List[int] = field(default_factory=list)
    labels: Optional[np.ndarray] = None
    output_dir: str = ""
    ini_path: str = ""
    n_contacts: int = 0
    cluster_summary: Optional[dict] = None


def run_extract_group_volumes(
    mesh,
    labels: np.ndarray,
    summary,  # LithologySummary
    local_origin: tuple = (0.0, 0.0, 0.0),
    z_datum: float = 0.0,
    output_dir: str = "lithoseed_output",
    export_formats: tuple = ("ts", "obj"),
    clean_min_cells: int = 25,
    gravity_xyz: Optional[np.ndarray] = None,
    gravity_obs: Optional[np.ndarray] = None,
    gravity_std: Optional[np.ndarray] = None,
    magnetic_xyz: Optional[np.ndarray] = None,
    magnetic_obs: Optional[np.ndarray] = None,
    magnetic_std: Optional[np.ndarray] = None,
    borehole_constraints: Optional[list] = None,
    igrf_f: float = 55000.0,
    igrf_i: float = -66.0,
    igrf_d: float = 0.0,
    dem_csv: Optional[str] = None,
    has_susceptibility: bool = False,
    magnetic_weight: float = 1.0,
    min_density_range: float = 0.2,
) -> PipelineResult:
    """Extract closed volumes per lithology group via voxel-face boundaries.

    Each group gets one or more closed volumes (disconnected spatial
    components become separate volumes).  Volumes are non-overlapping
    by construction (shared vertices, zero gaps, squared-off edges).
    No decimation is applied — shared boundary vertices remain identical.

    When *clean_min_cells* > 0, small disconnected label fragments
    are reassigned to neighbouring groups before extraction.  Set to
    0 to keep all fragments (useful for debugging).

    Writes a full resolved_config.ini with ``group_meshes`` and
    ``label_grid`` keys for the C++ volumetric loading path.
    """
    os.makedirs(output_dir, exist_ok=True)

    if clean_min_cells > 0:
        labels = cleanup_label_fragments(mesh, labels, min_cells=clean_min_cells)

    group_ids = sorted(set(labels[labels >= 0]))

    volumes = extract_group_volumes(
        mesh, labels, group_ids,
        local_origin=local_origin,
    )

    mesh_dir = os.path.join(output_dir, "meshes")
    if os.path.isdir(mesh_dir):
        import shutil
        shutil.rmtree(mesh_dir)
    os.makedirs(mesh_dir, exist_ok=True)

    # Group volumes by label and merge all disconnected components
    # into a single closed mesh per label.  The C++ side expects exactly
    # one mesh per lithology group (setGroupMesh), so we must merge.
    vols_by_label: dict[int, list] = {}
    for vol in volumes:
        vols_by_label.setdefault(vol.group_above, []).append(vol)

    # Sort components within each label by size (largest first) for
    # deterministic output, then merge.
    merged_volumes = {}
    for label, comps in vols_by_label.items():
        comps.sort(key=lambda c: len(c.faces), reverse=True)
        merged_volumes[label] = merge_volume_components(comps)

    # Write one TS/OBJ per merged group (label).  Filename uses the
    # original pipeline label so we don't conflate labels that map to
    # the same cluster.
    written_labels = []
    for label in sorted(merged_volumes):
        vol = merged_volumes[label]
        name = f"volume_group_{label}"
        if "ts" in export_formats:
            ts_path = os.path.join(mesh_dir, f"{name}.ts")
            write_volume_ts(ts_path, vol, 0)
        if "obj" in export_formats:
            obj_path = os.path.join(mesh_dir, f"{name}.obj")
            write_volume_obj(obj_path, vol, 0)
        written_labels.append(label)

    total_components = sum(len(c) for c in vols_by_label.values())
    print(f"  Extracted {total_components} closed volume(s) across {len(merged_volumes)} labels"
          f" ({len(written_labels)} written as merged meshes)")
    for label in sorted(merged_volumes):
        n_comp = len(vols_by_label[label])
        comp_str = f" ({n_comp} component{'s' if n_comp > 1 else ''})" if n_comp > 1 else ""
        print(f"    Label {label}: {len(merged_volumes[label].vertices)} verts,"
              f" {len(merged_volumes[label].faces)} faces{comp_str}")

    cluster_full = _summary_to_dict(summary)
    n_full = len(cluster_full["density_mean"])

    # Build cluster property dict and group_mesh_paths, preserving the
    # 1:1 label-to-cluster mapping from cluster_intersection (labels are
    # 1-indexed sequential IDs; cluster_full arrays are 0-indexed).
    # All labels are kept — negative SimPEG densities get a small positive
    # floor so they remain valid for the C++ gravity forward model.  The
    # C++ inversion applies its own density bounds from the INI.
    cluster_dict: dict = {
        "group_id": [],
        "density_mean": [],
        "density_std": [],
    }
    if "susc_mean" in cluster_full:
        cluster_dict["susc_mean"] = []
        cluster_dict["susc_std"] = []

    # label_grid → group_id remapping (1-indexed labels → 0-indexed groups)
    label_to_group: dict[int, int] = {}

    filtered_mesh_paths = []
    for label in sorted(merged_volumes):
        cluster_idx = label - 1  # 1-indexed label → 0-indexed cluster array
        if cluster_idx < 0 or cluster_idx >= n_full:
            print(f"  [filter] label {label}: cluster index {cluster_idx} out of"
                  f" range [0, {n_full}), excluded")
            continue
        d_mean = cluster_full["density_mean"][cluster_idx]

        new_gid = len(cluster_dict["group_id"])
        label_to_group[label] = new_gid
        cluster_dict["group_id"].append(new_gid)
        cluster_dict["density_mean"].append(d_mean)
        cluster_dict["density_std"].append(cluster_full["density_std"][cluster_idx])
        if "susc_mean" in cluster_full:
            cluster_dict["susc_mean"].append(cluster_full["susc_mean"][cluster_idx])
            cluster_dict["susc_std"].append(cluster_full["susc_std"][cluster_idx])

        name = f"volume_group_{label}"
        ts_path = os.path.join(mesh_dir, f"{name}.ts")
        filtered_mesh_paths.append(ts_path)

    n_filtered = len(cluster_dict["group_id"])
    if n_filtered < len(merged_volumes):
        print(f"  [filter] kept {n_filtered} / {len(merged_volumes)} labels"
              f" (excluded {len(merged_volumes) - n_filtered})")

    # Convert SimPEG density contrasts to absolute densities.
    # SimPEG works in contrast space (Δρ = ρ − ρ_background), but the
    # C++ litho-inversion expects absolute densities (ρ in g/cm³).
    # Scale contrasts to ensure at least min_density_range separation,
    # then shift by background density.
    BACKGROUND_DENSITY = 2.67
    d_arr = np.array(cluster_dict["density_mean"])
    d_range = d_arr.max() - d_arr.min()
    if d_range > 0 and d_range < min_density_range:
        scale = min_density_range / d_range
        cluster_dict["density_mean"] = list(d_arr * scale)
        cluster_dict["density_std"] = list(
            np.array(cluster_dict["density_std"]) * scale)
        print(f"  [density scaling] range {d_range:.4f} -> {min_density_range:.4f}"
              f" (scale={scale:.1f}x)")
        d_arr = np.array(cluster_dict["density_mean"])
    cluster_dict["density_mean"] = list(d_arr + BACKGROUND_DENSITY)
    print(f"  [density absolute] means: {[f'{v:.3f}' for v in cluster_dict['density_mean']]} g/cm³")

    cluster_path = os.path.join(output_dir, "cluster_properties.csv")
    write_cluster_csv(cluster_path, cluster_dict)

    # Export label grid — remap 1-indexed labels to 0-indexed group IDs
    # so classifyPoint's raw lookup matches the group index in the
    # cluster CSV / LithologyModel groups vector.
    labels_3d = mesh.to_3d_grid(labels, fill=-1).astype(np.int32)
    labels_remapped = np.full_like(labels_3d, -1)
    for old_label, new_gid in label_to_group.items():
        labels_remapped[labels_3d == old_label] = new_gid
    label_grid_path = os.path.join(mesh_dir, "label_grid.bin")
    export_label_grid(label_grid_path, labels_remapped, mesh)

    # Write observed data CSVs
    grav_path = None
    if gravity_xyz is not None and gravity_obs is not None:
        grav_path = os.path.join(output_dir, "observed_gravity.csv")
        write_observed_gravity_csv(grav_path, gravity_xyz, gravity_obs, gravity_std)

    mag_path = None
    if magnetic_xyz is not None and magnetic_obs is not None:
        mag_path = os.path.join(output_dir, "observed_magnetic.csv")
        write_observed_magnetic_csv(mag_path, magnetic_xyz, magnetic_obs, magnetic_std)

    bh_path = None
    if borehole_constraints:
        bh_path = os.path.join(output_dir, "borehole_constraints.csv")
        write_borehole_constraints_csv(bh_path, borehole_constraints)

    # Write INI config with group_meshes (filtered primary components) + label_grid
    ini_path = os.path.join(output_dir, "resolved_config.ini")
    write_ini_config(
        ini_path,
        cluster_csv=cluster_path,
        group_mesh_paths=filtered_mesh_paths,
        observed_gravity_csv=grav_path,
        observed_magnetic_csv=mag_path,
        borehole_csv=bh_path,
        label_grid=label_grid_path,
        local_origin_x=local_origin[0],
        local_origin_y=local_origin[1],
        z_datum=z_datum,
        igrf_f=igrf_f,
        igrf_i=igrf_i,
        igrf_d=igrf_d,
        dem_csv=dem_csv,
        has_susceptibility=has_susceptibility,
        magnetic_weight=magnetic_weight,
    )

    return PipelineResult(
        contacts=volumes,
        group_order=group_ids,
        labels=labels,
        output_dir=output_dir,
        ini_path=ini_path,
        n_contacts=len(volumes),
        cluster_summary=cluster_dict,
    )


def run_extract(
    mesh,
    density: np.ndarray,
    susceptibility: Optional[np.ndarray] = None,
    n_clusters: int = 4,
    cluster_random_state: int = 42,
    target_vertices: int = 500,
    local_origin: tuple = (0.0, 0.0, 0.0),
    z_datum: float = 0.0,
    output_dir: str = "lithoseed_output",
    export_formats: tuple = ("ts", "obj"),
    gravity_xyz: Optional[np.ndarray] = None,
    gravity_obs: Optional[np.ndarray] = None,
    gravity_std: Optional[np.ndarray] = None,
    magnetic_xyz: Optional[np.ndarray] = None,
    magnetic_obs: Optional[np.ndarray] = None,
    magnetic_std: Optional[np.ndarray] = None,
    borehole_constraints: Optional[list] = None,
    igrf_f: float = 55000.0,
    igrf_i: float = -66.0,
    igrf_d: float = 0.0,
    dem_csv: Optional[str] = None,
    min_density_range: float = 0.2,
) -> PipelineResult:
    """Run the extract pipeline: cluster → extract → decimate → export."""
    from cluster_api._cluster import cluster_lithology

    if susceptibility is None:
        susceptibility = np.zeros_like(density)

    labels, summary = cluster_lithology(
        density, susceptibility,
        n_clusters=n_clusters,
        random_state=cluster_random_state,
    )

    d_means = np.array(summary.density_mean)
    d_range = d_means.max() - d_means.min()
    if d_range > 0 and d_range < min_density_range:
        scale = min_density_range / d_range
        summary.density_mean = list(d_means * scale)
        summary.density_std = list(np.array(summary.density_std) * scale)
        print(f"  [density scaling] range {d_range:.4f} -> {min_density_range:.4f}"
              f" g/cc (scale={scale:.1f}x)")

    # Convert SimPEG density contrasts to absolute densities.
    BACKGROUND_DENSITY = 2.67
    summary.density_mean = [d + BACKGROUND_DENSITY for d in summary.density_mean]

    cluster_dict = _summary_to_dict(summary)
    return _extract_contacts_and_export(
        mesh, density, susceptibility, labels, cluster_dict,
        target_vertices, local_origin, z_datum, output_dir,
        export_formats, gravity_xyz, gravity_obs, gravity_std,
        magnetic_xyz, magnetic_obs, magnetic_std,
        borehole_constraints, igrf_f, igrf_i, igrf_d, dem_csv,
    )


def run_extract_from_labels(
    mesh,
    density: np.ndarray,
    susceptibility: np.ndarray,
    labels: np.ndarray,
    summary,  # LithologySummary from cluster_intersection or cluster_lithology
    target_vertices: int = 500,
    local_origin: tuple = (0.0, 0.0, 0.0),
    z_datum: float = 0.0,
    output_dir: str = "lithoseed_output",
    export_formats: tuple = ("ts", "obj"),
    gravity_xyz: Optional[np.ndarray] = None,
    gravity_obs: Optional[np.ndarray] = None,
    gravity_std: Optional[np.ndarray] = None,
    magnetic_xyz: Optional[np.ndarray] = None,
    magnetic_obs: Optional[np.ndarray] = None,
    magnetic_std: Optional[np.ndarray] = None,
    borehole_constraints: Optional[list] = None,
    igrf_f: float = 55000.0,
    igrf_i: float = -66.0,
    igrf_d: float = 0.0,
    dem_csv: Optional[str] = None,
) -> PipelineResult:
    """Extract contact surfaces from pre-computed labels and export.

    Use this when labels come from intersection clustering or any
    external clustering method.  Skips internal GMM clustering.
    """
    cluster_dict = _summary_to_dict(summary)
    return _extract_contacts_and_export(
        mesh, density, susceptibility, labels, cluster_dict,
        target_vertices, local_origin, z_datum, output_dir,
        export_formats, gravity_xyz, gravity_obs, gravity_std,
        magnetic_xyz, magnetic_obs, magnetic_std,
        borehole_constraints, igrf_f, igrf_i, igrf_d, dem_csv,
    )


def _summary_to_dict(summary) -> dict:
    """Convert LithologySummary to a plain dict for cluster_dict."""
    d = {
        "group_id": list(range(summary.n_groups)),
        "density_mean": list(summary.density_mean),
        "density_std": list(summary.density_std),
    }
    if hasattr(summary, "susc_mean") and summary.susc_mean:
        d["susc_mean"] = list(summary.susc_mean)
    if hasattr(summary, "susc_std") and summary.susc_std:
        d["susc_std"] = list(summary.susc_std)
    return d


def _estimate_contact_depth(ga: int, gb: int, existing: list) -> float:
    """Estimate depth for a missing contact pair from nearby contacts."""
    if existing:
        return float(np.median([cs.median_depth for cs in existing]))
    return -200.0


def _make_flat_contact(x0: float, y0: float, wx: float, wy: float,
                        z: float, n: int = 11) -> tuple:
    """Create a flat square contact surface spanning the mesh area."""
    x = np.linspace(x0, x0 + wx, n)
    y = np.linspace(y0, y0 + wy, n)
    xv, yv = np.meshgrid(x, y)
    verts = np.column_stack([xv.ravel(), yv.ravel(), np.full(n * n, z)])
    faces = []
    for iy in range(n - 1):
        for ix in range(n - 1):
            i0 = iy * n + ix
            i1 = i0 + 1
            i2 = (iy + 1) * n + ix
            i3 = i2 + 1
            faces.append([i0, i1, i2])
            faces.append([i1, i3, i2])
    return verts.astype(np.float64), np.array(faces, dtype=np.int32)


def _extract_contacts_and_export(
    mesh,
    density: np.ndarray,
    susceptibility: np.ndarray,
    labels: np.ndarray,
    cluster_dict: dict,
    target_vertices: int,
    local_origin: tuple,
    z_datum: float,
    output_dir: str,
    export_formats: tuple,
    gravity_xyz: Optional[np.ndarray],
    gravity_obs: Optional[np.ndarray],
    gravity_std: Optional[np.ndarray],
    magnetic_xyz: Optional[np.ndarray],
    magnetic_obs: Optional[np.ndarray],
    magnetic_std: Optional[np.ndarray],
    borehole_constraints: Optional[list],
    igrf_f: float,
    igrf_i: float,
    igrf_d: float,
    dem_csv: Optional[str],
) -> PipelineResult:
    """Core extraction + export shared by run_extract and run_extract_from_labels."""
    os.makedirs(output_dir, exist_ok=True)

    group_order = sorted(set(labels[labels >= 0]))

    contacts_raw = extract_contact_surfaces(
        mesh, labels, group_order,
        local_origin=local_origin,
    )

    contacts = []
    for cs in contacts_raw:
        if len(cs.vertices) > target_vertices:
            dec_v, dec_f = decimate_mesh(cs.vertices, cs.faces, target_vertices)
            cs = ContactSurface(
                group_above=cs.group_above,
                group_below=cs.group_below,
                vertices=dec_v,
                faces=dec_f,
                median_depth=cs.median_depth,
            )
        contacts.append(cs)

    # Filter small disconnected noise patches; warn about multi-component surfaces.
    contacts_filtered = []
    for cs in contacts:
        components = separate_connected_components(cs, min_faces=10)
        if len(components) > 1:
            print(f"  [components] contact {cs.group_above}_{cs.group_below}:"
                  f" {len(components)} connected components"
                  f" (sizes: {[len(c.faces) for c in components]}),"
                  f" keeping largest")
        contacts_filtered.append(components[0])
    contacts = contacts_filtered

    # Enforce stratigraphic ordering so surfaces don't cross in 3D.
    contacts = enforce_stratigraphic_ordering(contacts, margin=5.0)

    # Fill in missing contact pairs with flat fallback surfaces.
    # Some group pairs may have no spatial adjacency in the 3D mesh
    # (e.g. intersection clustering can produce non-adjacent groups).
    existing_pairs = {(cs.group_above, cs.group_below) for cs in contacts}
    for i in range(len(group_order) - 1):
        ga = group_order[i]
        gb = group_order[i + 1]
        if (ga, gb) not in existing_pairs:
            est_z = _estimate_contact_depth(ga, gb, contacts)
            verts, faces = _make_flat_contact(mesh.x0, mesh.y0,
                                              mesh.nx * mesh.dx,
                                              mesh.ny * mesh.dy, est_z)
            contacts.append(ContactSurface(
                group_above=ga, group_below=gb,
                vertices=verts, faces=faces,
                median_depth=est_z,
            ))
            print(f"  [fallback] contact {ga}_{gb}: no spatial adjacency,"
                  f" using flat surface at z={est_z:.0f}m")

    # Re-sort contacts by group pair
    contacts.sort(key=lambda c: (c.group_above, c.group_below))

    mesh_dir = os.path.join(output_dir, "meshes")
    os.makedirs(mesh_dir, exist_ok=True)

    mesh_paths = []
    for cs in contacts:
        name = f"contact_{cs.group_above}_{cs.group_below}"
        if "ts" in export_formats:
            ts_path = os.path.join(mesh_dir, f"{name}.ts")
            write_contact_ts(ts_path, cs)
            mesh_paths.append(ts_path)
        if "obj" in export_formats:
            obj_path = os.path.join(mesh_dir, f"{name}.obj")
            write_contact_obj(obj_path, cs)

    cluster_path = os.path.join(output_dir, "cluster_properties.csv")
    write_cluster_csv(cluster_path, cluster_dict)

    grav_path = None
    if gravity_xyz is not None and gravity_obs is not None:
        grav_path = os.path.join(output_dir, "observed_gravity.csv")
        write_observed_gravity_csv(grav_path, gravity_xyz, gravity_obs, gravity_std)

    mag_path = None
    if magnetic_xyz is not None and magnetic_obs is not None:
        mag_path = os.path.join(output_dir, "observed_magnetic.csv")
        write_observed_magnetic_csv(mag_path, magnetic_xyz, magnetic_obs, magnetic_std)

    bh_path = None
    if borehole_constraints:
        bh_path = os.path.join(output_dir, "borehole_constraints.csv")
        write_borehole_constraints_csv(bh_path, borehole_constraints)

    has_susc = susceptibility is not None and np.any(np.abs(susceptibility) > 1e-12)

    ini_path = os.path.join(output_dir, "resolved_config.ini")
    write_ini_config(
        ini_path,
        cluster_csv=cluster_path,
        group_mesh_paths=mesh_paths,
        observed_gravity_csv=grav_path,
        observed_magnetic_csv=mag_path,
        borehole_csv=bh_path,
        local_origin_x=local_origin[0],
        local_origin_y=local_origin[1],
        z_datum=z_datum,
        igrf_f=igrf_f,
        igrf_i=igrf_i,
        igrf_d=igrf_d,
        dem_csv=dem_csv,
        has_susceptibility=has_susc,
    )

    return PipelineResult(
        contacts=contacts,
        group_order=group_order,
        labels=labels,
        output_dir=output_dir,
        ini_path=ini_path,
        n_contacts=len(contacts),
        cluster_summary=cluster_dict,
    )

