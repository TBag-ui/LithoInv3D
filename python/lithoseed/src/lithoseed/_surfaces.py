"""3D contact surface extraction + mesh decimation for LithoSeed."""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional

from skimage.measure import marching_cubes
from scipy import ndimage


@dataclass
class ContactSurface:
    """A triangulated contact surface between two adjacent lithology groups."""
    group_above: int
    group_below: int
    vertices: np.ndarray   # (N, 3) world coordinates
    faces: np.ndarray      # (M, 3) 0-based vertex indices
    median_depth: float


def extract_contact_surfaces(
    mesh,
    labels: np.ndarray,
    group_order: List[int],
    local_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    min_faces: int = 10,
) -> List[ContactSurface]:
    """Extract 3D contact surfaces between adjacent lithology groups.

    For each consecutive pair in group_order, runs marching cubes on a
    binary volume and filters faces to keep only the true contact between
    the two groups.

    Parameters
    ----------
    mesh : MeshData
        Voxel mesh with nx/ny/nz, cell sizes, and origin.
    labels : np.ndarray
        Per-active-cell group labels (length n_active).
    group_order : list of int
        Ordered group IDs from top (shallowest) to bottom (deepest).
    local_origin : tuple
        (x, y, z) offset subtracted from world coordinates.
    min_faces : int
        Contacts with fewer faces than this are discarded.

    Returns
    -------
    list of ContactSurface
    """
    labels_3d = mesh.to_3d_grid(labels, fill=-1).astype(int)
    contacts = []

    for i in range(len(group_order) - 1):
        g_above = group_order[i]
        g_below = group_order[i + 1]

        binary = (labels_3d == g_below).astype(np.float64)

        if binary.sum() < 2:
            continue

        try:
            verts, faces, normals, _ = marching_cubes(binary, level=0.5)
        except (ValueError, RuntimeError):
            continue

        if len(faces) < min_faces:
            continue

        keep = _neighbor_filter(verts, faces, normals, labels_3d, g_above)
        faces = faces[keep]

        if len(faces) < min_faces:
            continue

        used = np.unique(faces)
        remap = np.full(len(verts), -1, dtype=int)
        remap[used] = np.arange(len(used))
        verts = verts[used]
        faces = remap[faces]
        if normals is not None and len(normals) > 0:
            normals = normals[used]

        world_verts = _voxel_to_world(verts, mesh, local_origin)

        median_z = np.median(world_verts[:, 2])

        contacts.append(ContactSurface(
            group_above=g_above,
            group_below=g_below,
            vertices=world_verts,
            faces=faces,
            median_depth=-median_z,
        ))

    return contacts


def _neighbor_filter(
    verts: np.ndarray,
    faces: np.ndarray,
    normals: np.ndarray,
    labels_3d: np.ndarray,
    group_above: int,
) -> np.ndarray:
    """Keep faces whose non-group_below neighbor voxel is group_above.

    For each face, samples the 8-connected corner voxels of each vertex
    (floor/ceil in all three axes). If any sampled voxel has label ==
    group_above, the face is kept. This is robust to faces sitting exactly
    on voxel boundaries.
    """
    nx, ny, nz = labels_3d.shape
    keep = np.zeros(len(faces), dtype=bool)

    for v_idx in range(3):
        coords = verts[faces[:, v_idx]]  # (M, 3)
        for di in [0, 1]:
            for dj in [0, 1]:
                for dk in [0, 1]:
                    if di == 0:
                        vi = np.floor(coords[:, 0]).astype(int)
                    else:
                        vi = np.minimum(np.floor(coords[:, 0]).astype(int) + 1, nx - 1)
                    if dj == 0:
                        vj = np.floor(coords[:, 1]).astype(int)
                    else:
                        vj = np.minimum(np.floor(coords[:, 1]).astype(int) + 1, ny - 1)
                    if dk == 0:
                        vk = np.floor(coords[:, 2]).astype(int)
                    else:
                        vk = np.minimum(np.floor(coords[:, 2]).astype(int) + 1, nz - 1)

                    vi = np.clip(vi, 0, nx - 1)
                    vj = np.clip(vj, 0, ny - 1)
                    vk = np.clip(vk, 0, nz - 1)

                    keep |= (labels_3d[vi, vj, vk] == group_above)

    return keep


def _voxel_to_world(
    verts: np.ndarray,
    mesh,
    local_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    """Convert marching-cubes voxel coordinates to world coordinates.

    Marching cubes returns vertices in index space [0, nx] × [0, ny] × [0, nz].
    MeshData convention: x0/y0 is the grid origin, z0 is the top surface,
    z increases downward in index space.

    Output is world coordinates (no origin subtraction). The C++ inversion
    applies the local-origin transform when loading meshes.
    """
    world = np.empty_like(verts)
    world[:, 0] = mesh.x0 + verts[:, 0] * mesh.dx
    world[:, 1] = mesh.y0 + verts[:, 1] * mesh.dy
    world[:, 2] = mesh.z0 - verts[:, 2] * mesh.dz
    return world


def enforce_stratigraphic_ordering(
    surfaces: List[ContactSurface],
    margin: float = 1.0,
) -> List[ContactSurface]:
    """Ensure contact surfaces do not cross in 3D.

    For each consecutive pair (upper, lower), clamps upper vertices that
    dip below the lower surface and lower vertices that rise above the
    upper surface. Uses nearest-neighbour lookup in the XY plane.

    Parameters
    ----------
    surfaces : list of ContactSurface
        Ordered from shallowest (top) to deepest (bottom).
    margin : float
        Minimum vertical separation enforced between surfaces (metres).

    Returns
    -------
    list of ContactSurface (modified in place).
    """
    if len(surfaces) < 2:
        return surfaces

    from scipy.spatial import cKDTree

    for i in range(len(surfaces) - 1):
        upper = surfaces[i]
        lower = surfaces[i + 1]

        lower_xy = lower.vertices[:, :2]
        tree_lower = cKDTree(lower_xy)

        n_clamped = 0
        for j in range(len(upper.vertices)):
            ux, uy, uz = upper.vertices[j]
            _, idx = tree_lower.query([ux, uy])
            lower_z = lower.vertices[idx, 2]
            if uz < lower_z + margin:
                upper.vertices[j, 2] = lower_z + margin
                n_clamped += 1

        upper_xy = upper.vertices[:, :2]
        tree_upper = cKDTree(upper_xy)

        for j in range(len(lower.vertices)):
            lx, ly, lz = lower.vertices[j]
            _, idx = tree_upper.query([lx, ly])
            upper_z = upper.vertices[idx, 2]
            if lz > upper_z - margin:
                lower.vertices[j, 2] = upper_z - margin
                n_clamped += 1

        if n_clamped > 0:
            import sys
            print(f"  [ordering] contact {upper.group_above}_{upper.group_below}"
                  f" / {lower.group_above}_{lower.group_below}:"
                  f" clamped {n_clamped} vertices",
                  file=sys.stderr)

    return surfaces


def separate_connected_components(
    surface: ContactSurface,
    min_faces: int = 5,
) -> List[ContactSurface]:
    """Separate a contact surface into its connected components.

    Parameters
    ----------
    surface : ContactSurface
    min_faces : int
        Components with fewer faces than this are discarded.

    Returns
    -------
    list of ContactSurface, one per connected component that meets the
    minimum face count. If the surface has only one component it is
    returned as a single-element list.
    """
    n_verts = len(surface.vertices)
    n_faces = len(surface.faces)

    if n_faces == 0:
        return []

    # Build face adjacency graph (two faces are connected if they share an edge)
    # Map each undirected edge to the face(s) that contain it
    edge_to_faces = {}
    for fi, (v0, v1, v2) in enumerate(surface.faces):
        for edge in [(v0, v1), (v1, v2), (v2, v0)]:
            key = (min(edge), max(edge))
            edge_to_faces.setdefault(key, []).append(fi)

    # Union-Find over faces
    parent = list(range(n_faces))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for edge_faces in edge_to_faces.values():
        for fi in edge_faces[1:]:
            union(edge_faces[0], fi)

    # Group faces by component
    comp_faces = {}
    for fi in range(n_faces):
        root = find(fi)
        comp_faces.setdefault(root, []).append(fi)

    # Build separate surfaces for components with enough faces
    components = []
    for root, fi_list in comp_faces.items():
        if len(fi_list) < min_faces:
            continue
        faces_sub = surface.faces[fi_list]
        used = np.unique(faces_sub)
        remap = np.full(n_verts, -1, dtype=int)
        remap[used] = np.arange(len(used))
        verts_sub = surface.vertices[used].copy()
        faces_remapped = remap[faces_sub]

        median_z = np.median(verts_sub[:, 2])
        components.append(ContactSurface(
            group_above=surface.group_above,
            group_below=surface.group_below,
            vertices=verts_sub,
            faces=faces_remapped,
            median_depth=-median_z,
        ))

    # Sort by face count descending so largest component is first.
    components.sort(key=lambda c: len(c.faces), reverse=True)
    return components if components else [surface]


def cleanup_label_fragments(
    mesh,
    labels: np.ndarray,
    min_cells: int = 25,
    max_passes: int = 5,
) -> np.ndarray:
    """Reassign cells in small disconnected components to neighbouring groups.

    For each lithology group, finds connected components in 3D.  Components
    with fewer than *min_cells* cells are treated as noise fragments: their
    cells are reassigned to the most common label among 6-connected neighbours.

    Parameters
    ----------
    mesh : MeshData
    labels : (N_active,) int array
    min_cells : int
        Minimum cells for a component to survive.
    max_passes : int
        Maximum iterations of cleanup (fragments can cascade).

    Returns
    -------
    labels_clean : (N_active,) int array
    """
    labels_3d = mesh.to_3d_grid(labels, fill=-1).astype(int)
    nx, ny, nz = labels_3d.shape

    total_changed = 0
    for _pass in range(max_passes):
        pass_changed = 0
        present = set(np.unique(labels_3d))
        present.discard(-1)
        if not present:
            break

        for gid in sorted(present):
            binary = (labels_3d == gid).astype(np.int32)
            labeled, n_comp = ndimage.label(binary)
            if n_comp <= 1:
                continue

            for comp_id in range(1, n_comp + 1):
                comp_mask = labeled == comp_id
                comp_count = int(comp_mask.sum())
                if comp_count >= min_cells:
                    continue

                comp_inds = np.argwhere(comp_mask)
                for ci, cj, ck in comp_inds:
                    nb_labels = []
                    for di, dj, dk in [(-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1)]:
                        ni, nj, nk = ci + di, cj + dj, ck + dk
                        if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz:
                            nb_lab = labels_3d[ni, nj, nk]
                            if nb_lab >= 0 and nb_lab != gid:
                                nb_labels.append(nb_lab)
                    if nb_labels:
                        uniq, counts = np.unique(nb_labels, return_counts=True)
                        new_label = uniq[np.argmax(counts)]
                        labels_3d[ci, cj, ck] = new_label
                        pass_changed += 1

        total_changed += pass_changed
        if pass_changed == 0:
            break

    if total_changed > 0:
        import sys
        print(f"  [cleanup] reassigned {total_changed} cells from small fragments"
              f" (min_cells={min_cells}, passes={_pass + 1})", file=sys.stderr)

    # Map back to active-cell labels
    labels_clean = labels.copy()
    for i in range(len(labels)):
        vi = mesh.ix[i] - 1
        vj = mesh.iy[i] - 1
        vk = mesh.iz[i] - 1
        if 0 <= vi < nx and 0 <= vj < ny and 0 <= vk < nz:
            labels_clean[i] = labels_3d[vi, vj, vk]

    return labels_clean


def extract_group_volumes(
    mesh,
    labels: np.ndarray,
    group_ids: List[int],
    local_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> List[ContactSurface]:
    """Extract non-overlapping closed volumes via voxel-face boundary extraction.

    Iterates every face in the 3D label grid.  When adjacent voxels have
    different labels, the shared quad face is assigned to both label
    volumes (with opposite winding).  Boundary faces close the volume at
    mesh edges.

    Every connected component is kept — no faces are discarded.  Neighbouring
    volumes share identical vertices at their contacts (true zero-gap, no
    decimation that would pull shared vertices apart).

    Returns one ContactSurface per connected component of each group.
    """
    labels_3d = mesh.to_3d_grid(labels, fill=-1).astype(int)
    nx, ny, nz = labels_3d.shape
    dx, dy, dz = mesh.dx, mesh.dy, mesh.dz
    x0, y0, z0 = mesh.x0, mesh.y0, mesh.z0

    # Per-group: list of triangle vertex triples (world coords)
    tris_by_group: dict[int, list] = {gid: [] for gid in group_ids}

    def _world_vertex(vi: float, vj: float, vk: float) -> tuple:
        """Voxel coords → world coords (same convention as _voxel_to_world)."""
        wx = x0 + vi * dx
        wy = y0 + vj * dy
        wz = z0 - vk * dz
        return (wx, wy, wz)

    def _emit_quad(v0, v1, v2, v3, gid_owner, gid_other):
        """Add a quad (2 tris) to gid_ownerʼs volume.

        v0..v3 are tuples (wx, wy, wz) going CCW around the quad
        when viewed from gid_owner toward gid_other.
        """
        tris_by_group.setdefault(gid_owner, []).append((v0, v2, v1))
        tris_by_group.setdefault(gid_owner, []).append((v0, v3, v2))

    # ── Interior faces (shared between two voxels) ─────────────────
    # X-faces (normal ±x)
    for i in range(1, nx):
        for j in range(ny):
            for k in range(nz):
                la = labels_3d[i - 1, j, k]
                lb = labels_3d[i, j, k]
                if la == lb:
                    continue
                vx = float(i)
                v0 = _world_vertex(vx, float(j), float(k))
                v1 = _world_vertex(vx, float(j + 1), float(k))
                v2 = _world_vertex(vx, float(j + 1), float(k + 1))
                v3 = _world_vertex(vx, float(j), float(k + 1))
                # Quad CCW viewed from +x side: (v0, v1, v2, v3)
                if la >= 0:
                    _emit_quad(v0, v1, v2, v3, la, lb)  # outward +x for la
                if lb >= 0:
                    _emit_quad(v0, v3, v2, v1, lb, la)  # outward -x for lb

    # Y-faces (normal ±y)
    for i in range(nx):
        for j in range(1, ny):
            for k in range(nz):
                la = labels_3d[i, j - 1, k]
                lb = labels_3d[i, j, k]
                if la == lb:
                    continue
                vy = float(j)
                v0 = _world_vertex(float(i), vy, float(k))
                v1 = _world_vertex(float(i + 1), vy, float(k))
                v2 = _world_vertex(float(i + 1), vy, float(k + 1))
                v3 = _world_vertex(float(i), vy, float(k + 1))
                if la >= 0:
                    _emit_quad(v0, v1, v2, v3, la, lb)
                if lb >= 0:
                    _emit_quad(v0, v3, v2, v1, lb, la)

    # Z-faces (normal ±z)
    for i in range(nx):
        for j in range(ny):
            for k in range(1, nz):
                la = labels_3d[i, j, k - 1]  # upper voxel
                lb = labels_3d[i, j, k]       # lower voxel
                if la == lb:
                    continue
                vz = float(k)
                v0 = _world_vertex(float(i), float(j), vz)
                v1 = _world_vertex(float(i + 1), float(j), vz)
                v2 = _world_vertex(float(i + 1), float(j + 1), vz)
                v3 = _world_vertex(float(i), float(j + 1), vz)
                if la >= 0:
                    _emit_quad(v0, v1, v2, v3, la, lb)  # outward +z (down) for upper
                if lb >= 0:
                    _emit_quad(v0, v3, v2, v1, lb, la)  # outward -z (up) for lower

    # ── Boundary faces (label vs outside mesh) ────────────────────
    # X = 0
    for j in range(ny):
        for k in range(nz):
            lb = labels_3d[0, j, k]
            if lb < 0:
                continue
            vx = 0.0
            v0 = _world_vertex(vx, float(j), float(k))
            v1 = _world_vertex(vx, float(j + 1), float(k))
            v2 = _world_vertex(vx, float(j + 1), float(k + 1))
            v3 = _world_vertex(vx, float(j), float(k + 1))
            _emit_quad(v0, v3, v2, v1, lb, -1)  # outward -x

    # X = nx
    for j in range(ny):
        for k in range(nz):
            la = labels_3d[nx - 1, j, k]
            if la < 0:
                continue
            vx = float(nx)
            v0 = _world_vertex(vx, float(j), float(k))
            v1 = _world_vertex(vx, float(j + 1), float(k))
            v2 = _world_vertex(vx, float(j + 1), float(k + 1))
            v3 = _world_vertex(vx, float(j), float(k + 1))
            _emit_quad(v0, v1, v2, v3, la, -1)  # outward +x

    # Y = 0
    for i in range(nx):
        for k in range(nz):
            lb = labels_3d[i, 0, k]
            if lb < 0:
                continue
            vy = 0.0
            v0 = _world_vertex(float(i), vy, float(k))
            v1 = _world_vertex(float(i + 1), vy, float(k))
            v2 = _world_vertex(float(i + 1), vy, float(k + 1))
            v3 = _world_vertex(float(i), vy, float(k + 1))
            _emit_quad(v0, v3, v2, v1, lb, -1)

    # Y = ny
    for i in range(nx):
        for k in range(nz):
            la = labels_3d[i, ny - 1, k]
            if la < 0:
                continue
            vy = float(ny)
            v0 = _world_vertex(float(i), vy, float(k))
            v1 = _world_vertex(float(i + 1), vy, float(k))
            v2 = _world_vertex(float(i + 1), vy, float(k + 1))
            v3 = _world_vertex(float(i), vy, float(k + 1))
            _emit_quad(v0, v1, v2, v3, la, -1)

    # Z = 0 (top surface — flat top)
    for i in range(nx):
        for j in range(ny):
            lb = labels_3d[i, j, 0]
            if lb < 0:
                continue
            vz = 0.0
            v0 = _world_vertex(float(i), float(j), vz)
            v1 = _world_vertex(float(i + 1), float(j), vz)
            v2 = _world_vertex(float(i + 1), float(j + 1), vz)
            v3 = _world_vertex(float(i), float(j + 1), vz)
            _emit_quad(v0, v3, v2, v1, lb, -1)  # outward -z (up)

    # Z = nz (bottom)
    for i in range(nx):
        for j in range(ny):
            la = labels_3d[i, j, nz - 1]
            if la < 0:
                continue
            vz = float(nz)
            v0 = _world_vertex(float(i), float(j), vz)
            v1 = _world_vertex(float(i + 1), float(j), vz)
            v2 = _world_vertex(float(i + 1), float(j + 1), vz)
            v3 = _world_vertex(float(i), float(j + 1), vz)
            _emit_quad(v0, v1, v2, v3, la, -1)  # outward +z (down)

    # ── Deduplicate vertices, build ContactSurface per group ─────
    volumes: List[ContactSurface] = []
    for gid in sorted(tris_by_group):
        tris = tris_by_group[gid]
        if not tris:
            continue

        # Deduplicate vertices (world-coord tuples → index)
        vert_dict: dict[tuple, int] = {}
        vert_list: list[tuple] = []
        faces: list[tuple] = []

        for tri in tris:
            fi = []
            for v in tri:
                idx = vert_dict.get(v)
                if idx is None:
                    idx = len(vert_list)
                    vert_dict[v] = idx
                    vert_list.append(v)
                fi.append(idx)
            faces.append(tuple(fi))

        vertices = np.array(vert_list, dtype=np.float64)
        faces_arr = np.array(faces, dtype=np.int32)

        median_wz = float(np.median(vertices[:, 2]))

        cs = ContactSurface(
            group_above=gid,
            group_below=-1,
            vertices=vertices,
            faces=faces_arr,
            median_depth=-median_wz,
        )

        # Split into connected components — keep every component
        # (cleanup_label_fragments already removed noise at the label level)
        components = separate_connected_components(cs, min_faces=1)
        for comp in components:
            volumes.append(comp)

    return volumes


def merge_volume_components(components: List[ContactSurface]) -> ContactSurface:
    """Merge multiple closed volume components into a single ContactSurface.

    All components must belong to the same group (same group_above).
    Vertices and faces are concatenated with face index offsetting.
    Returns a single ContactSurface, or the sole component unchanged if
    there is only one.
    """
    if not components:
        raise ValueError("merge_volume_components: empty components list")
    if len(components) == 1:
        return components[0]

    gid = components[0].group_above
    all_verts = []
    all_faces = []
    vertex_offset = 0

    for comp in components:
        if comp.group_above != gid:
            raise ValueError(
                f"merge_volume_components: all components must share group_above "
                f"(expected {gid}, got {comp.group_above})"
            )
        all_verts.append(comp.vertices)
        all_faces.append(comp.faces + vertex_offset)
        vertex_offset += len(comp.vertices)

    merged_verts = np.vstack(all_verts).astype(np.float64)
    merged_faces = np.vstack(all_faces).astype(np.int32)

    median_wz = float(np.median(merged_verts[:, 2]))

    return ContactSurface(
        group_above=gid,
        group_below=-1,
        vertices=merged_verts,
        faces=merged_faces,
        median_depth=-median_wz,
    )


def export_label_grid(path: str, labels_3d: np.ndarray, mesh) -> None:
    """Write a 3D label grid as binary for the C++ label grid loader.

    Format:
      ASCII header: "nx ny nz x0 y0 z0 dx dy dz\\n"
      Then flat int32 array in k-major order (k*ny*nx + j*nx + i).
    """
    nx, ny, nz = labels_3d.shape
    with open(path, "wb") as f:
        header = f"{nx} {ny} {nz} {mesh.x0} {mesh.y0} {mesh.z0} {mesh.dx} {mesh.dy} {mesh.dz}\n"
        f.write(header.encode("ascii"))
        flat = np.asfortranarray(labels_3d).ravel(order='F')
        f.write(flat.astype(np.int32).tobytes())


def decimate_mesh(
    vertices: np.ndarray,
    faces: np.ndarray,
    target_vertices: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Vertex-clustering mesh decimation.

    Divides the bounding box into a grid of bins, merges vertices within
    each bin to their centroid, and remaps faces. Degenerate faces
    (collapsed edges) are removed.

    Parameters
    ----------
    vertices : (N, 3) array
    faces : (M, 3) array of int
    target_vertices : int
        Approximate number of vertices after decimation.

    Returns
    -------
    (decimated_vertices, decimated_faces)
    """
    n_verts = len(vertices)
    if n_verts <= target_vertices:
        return vertices.copy(), faces.copy()

    vmin = vertices.min(axis=0)
    vmax = vertices.max(axis=0)
    extent = vmax - vmin
    extent[extent < 1e-12] = 1.0

    n_per_axis = max(2, int(np.ceil(target_vertices ** (1.0 / 3.0))))
    bin_size = extent / n_per_axis

    bin_idx = np.floor((vertices - vmin) / bin_size).astype(int)
    bin_idx = np.clip(bin_idx, 0, n_per_axis - 1)

    flat_bin = bin_idx[:, 0] + bin_idx[:, 1] * n_per_axis + bin_idx[:, 2] * n_per_axis * n_per_axis

    unique_bins, inverse = np.unique(flat_bin, return_inverse=True)
    n_new = len(unique_bins)

    new_verts = np.zeros((n_new, 3))
    counts = np.zeros(n_new)
    for i in range(n_verts):
        new_verts[inverse[i]] += vertices[i]
        counts[inverse[i]] += 1
    new_verts /= counts[:, None]

    new_faces = inverse[faces]

    non_degenerate = (
        (new_faces[:, 0] != new_faces[:, 1]) &
        (new_faces[:, 1] != new_faces[:, 2]) &
        (new_faces[:, 0] != new_faces[:, 2])
    )
    new_faces = new_faces[non_degenerate]

    used = np.unique(new_faces)
    remap = np.full(n_new, -1, dtype=int)
    remap[used] = np.arange(len(used))
    new_verts = new_verts[used]
    new_faces = remap[new_faces]

    return new_verts, new_faces
