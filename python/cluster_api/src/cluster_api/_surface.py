"""
Closed triangulated surface extraction from clustered lithology labels.

Extracts an isosurface (marching cubes) for each lithology group, then
exports to OBJ and/or GOCAD TSurf (.ts) formats.
"""

import numpy as np
import os
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from ._io import MeshData


@dataclass
class SurfaceResult:
    """One extracted surface."""

    group: int
    n_cells: int
    vertices: np.ndarray   # (M, 3)
    faces: np.ndarray      # (K, 3)  0-based indices
    obj_path: Optional[str] = None
    ts_path: Optional[str] = None


# ---------------------------------------------------------------------------
# Public
# ---------------------------------------------------------------------------


def extract_surfaces(
    mesh: MeshData,
    labels: np.ndarray,
    summary,
    output_dir: str = "surfaces",
    formats: Sequence[str] = ("obj", "ts"),
    min_cells: int = 8,
) -> List[SurfaceResult]:
    """
    Extract closed triangulated surfaces for each lithology group.

    Parameters
    ----------
    mesh : MeshData
        The tensor mesh (must have ix/iy/iz active-cell indices).
    labels : (N,) int array
        Per-active-cell lithology label (-1 = sentinel, 1..K = group).
    summary : LithologySummary
        Group property statistics (embedded in OBJ/TS headers).
    output_dir : str
        Directory to write surface files.
    formats : sequence of str
        Subset of {"obj", "ts"}.
    min_cells : int
        Skip groups with fewer than this many cells.

    Returns
    -------
    surfaces : list of SurfaceResult
    """
    from skimage import measure
    from scipy.ndimage import binary_closing

    os.makedirs(output_dir, exist_ok=True)

    # Map labels to 3-D grid
    label_grid = mesh.to_3d_grid(labels.astype(np.int16), fill=-1)

    # Node coordinates for marching cubes
    xn, yn, zn = mesh.node_arrays()

    surfaces = []

    for k in range(1, summary.n_groups + 1):
        binary = (label_grid == k).astype(np.float32)
        n_vox = int(binary.sum())

        if n_vox < min_cells:
            continue

        # Morphological close to bridge single-cell gaps
        try:
            binary = binary_closing(binary, structure=np.ones((3, 3, 3)), iterations=1)
        except Exception:
            pass

        # Marching cubes
        spacing = (mesh.dx, mesh.dy, mesh.dz)
        try:
            verts, faces, _, _ = measure.marching_cubes(
                binary, level=0.5, spacing=spacing,
                allow_degenerate=False,
            )
        except Exception:
            # Pad with a layer of zeros to ensure closure
            padded = np.pad(binary, 1, mode='constant', constant_values=0)
            verts, faces, _, _ = measure.marching_cubes(
                padded, level=0.5, spacing=spacing,
                allow_degenerate=False,
            )
            verts[:, 0] += xn[0] - mesh.dx
            verts[:, 1] += yn[0] - mesh.dy
            verts[:, 2] += zn[-1] - mesh.dz
            result = _build_result(k, n_vox, verts, faces, summary, output_dir, formats)
            surfaces.append(result)
            continue

        # Shift to world coordinates
        verts[:, 0] += mesh.x0
        verts[:, 1] += mesh.y0
        verts[:, 2] += mesh.z0 - mesh.nz * mesh.dz

        result = _build_result(k, n_vox, verts, faces, summary, output_dir, formats)
        surfaces.append(result)

    return surfaces


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------


def _build_result(
    group: int,
    n_cells: int,
    vertices: np.ndarray,
    faces: np.ndarray,
    summary,
    output_dir: str,
    formats: Sequence[str],
) -> SurfaceResult:
    result = SurfaceResult(
        group=group, n_cells=n_cells,
        vertices=vertices, faces=faces,
    )

    # Look up group properties
    idx = summary.labels.index(group)
    d_mean = summary.density_mean[idx]
    s_mean = summary.susc_mean[idx]

    base = os.path.join(output_dir, f"litho_group_{group}")

    if "obj" in formats:
        p = base + ".obj"
        _write_obj(p, vertices, faces, group, d_mean, s_mean)
        result.obj_path = p

    if "ts" in formats:
        p = base + ".ts"
        _write_ts(p, vertices, faces, group, d_mean, s_mean)
        result.ts_path = p

    return result


def _write_obj(
    path: str,
    verts: np.ndarray,
    faces: np.ndarray,
    group: int,
    d_mean: float,
    s_mean: float,
) -> None:
    with open(path, 'w') as f:
        f.write(f"# Lithology group {group}\n")
        f.write(f"# density_mean={d_mean:.6f} g/cc  susc_mean={s_mean:.6f} SI\n")
        f.write(f"# vertices={len(verts)} faces={len(faces)}\n")
        for v in verts:
            f.write(f"v {v[0]:.3f} {v[1]:.3f} {v[2]:.3f}\n")
        for tri in faces:
            f.write(f"f {tri[0]+1} {tri[1]+1} {tri[2]+1}\n")


def _write_ts(
    path: str,
    verts: np.ndarray,
    faces: np.ndarray,
    group: int,
    d_mean: float,
    s_mean: float,
) -> None:
    with open(path, 'w') as f:
        f.write("GOCAD TSurf 1\n")
        f.write("HEADER {\n")
        f.write(f"name: litho_group_{group}\n")
        f.write(f"*density_mean: {d_mean:.6f} g/cc\n")
        f.write(f"*susc_mean: {s_mean:.6f} SI\n")
        f.write("}\n")
        f.write("TRIANGLES\n")
        for i, v in enumerate(verts):
            f.write(f"VRTX {i+1} {v[0]:.3f} {v[1]:.3f} {v[2]:.3f}\n")
        for tri in faces:
            f.write(f"TRGL {tri[0]+1} {tri[1]+1} {tri[2]+1}\n")
        f.write("END\n")
