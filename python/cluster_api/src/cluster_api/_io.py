"""
Mesh and model file I/O for SimPEG inversion outputs.

Supports:
  - SimPEG tensor-mesh .msh files with active-cell lists
  - Per-cell .xyz property files (X Y Z value)
  - Per-cell .mod  property files (value only, aligned with mesh active cells)
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional
import os


@dataclass
class MeshData:
    """A 3-D tensor mesh with active-cell indexing."""

    nx: int
    ny: int
    nz: int
    x0: float
    y0: float
    z0: float
    dx: float
    dy: float
    dz: float
    n_active: int
    ix: np.ndarray   # 1-based x-indices of active cells  [n_active]
    iy: np.ndarray   # 1-based y-indices
    iz: np.ndarray   # 1-based z-indices
    x_center: np.ndarray  # cell-centre X coordinates  [n_active]
    y_center: np.ndarray  # cell-centre Y coordinates
    z_center: np.ndarray  # cell-centre Z coordinates

    @property
    def n_total(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def cell_volume(self) -> float:
        return self.dx * self.dy * self.dz

    def to_3d_grid(self, values: np.ndarray, fill: float = np.nan) -> np.ndarray:
        """Map per-active-cell *values* into a dense [nx, ny, nz] grid."""
        grid = np.full((self.nx, self.ny, self.nz), fill, dtype=values.dtype)
        for i in range(self.n_active):
            grid[self.ix[i] - 1, self.iy[i] - 1, self.iz[i] - 1] = values[i]
        return grid

    def node_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (xn, yn, zn) node-coordinate arrays for marching cubes."""
        xn = self.x0 + np.arange(self.nx + 1) * self.dx
        yn = self.y0 + np.arange(self.ny + 1) * self.dy
        zn = self.z0 - np.arange(self.nz + 1) * self.dz
        return xn, yn, zn


# ---------------------------------------------------------------------------
# Public loaders
# ---------------------------------------------------------------------------


def load_inversion(
    mesh_path: str,
    density_path: Optional[str] = None,
    susceptibility_path: Optional[str] = None,
) -> Tuple[MeshData, np.ndarray, np.ndarray]:
    """
    Load a SimPEG joint-inversion result.

    Parameters
    ----------
    mesh_path : str
        Path to the .msh tensor-mesh file.
    density_path : str, optional
        Path to density values.  One of:
          - .xyz file (X Y Z value per line, header line skipped)
          - .mod file (one value per active cell, same order as mesh)
    susceptibility_path : str, optional
        Same format as *density_path*.

    Returns
    -------
    mesh : MeshData
    density : np.ndarray [n_active]
    susceptibility : np.ndarray [n_active]
    """
    mesh = _read_msh(mesh_path)

    density = _read_property(density_path, mesh.n_active)
    susceptibility = _read_property(susceptibility_path, mesh.n_active)

    # Compute cell-centre coordinates from xyz file if available, else from mesh
    if density_path and os.path.splitext(density_path)[1] == '.xyz':
        xyz = np.loadtxt(density_path, skiprows=1)
        mesh.x_center = xyz[:, 0]
        mesh.y_center = xyz[:, 1]
        mesh.z_center = xyz[:, 2]
    else:
        mesh.x_center = mesh.x0 + (mesh.ix - 0.5) * mesh.dx
        mesh.y_center = mesh.y0 + (mesh.iy - 0.5) * mesh.dy
        mesh.z_center = mesh.z0 - (mesh.iz - 0.5) * mesh.dz

    return mesh, density, susceptibility


# ---------------------------------------------------------------------------
# Internal readers
# ---------------------------------------------------------------------------


def _read_msh(path: str) -> MeshData:
    """Read a SimPEG tensor-mesh .msh file with active-cell list."""
    with open(path) as f:
        lines = f.readlines()

    nx, ny, nz = map(int, lines[0].split())
    x0, y0, z0 = map(float, lines[1].split())
    dx, dy, dz = map(float, lines[2].split())
    n_active = int(lines[3])

    ix = np.zeros(n_active, dtype=int)
    iy = np.zeros(n_active, dtype=int)
    iz = np.zeros(n_active, dtype=int)

    for i, line in enumerate(lines[4 : 4 + n_active]):
        parts = line.strip().split()
        ix[i] = int(parts[0])
        iy[i] = int(parts[1])
        iz[i] = int(parts[2])

    return MeshData(
        nx=nx, ny=ny, nz=nz,
        x0=x0, y0=y0, z0=z0,
        dx=dx, dy=dy, dz=dz,
        n_active=n_active,
        ix=ix, iy=iy, iz=iz,
        x_center=np.empty(n_active),
        y_center=np.empty(n_active),
        z_center=np.empty(n_active),
    )


def _read_property(path: Optional[str], n_expected: int) -> np.ndarray:
    """Read a property file, returning a numpy array of length *n_expected*.

    Supports:
      - .xyz  : whitespace-delimited, first line is a header → column 4
      - .mod  : one float per line, no header
    """
    if path is None:
        return np.full(n_expected, np.nan)

    ext = os.path.splitext(path)[1].lower()

    if ext == '.xyz':
        data = np.loadtxt(path, skiprows=1)
        if data.ndim == 2 and data.shape[1] >= 4:
            return data[:, 3].astype(np.float64)
        raise ValueError(f"{path}: expected >=4 columns in .xyz file")

    if ext == '.mod':
        with open(path) as f:
            first = f.readline().split()
        skip = 1 if len(first) == 3 else 0
        data = np.loadtxt(path, skiprows=skip)
        if data.ndim == 1:
            return data.astype(np.float64)
        raise ValueError(f"{path}: expected 1-D .mod file")

    raise ValueError(f"Unsupported property file extension: {ext}")

