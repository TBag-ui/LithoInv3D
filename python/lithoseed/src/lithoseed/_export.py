"""Export writers for LithoSeed: TS, OBJ, CSV formats."""

import numpy as np
from typing import Optional, Dict, List, Any


def write_contact_ts(path: str, surface) -> None:
    """Write a contact surface as GOCAD TSurf (.ts) format.

    GOCAD TSurf is the standard interchange format for geological
    visualization software (Geoscience ANALYST, Leapfrog, etc.).
    """
    with open(path, "w") as f:
        f.write("GOCAD TSurf 1\n")
        f.write("HEADER {\n")
        f.write(f"name:contact_{surface.group_above}_{surface.group_below}\n")
        f.write("}\n")
        f.write("TFACE\n")
        for i, v in enumerate(surface.vertices):
            f.write(f"VRTX {i + 1} {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in surface.faces:
            f.write(f"TRGL {face[0] + 1} {face[1] + 1} {face[2] + 1}\n")
        f.write("END\n")


def write_contact_obj(path: str, surface) -> None:
    """Write a contact surface as Wavefront OBJ format."""
    with open(path, "w") as f:
        f.write(f"# Contact surface: group {surface.group_above} / {surface.group_below}\n")
        f.write(f"# Vertices: {len(surface.vertices)}, Faces: {len(surface.faces)}\n")
        for v in surface.vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in surface.faces:
            f.write(f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}\n")


def write_volume_ts(path: str, surface, comp_idx: int = 0) -> None:
    """Write a closed group volume as GOCAD TSurf (.ts) format."""
    with open(path, "w") as f:
        f.write("GOCAD TSurf 1\n")
        f.write("HEADER {\n")
        f.write(f"name:volume_group_{surface.group_above}_comp_{comp_idx}\n")
        f.write("}\n")
        f.write("TFACE\n")
        for i, v in enumerate(surface.vertices):
            f.write(f"VRTX {i + 1} {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in surface.faces:
            f.write(f"TRGL {face[0] + 1} {face[1] + 1} {face[2] + 1}\n")
        f.write("END\n")


def write_volume_obj(path: str, surface, comp_idx: int = 0) -> None:
    """Write a closed group volume as Wavefront OBJ format."""
    with open(path, "w") as f:
        f.write(f"# Volume: group {surface.group_above}, component {comp_idx}\n")
        f.write(f"# Vertices: {len(surface.vertices)}, Faces: {len(surface.faces)}\n")
        for v in surface.vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in surface.faces:
            f.write(f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}\n")


def write_cluster_csv(path: str, summary: Dict[str, list]) -> None:
    """Write cluster property summary as CSV in C++ loadClusterProperties format.

    Writes the 11-column format expected by the C++ cluster_loader:
    cluster_id, working_name, sample_count, density_median_gcc,
    density_p10, density_p90, susceptibility_median_SI,
    susceptibility_p10, susceptibility_p90,
    has_measured_density, has_measured_susceptibility
    """
    group_ids = summary.get("group_id", [])
    density_mean = summary.get("density_mean", [0.0] * len(group_ids))
    density_std = summary.get("density_std", [0.0] * len(group_ids))
    susc_mean = summary.get("susc_mean", [0.0] * len(group_ids))
    susc_std = summary.get("susc_std", [0.0] * len(group_ids))
    has_susc = any(abs(s) > 1e-12 for s in susc_mean)

    header = ("cluster_id,working_name,sample_count,"
              "density_median_gcc,density_p10,density_p90,"
              "susceptibility_median_SI,susceptibility_p10,susceptibility_p90,"
              "has_measured_density,has_measured_susceptibility")

    with open(path, "w") as f:
        f.write(header + "\n")
        for i, gid in enumerate(group_ids):
            d_mean = density_mean[i]
            d_std = density_std[i]
            s_mean = susc_mean[i] if i < len(susc_mean) else 0.0
            s_std = susc_std[i] if i < len(susc_std) else 0.0

            d_p10 = d_mean - d_std
            d_p90 = d_mean + d_std
            s_p10 = max(0.0, s_mean - s_std) if has_susc else 0.0
            s_p90 = s_mean + s_std if has_susc else 0.0

            row = (f"{gid},group_{gid},0,"
                   f"{d_mean},{d_p10},{d_p90},"
                   f"{s_mean},{s_p10},{s_p90},"
                   f"yes,{'yes' if has_susc else 'no'}")
            f.write(row + "\n")


def write_observed_gravity_csv(
    path: str,
    xyz: np.ndarray,
    g_obs: np.ndarray,
    g_std: Optional[np.ndarray] = None,
) -> None:
    """Write observed gravity data as CSV (compatible with C++ loader)."""
    has_std = g_std is not None
    with open(path, "w") as f:
        header = "x,y,z,g_obs"
        if has_std:
            header += ",g_std"
        f.write(header + "\n")
        for i in range(len(g_obs)):
            line = f"{xyz[i, 0]:.6f},{xyz[i, 1]:.6f},{xyz[i, 2]:.6f},{g_obs[i]:.6f}"
            if has_std:
                line += f",{g_std[i]:.6f}"
            f.write(line + "\n")


def write_observed_magnetic_csv(
    path: str,
    xyz: np.ndarray,
    t_obs: np.ndarray,
    t_std: Optional[np.ndarray] = None,
) -> None:
    """Write observed magnetic data as CSV (compatible with C++ loader)."""
    has_std = t_std is not None
    with open(path, "w") as f:
        header = "x,y,z,t_obs"
        if has_std:
            header += ",t_std"
        f.write(header + "\n")
        for i in range(len(t_obs)):
            line = f"{xyz[i, 0]:.6f},{xyz[i, 1]:.6f},{xyz[i, 2]:.6f},{t_obs[i]:.6f}"
            if has_std:
                line += f",{t_std[i]:.6f}"
            f.write(line + "\n")


def write_borehole_constraints_csv(
    path: str,
    constraints: List[Dict[str, Any]],
) -> None:
    """Write borehole constraints CSV (compatible with C++ CSVConstraintLoader).

    Each constraint is a dict with keys: x, y, z_top, z_bottom, litho_group_id.
    z_top and z_bottom are positive-down depths.
    """
    with open(path, "w") as f:
        f.write("x,y,z_top,z_bottom,litho_group_id\n")
        for c in constraints:
            f.write(f"{c['x']:.6f},{c['y']:.6f},{c['z_top']:.6f},"
                    f"{c['z_bottom']:.6f},{c['litho_group_id']}\n")
