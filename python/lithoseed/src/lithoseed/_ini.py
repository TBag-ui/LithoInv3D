"""INI config generator for the C++ litho-constrained inversion."""

import os
from typing import List, Optional


def write_ini_config(
    path: str,
    cluster_csv: str,
    group_mesh_paths: List[str],
    observed_gravity_csv: Optional[str] = None,
    observed_magnetic_csv: Optional[str] = None,
    borehole_csv: Optional[str] = None,
    local_origin_x: float = 0.0,
    local_origin_y: float = 0.0,
    z_datum: float = 0.0,
    # ── inversion section ─────────────────────────────────────────
    vertex_freedom: str = "xyz_free",
    max_iterations: int = 50,
    tolerance: float = 1e-5,
    lambda_reg: float = 1.0,
    omega: float = 10.0,
    lbfgs_history: int = 20,
    fd_step: float = 10.0,
    control_point_stride: int = 5,
    armijo_c1: float = 0.001,
    line_search_max_iter: int = 30,
    solver: str = "lbfgsb",
    enable_eigenvalue_scaling: bool = True,
    enable_reference_model: bool = False,
    lambda_ref: float = 0.1,
    enable_depth_bounds: bool = False,
    depth_bound_margin: float = 100.0,
    enable_property_inversion: bool = True,
    property_inversion_interval: int = 5,
    property_inversion_max_iter: int = 10,
    # ── gravity section ───────────────────────────────────────────
    gravity_uncertainty: float = 0.0,
    density_min: float = -0.5,
    density_max: float = 0.5,
    # ── magnetic section ──────────────────────────────────────────
    magnetic_weight: float = 1.0,
    magnetic_uncertainty: float = 0.0,
    susceptibility_min: float = 0.0,
    susceptibility_max: float = 0.1,
    remanence_min: float = 0.0,
    remanence_max: float = 10.0,
    remanence_mode: str = "effective_susceptibility",
    ## IGRF ────────────────────────────────────────────────────────
    igrf_f: float = 55000.0,
    igrf_i: float = -66.0,
    igrf_d: float = 0.0,
    has_susceptibility: bool = False,
    ## topography ──────────────────────────────────────────────────
    label_grid: Optional[str] = None,
    dem_csv: Optional[str] = None,
    topo_mode: str = "none",
    bouguer_density: float = 2.67,
    padding_rings: int = 0,
    padding_cell_size: float = 100.0,
    invert_halfspace_properties: bool = False,
    ## padding group ───────────────────────────────────────────────
    enable_padding_group: bool = True,
    padding_density_initial: float = 2.68,
    padding_density_lower: float = 1.5,
    padding_density_upper: float = 5.0,
    padding_depth: float = -100000.0,
) -> None:
    """Generate a resolved INI config for csv_invert_fixed.exe.

    All paths are written relative to the INI file's directory.
    """
    ini_dir = os.path.dirname(os.path.abspath(path))

    def rel(p):
        if p is None:
            return ""
        return os.path.relpath(p, ini_dir).replace("\\", "/")

    meshes = ",".join(rel(m) for m in group_mesh_paths)

    sections: List[str] = []

    # ── [data] ──────────────────────────────────────────────────────
    sections.append("[data]")
    sections.append("real_data = true")
    sections.append("group_column = cluster_id")
    sections.append(f"cluster_csv = {rel(cluster_csv)}")
    if observed_gravity_csv:
        sections.append(f"observed_gravity_csv = {rel(observed_gravity_csv)}")
    if observed_magnetic_csv:
        sections.append(f"observed_magnetic_csv = {rel(observed_magnetic_csv)}")
        if not has_susceptibility:
            magnetic_weight = 0.0
    sections.append(f"group_meshes = {meshes}")
    if label_grid:
        sections.append(f"label_grid = {rel(label_grid)}")
    if borehole_csv:
        sections.append(f"borehole_classified_csvs = {rel(borehole_csv)}")
    sections.append(f"local_origin_x = {local_origin_x}")
    sections.append(f"local_origin_y = {local_origin_y}")
    sections.append(f"z_datum = {z_datum}")
    sections.append("fail_on_missing_meshes = true")

    # ── [inversion] ─────────────────────────────────────────────────
    sections.append("")
    sections.append("[inversion]")
    sections.append(f"solver = {solver}")
    sections.append(f"max_iterations = {max_iterations}")
    sections.append(f"tolerance = {tolerance}")
    sections.append(f"lambda = {lambda_reg}")
    sections.append(f"omega = {omega}")
    sections.append(f"lbfgs_history = {lbfgs_history}")
    sections.append(f"fd_step = {fd_step}")
    sections.append(f"control_point_stride = {control_point_stride}")
    sections.append(f"armijo_c1 = {armijo_c1}")
    sections.append(f"line_search_max_iter = {line_search_max_iter}")
    sections.append(f"vertex_freedom = {vertex_freedom}")
    sections.append(f"enable_eigenvalue_scaling = {_bool_cfg(enable_eigenvalue_scaling)}")
    sections.append(f"enable_property_inversion = {_bool_cfg(enable_property_inversion)}")
    sections.append(f"property_inversion_interval = {property_inversion_interval}")
    sections.append(f"property_inversion_max_iter = {property_inversion_max_iter}")
    if enable_reference_model:
        sections.append(f"enable_reference_model = true")
        sections.append(f"lambda_ref = {lambda_ref}")
    if enable_depth_bounds:
        sections.append(f"enable_depth_bounds = true")
        sections.append(f"depth_bound_margin = {depth_bound_margin}")

    # ── [gravity] ───────────────────────────────────────────────────
    if observed_gravity_csv:
        sections.append("")
        sections.append("[gravity]")
        sections.append(f"density_min = {density_min}")
        sections.append(f"density_max = {density_max}")
        sections.append(f"gravity_uncertainty = {gravity_uncertainty}")

    # ── [magnetic] ──────────────────────────────────────────────────
    if observed_magnetic_csv:
        sections.append("")
        sections.append("[magnetic]")
        sections.append(f"magnetic_weight = {magnetic_weight}")
        sections.append(f"magnetic_uncertainty = {magnetic_uncertainty}")
        sections.append(f"field_nT = {igrf_f}")
        sections.append(f"inclination = {igrf_i}")
        sections.append(f"declination = {igrf_d}")
        sections.append(f"susceptibility_min = {susceptibility_min}")
        sections.append(f"susceptibility_max = {susceptibility_max}")
        sections.append(f"remanence_min = {remanence_min}")
        sections.append(f"remanence_max = {remanence_max}")
        sections.append(f"remanence_mode = {remanence_mode}")

    # ── [topography] ────────────────────────────────────────────────
    sections.append("")
    sections.append("[topography]")
    sections.append(f"mode = {topo_mode}")
    if topo_mode == "raw" and dem_csv:
        sections.append(f"dem_file = {rel(dem_csv)}")
    else:
        sections.append("dem_file =")
    sections.append(f"datum_elevation = {z_datum}")
    sections.append(f"bouguer_density = {bouguer_density}")
    sections.append(f"padding_rings = {padding_rings}")
    sections.append(f"padding_cell_size = {padding_cell_size}")
    sections.append(f"invert_halfspace_properties = {_bool_cfg(invert_halfspace_properties)}")

    # ── [padding] ───────────────────────────────────────────────────
    sections.append("")
    sections.append("[padding]")
    sections.append(f"enable_padding_group = {_bool_cfg(enable_padding_group)}")
    sections.append(f"padding_density_initial = {padding_density_initial}")
    sections.append(f"padding_density_lower = {padding_density_lower}")
    sections.append(f"padding_density_upper = {padding_density_upper}")
    sections.append(f"padding_depth = {padding_depth}")

    with open(path, "w") as f:
        f.write("\n".join(sections) + "\n")


def _bool_cfg(v: bool) -> str:
    return "true" if v else "false"

