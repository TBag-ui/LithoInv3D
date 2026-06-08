#!/usr/bin/env python3
"""
Forrestania Post-Inversion Script — launches forrestania_invert.exe and
collects outputs into the presentation directory.

Usage:
    python forrestania_post_inversion.py combination    # GMM intersection
    python forrestania_post_inversion.py user            # lithology gates
    python forrestania_post_inversion.py both            # both in sequence
"""

import os
import sys
import shutil
import subprocess
import time
from pathlib import Path

PROJ_ROOT = Path(__file__).resolve().parent.parent
PRESENTATION = PROJ_ROOT / "apps" / "forrestania" / "presentation"
EXE = PROJ_ROOT / "build" / "release" / "forrestania_invert.exe"

# ── Per-cluster config overrides ────────────────────────────────────────────

CLUSTER_CONFIGS = {
    "combination": {
        "ini_name": "config_combination.ini",
        "volumes_rel": "../volumes/combination_clusters",
        "output_subdir": "./inprogress_combination",
        "enable_property_inversion": True,
        "property_inversion_interval": 100,
        "property_inversion_max_iter": 10,
    },
    "user": {
        "ini_name": "config_user.ini",
        "volumes_rel": "../volumes/user_clusters",
        "output_subdir": "./inprogress_user",
        "enable_property_inversion": False,
        "property_inversion_interval": 100,
        "property_inversion_max_iter": 10,
    },
}


def build_ini(cluster_key):
    """Build an INI config for a specific cluster type."""
    cfg = CLUSTER_CONFIGS[cluster_key]
    volumes_rel = cfg["volumes_rel"]

    # Discover volume group mesh files
    meshes_abs = str(PRESENTATION / "configs" / volumes_rel / "meshes")
    group_meshes = []
    if os.path.isdir(meshes_abs):
        for fname in sorted(os.listdir(meshes_abs)):
            if fname.startswith("volume_group_") and fname.endswith(".ts"):
                group_meshes.append(f"{volumes_rel}/meshes/{fname}")
    group_meshes_str = ",".join(group_meshes) if group_meshes else ""

    prop_inv = "true" if cfg["enable_property_inversion"] else "false"

    ini = f"""[bounds]
depth_bound_margin = 500.000000
enable_depth_bounds = true

[data]
real_data = true
group_column = cluster_id
cluster_csv = {volumes_rel}/cluster_properties.csv
group_meshes = {group_meshes_str}
label_grid = {volumes_rel}/meshes/label_grid.bin
local_origin_x = 747019.390000
local_origin_y = 6417259.750000
observed_gravity_csv = ../volumes/observed_gravity.csv
observed_magnetic_csv = ../volumes/observed_magnetic.csv
z_datum = 381.232137
fail_on_missing_meshes = false

[gravity]
density_max = 0.500000
density_min = -0.500000
gravity_uncertainty = 0.000000

[inversion]
solver = lbfgsb
max_iterations = 250
tolerance = 1e-10
lambda = 1.0
omega = 10.0
lbfgs_history = 20
fd_step = 10.0
control_point_stride = 0
armijo_c1 = 0.001
line_search_max_iter = 30
vertex_freedom = xyz_free
mesh_boundary_mode = free
enable_eigenvalue_scaling = true
enable_geometry_inversion = true
enable_property_inversion = {prop_inv}
property_inversion_interval = {cfg['property_inversion_interval']}
property_inversion_max_iter = {cfg['property_inversion_max_iter']}
property_damping = 0.010000
gncg_cg_max_iter = 50
gncg_cg_tolerance = 0.000001
disable_line_search = true

[magnetic]
declination = -0.100000
field_nT = 58874.000000
inclination = -66.200000
magnetic_uncertainty = 0.000000
magnetic_weight = 1.000000
remanence_component_max = 10.000000
remanence_component_min = -10.000000
remanence_max = 10.000000
remanence_min = 0.000000
remanence_mode = effective_susceptibility
susceptibility_max = 0.100000
susceptibility_min = 0.000000

[output]
iteration_export_dir = {cfg['output_subdir']}

[padding]
enable_padding_group = false
padding_conductivity_initial = 0.000100
padding_conductivity_lower = 0.000010
padding_conductivity_upper = 100.000000
padding_density_initial = 2.680000
padding_density_lower = 1.500000
padding_density_upper = 5.000000
padding_depth = -100000.000000
padding_susceptibility_initial = 0.000000
padding_susceptibility_lower = 0.000000
padding_susceptibility_upper = 1.000000

[regularization]
enable_reference_model = false
gradient_smoothing_weight = 0.300000
lambda_ref = 0.100000

[topography]
mode = none
dem_file =
datum_elevation = 381.232137
bouguer_density = 2.67
padding_rings = 0
padding_cell_size = 100.0
invert_halfspace_properties = false
projected_beyond_survey_area = independent_flat

[true_model]
dip_angle = 0.000000
dip_direction = 90.000000
mode = voiseys
"""
    return ini


# Coordinate transform constants (from pre-pipeline step 1)
LOCAL_ORIGIN_X = 747019.39
LOCAL_ORIGIN_Y = 6417259.75
Z_DATUM = 381.232137


def parse_convergence(stdout_text, output_dir):
    """Parse iteration lines from C++ stdout, write convergence CSV."""
    import re, csv

    lines = [l.strip() for l in stdout_text.splitlines()
             if l.strip().startswith("Iter ")]
    if not lines:
        print("  WARNING: No convergence lines found in stdout")
        return

    csv_path = os.path.join(output_dir, "convergence.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iter", "rms", "objective",
                     "dw_g_x", "dw_g_y", "mag_misfit", "dw_m_x", "dw_m_y"])
        pat = re.compile(
            r"Iter\s+(\d+):\s+RMS=([\d.\-nan(ind)]+)\s+"
            r"obj=([\d.e+\-]+)\s+DW_g=\(([\d.]+),([\d.]+)\)\s+"
            r"mag=([\d.e+\-]+)\s+DW_m=\(([\d.]+),([\d.]+)\)")
        for line in lines:
            m = pat.match(line)
            if m:
                w.writerow(list(m.groups()))
    print(f"  Convergence: {len(lines)} iterations -> {csv_path}")


def translate_surface_to_utm(ts_path, dst_path):
    """Translate a model-space .ts surface to UTM data coordinates.

    Model space: local origin at (LOCAL_ORIGIN_X, LOCAL_ORIGIN_Y, Z_DATUM)
    Data space:  UTM MGA Zone 50
    """
    with open(ts_path) as f:
        content = f.read()

    lines_out = []
    for line in content.splitlines():
        parts = line.split()
        # GOCAD TSurf format: VRTX <id> <x> <y> <z>
        if line.startswith("VRTX") and len(parts) >= 5:
            try:
                vid = parts[1]
                x = float(parts[2]) + LOCAL_ORIGIN_X
                y = float(parts[3]) + LOCAL_ORIGIN_Y
                z = float(parts[4]) + Z_DATUM
                lines_out.append(f"VRTX {vid} {x:.6f} {y:.6f} {z:.6f}")
            except ValueError:
                lines_out.append(line)
        else:
            lines_out.append(line)

    os.makedirs(os.path.dirname(dst_path) or ".", exist_ok=True)
    with open(dst_path, "w") as f:
        f.write("\n".join(lines_out) + "\n")


def run_inversion(cluster_key):
    """Run forrestania_invert.exe and collect all outputs."""
    cfg = CLUSTER_CONFIGS[cluster_key]

    if not EXE.is_file():
        print(f"ERROR: forrestania_invert.exe not found at {EXE}")
        print("  Build it first:  .\\build\\build.bat")
        return 1

    # ── Write config ────────────────────────────────────────────────────
    configs_dir = PRESENTATION / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)
    ini_path = configs_dir / cfg["ini_name"]

    print(f"Writing config: {ini_path}")
    ini_path.write_text(build_ini(cluster_key))

    # ── Output directory ─────────────────────────────────────────────────
    output_dir = PRESENTATION / "inversion_output" / cluster_key
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    # Iteration export dir (relative to config location)
    iter_dir = configs_dir / cfg["output_subdir"].lstrip("./")
    iter_dir.mkdir(parents=True, exist_ok=True)

    # ── Launch ───────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"RUNNING: {cluster_key}")
    print(f"  Config:     {ini_path}")
    print(f"  Output:     {output_dir}")
    print(f"  Iterations: {iter_dir}")
    print(f"  Property inversion: {cfg['enable_property_inversion']}")
    print(f"{'='*70}\n")

    t0 = time.time()

    result = subprocess.run(
        [str(EXE), str(ini_path)],
        cwd=str(configs_dir),
        capture_output=True, text=True,
    )
    # Stream captured output to terminal
    print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    elapsed = time.time() - t0

    if result.returncode != 0:
        print(f"\nERROR: forrestania_invert.exe exited with code {result.returncode}")
        return result.returncode

    print(f"\nCompleted in {elapsed/60:.1f} minutes")

    # ── 1. Parse convergence ─────────────────────────────────────────────
    parse_convergence(result.stdout, str(output_dir))

    # ── 2. Copy iteration exports ────────────────────────────────────────
    if iter_dir.exists():
        dst_iter = output_dir / "iterations"
        if dst_iter.exists():
            shutil.rmtree(dst_iter)
        shutil.copytree(iter_dir, dst_iter)
        n_dirs = len([d for d in dst_iter.iterdir() if d.is_dir()])
        print(f"  Iteration exports: {n_dirs} directories -> {dst_iter}")

    # ── 3. Translate final surfaces to UTM ───────────────────────────────
    # Final iter dir is the last one by name
    iter_dirs = sorted([d for d in dst_iter.iterdir() if d.is_dir()],
                        key=lambda x: x.name)
    if iter_dirs:
        final_iter = iter_dirs[-1]
        utm_dir = output_dir / "final_utm"
        utm_dir.mkdir(exist_ok=True)
        n_translated = 0
        for ts_file in final_iter.glob("*.ts"):
            dst = utm_dir / ts_file.name
            translate_surface_to_utm(str(ts_file), str(dst))
            n_translated += 1
        print(f"  UTM translations: {n_translated} surfaces -> {utm_dir}")

        # Also copy model-space final
        model_dir = output_dir / "final_model_space"
        if model_dir.exists():
            shutil.rmtree(model_dir)
        shutil.copytree(str(final_iter), str(model_dir))
        print(f"  Model-space final: -> {model_dir}")

    # ── 4. Copy the config used ──────────────────────────────────────────
    shutil.copy2(str(ini_path), str(output_dir / cfg["ini_name"]))

    print(f"\nOutput: {output_dir}")
    return 0


def main():
    if len(sys.argv) < 2:
        print("Usage: python forrestania_post_inversion.py [combination|user|both]")
        return 1

    mode = sys.argv[1]

    if mode == "both":
        for key in ["combination", "user"]:
            rc = run_inversion(key)
            if rc != 0:
                print(f"\nABORTED: {key} failed with code {rc}")
                return rc
    elif mode in CLUSTER_CONFIGS:
        rc = run_inversion(mode)
        if rc != 0:
            return rc
    else:
        print(f"Unknown mode: {mode}. Use: combination, user, or both")
        return 1

    print("\nDONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
