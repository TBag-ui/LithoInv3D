#!/usr/bin/env python3
"""
Generate density-only and magnetic-only cluster closed volumes for diagnosis.

Loads existing SimPEG inversion results from the Forrestania e2e pipeline,
runs separate 1D GMM clusterings (density-only, magnetic-only), extracts
contact surfaces, and runs the C++ binary to produce starting-model closed
volumes in separate directories so the user can inspect why intersection
cluster volumes overlap in 3D space.

Output:
  inversion_output/density_clusters/   — closed volumes for density-only clusters
  inversion_output/magnetic_clusters/  — closed volumes for magnetic-only clusters
"""
import os
import sys
import subprocess
import numpy as np

PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ_ROOT, "cluster_api", "src"))
sys.path.insert(0, os.path.join(PROJ_ROOT, "lithoseed", "src"))

FORRESTANIA_APP = os.path.join(PROJ_ROOT, "..", "apps", "forrestania")
OUTPUT_DIR = os.path.join(FORRESTANIA_APP, "inversion_output")

# C++ binary from the active project (copy-17), which has the
# closed-volume resampling fixes from exporters.cpp.
_THIS_PROJ = os.path.abspath(os.path.join(PROJ_ROOT, "..", "..",
    "Litho_constrained_inversion_pro - Copy (17)", "build", "release",
    "csv_invert_fixed.exe"))
_ALT_PROJ = os.path.join(PROJ_ROOT, "..", "build", "release", "csv_invert_fixed.exe")
CPP_BINARY = _THIS_PROJ if os.path.isfile(_THIS_PROJ) else _ALT_PROJ

N_CLUSTERS = 4
TARGET_VERTICES = 300
RNG = np.random.default_rng(42)


def load_existing_results():
    """Load the SimPEG mesh, density, and susceptibility from the e2e output."""
    from discretize import TensorMesh
    from cluster_api._io import MeshData
    from run_forrestania_e2e import simpeg_to_meshdata

    msh_path = os.path.join(OUTPUT_DIR, "simpeg_mesh.msh")
    dens_path = os.path.join(OUTPUT_DIR, "simpeg_density.mod")
    susc_path = os.path.join(OUTPUT_DIR, "simpeg_susceptibility.mod")

    # Read the UBC-format tensor mesh using discretize
    simpeg_mesh = TensorMesh.read_UBC(msh_path)

    # Read density and susceptibility (full grid, inactive = 0)
    dens_full = np.loadtxt(dens_path, skiprows=1).ravel()
    susc_full = np.loadtxt(susc_path, skiprows=1).ravel()

    # Reconstruct active mask: cells with non-zero density (SimPEG stores
    # recovered model only in active cells; inactive = exactly 0.0).
    active = np.abs(dens_full) > 1e-12
    n_act = int(active.sum())
    rec_dens = dens_full[active]
    rec_susc = susc_full[active]

    # Build MeshData using the same conversion as the e2e pipeline
    mesh_data = simpeg_to_meshdata(simpeg_mesh, active, rec_dens, n_act)

    print(f"Loaded mesh: {mesh_data.nx}x{mesh_data.ny}x{mesh_data.nz}, {mesh_data.n_active} active cells")
    print(f"  Density range: [{rec_dens.min():.4f}, {rec_dens.max():.4f}]")
    print(f"  Susceptibility range: [{rec_susc.min():.6f}, {rec_susc.max():.6f}]")

    return mesh_data, rec_dens, rec_susc


def extract_clusters(mesh, density, susceptibility, cluster_name, n_density, n_susc,
                     local_origin, z_datum, output_subdir, is_density_only):
    """Run intersection clustering + extract contacts + write INI for one variant."""
    from cluster_api._cluster import cluster_intersection
    from lithoseed._pipeline import run_extract_from_labels

    print(f"\n{'=' * 60}")
    print(f"  {cluster_name}: n_density={n_density}, n_susc={n_susc}")
    print(f"{'=' * 60}")

    labels, summary = cluster_intersection(
        density, susceptibility,
        n_density=n_density, n_susc=n_susc,
        random_state=42,
    )
    print(f"  Groups: {summary.n_groups} (of {n_density * n_susc} possible)")
    for i in range(summary.n_groups):
        print(f"    Group {summary.labels[i]}: "
              f"d={summary.density_mean[i]:.4f}+/-{summary.density_std[i]:.4f}, "
              f"s={summary.susc_mean[i]:.6f}+/-{summary.susc_std[i]:.6f}, "
              f"n={summary.counts[i]}")

    # Use the real observed data CSVs from the main e2e output
    grav_csv = os.path.join(OUTPUT_DIR, "observed_gravity.csv")
    mag_csv = os.path.join(OUTPUT_DIR, "observed_magnetic.csv")
    if os.path.isfile(grav_csv) and os.path.isfile(mag_csv):
        g_data = np.loadtxt(grav_csv, delimiter=",", skiprows=1)
        m_data = np.loadtxt(mag_csv, delimiter=",", skiprows=1)
        g_xyz = g_data[:, :3]; g_obs = g_data[:, 3]
        m_xyz = m_data[:, :3]; m_obs = m_data[:, 3]
    else:
        g_xyz = np.array([[0.0, 0.0, 0.0]]); g_obs = np.array([0.0])
        m_xyz = np.array([[0.0, 0.0, 0.0]]); m_obs = np.array([0.0])

    result = run_extract_from_labels(
        mesh, density, susceptibility, labels, summary,
        target_vertices=TARGET_VERTICES,
        local_origin=local_origin,
        z_datum=z_datum,
        output_dir=output_subdir,
        export_formats=("ts", "obj"),
        gravity_xyz=g_xyz,
        gravity_obs=g_obs,
        magnetic_xyz=m_xyz,
        magnetic_obs=m_obs,
    )

    print(f"  Extracted {result.n_contacts} contact surface(s)")
    for cs in result.contacts:
        print(f"    Contact {cs.group_above}-{cs.group_below}: "
              f"{len(cs.vertices)} verts, depth={cs.median_depth:.0f}m")

    # Patch INI: set enable_property_inversion=false for diagnostic runs,
    # and use magnetic_weight=0 for density-only
    _patch_ini(os.path.join(output_subdir, "resolved_config.ini"),
               is_density_only=is_density_only)

    return result


def _patch_ini(ini_path, is_density_only=False):
    """Patch a generated INI for diagnostic closed-volume runs."""
    with open(ini_path, "r") as f:
        content = f.read()

    lines = content.splitlines()
    new_lines = []
    for line in lines:
        if line.startswith("enable_property_inversion"):
            new_lines.append("enable_property_inversion = false")
        elif line.startswith("max_iterations"):
            new_lines.append("max_iterations = 1")
        elif line.startswith("magnetic_weight") and is_density_only:
            new_lines.append("magnetic_weight = 0.0")
        else:
            new_lines.append(line)

    with open(ini_path, "w") as f:
        f.write("\n".join(new_lines) + "\n")


def run_cpp_inversion(output_subdir, label, max_iter=1):
    """Run the C++ binary to produce closed volumes."""
    ini_path = os.path.join(output_subdir, "resolved_config.ini")
    if not os.path.isfile(ini_path):
        print(f"  ERROR: INI not found: {ini_path}")
        return False

    print(f"\n  Running C++ inversion for {label} (max_iter={max_iter})...")
    print(f"  Binary: {CPP_BINARY}")
    print(f"  Config: {ini_path}")

    try:
        result = subprocess.run(
            [CPP_BINARY, ini_path, f"--max-iter={max_iter}"],
            capture_output=True, text=True, timeout=300,
            cwd=output_subdir,
        )
        if result.returncode != 0:
            print(f"  C++ process exited with code {result.returncode}")
            # Print last 40 lines of stderr for diagnosis
            stderr_lines = result.stderr.strip().splitlines()
            for line in stderr_lines[-40:]:
                print(f"    [stderr] {line}")
            return False

        # Print last 20 lines of stdout
        stdout_lines = result.stdout.strip().splitlines()
        for line in stdout_lines[-20:]:
            print(f"    {line}")
        return True
    except subprocess.TimeoutExpired:
        print(f"  C++ process timed out after 300s")
        return False
    except FileNotFoundError:
        print(f"  ERROR: C++ binary not found at {CPP_BINARY}")
        return False


def main():
    print("=" * 60)
    print("DIAGNOSTIC: Density-Only & Magnetic-Only Cluster Volumes")
    print("=" * 60)

    mesh, density, susceptibility = load_existing_results()

    # Determine local origin and z_datum from the original e2e INI
    ini_path = os.path.join(OUTPUT_DIR, "resolved_config.ini")
    local_origin = (747019.39, 6417259.75, 0.0)  # fallback
    z_datum = 381.232137  # fallback
    if os.path.isfile(ini_path):
        import configparser
        cfg = configparser.ConfigParser()
        cfg.read(ini_path)
        if cfg.has_option("data", "local_origin_x"):
            local_origin = (float(cfg.get("data", "local_origin_x")),
                            float(cfg.get("data", "local_origin_y")),
                            0.0)
        if cfg.has_option("data", "z_datum"):
            z_datum = float(cfg.get("data", "z_datum"))
    print(f"  local_origin=({local_origin[0]:.2f}, {local_origin[1]:.2f}), z_datum={z_datum:.2f}")

    # ── Density-only clustering ──────────────────────────────────────
    dens_dir = os.path.join(OUTPUT_DIR, "density_clusters")
    os.makedirs(dens_dir, exist_ok=True)

    extract_clusters(mesh, density, susceptibility,
                     cluster_name="DENSITY-ONLY",
                     n_density=N_CLUSTERS, n_susc=1,
                     local_origin=local_origin, z_datum=z_datum,
                     output_subdir=dens_dir,
                     is_density_only=True)

    # ── Magnetic-only clustering ─────────────────────────────────────
    mag_dir = os.path.join(OUTPUT_DIR, "magnetic_clusters")
    os.makedirs(mag_dir, exist_ok=True)

    extract_clusters(mesh, density, susceptibility,
                     cluster_name="MAGNETIC-ONLY",
                     n_density=1, n_susc=N_CLUSTERS,
                     local_origin=local_origin, z_datum=z_datum,
                     output_subdir=mag_dir,
                     is_density_only=False)

    # ── Run C++ to produce closed volumes ────────────────────────────
    print(f"\n{'=' * 60}")
    print("RUNNING C++ TO PRODUCE CLOSED VOLUMES")
    print(f"{'=' * 60}")

    dens_ok = run_cpp_inversion(dens_dir, "density clusters", max_iter=1)
    mag_ok = run_cpp_inversion(mag_dir, "magnetic clusters", max_iter=1)

    print(f"\n{'=' * 60}")
    if dens_ok and mag_ok:
        print("DONE. Closed volumes produced in:")
        print(f"  Density:  {dens_dir}")
        print(f"  Magnetic: {mag_dir}")
    else:
        print("PARTIAL FAILURE:")
        print(f"  Density:  {'OK' if dens_ok else 'FAILED'} -> {dens_dir}")
        print(f"  Magnetic: {'OK' if mag_ok else 'FAILED'} -> {mag_dir}")
    print(f"{'=' * 60}")

    return 0 if (dens_ok and mag_ok) else 1


if __name__ == "__main__":
    exit(main())
