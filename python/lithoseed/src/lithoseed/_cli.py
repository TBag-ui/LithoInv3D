"""CLI subcommands for LithoSeed."""

import sys
import os


def main():
    """Entry point for 'python -m LithoSeed'."""
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("LithoSeed - SimPEG starting-model pipeline for litho inversion")
        print()
        print("Subcommands:")
        print("  extract   Extract contact surfaces from existing inversion")
        print("  pipeline  Full pipeline: invert -> cluster -> extract -> export")
        print()
        print("Usage:")
        print("  python -m LithoSeed extract --mesh model.msh --density dens.mod ...")
        print("  python -m LithoSeed pipeline --mode joint --gravity-csv data.csv ...")
        return

    cmd = sys.argv[1]
    sys.argv = [sys.argv[0]] + sys.argv[2:]

    if cmd == "extract":
        _run_extract()
    elif cmd == "pipeline":
        _run_pipeline()
    else:
        print(f"Unknown subcommand: {cmd}")
        sys.exit(1)


def _run_extract():
    """CLI for extract-only mode (no SimPEG needed)."""
    from ._config import build_extract_parser
    from ._pipeline import run_extract

    p = build_extract_parser()
    args = p.parse_args()

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from cluster_api._io import load_inversion

    mesh, density = load_inversion(args.mesh, args.density)
    susc = None
    if args.susceptibility:
        _, susc = load_inversion(args.mesh, args.susceptibility)

    result = run_extract(
        mesh, density, susc,
        n_clusters=args.n_clusters,
        cluster_random_state=args.cluster_seed,
        target_vertices=args.target_vertices,
        local_origin=(args.origin_x, args.origin_y, 0.0),
        z_datum=args.z_datum,
        output_dir=args.output_dir,
    )

    print(f"Extracted {result.n_contacts} contact surface(s)")
    print(f"Output: {result.output_dir}")
    print(f"INI:    {result.ini_path}")


def _run_pipeline():
    """CLI for full pipeline (requires SimPEG)."""
    from ._config import build_parser, SimPEGConfig

    p = build_parser()
    args = p.parse_args()

    cfg = SimPEGConfig.from_argparse(args)
    if hasattr(args, 'no_cross_grad') and args.no_cross_grad:
        cfg.enable_cross_gradient = False
    if hasattr(args, 'no_plot') and args.no_plot:
        cfg.plot = False

    cfg.validate()

    print(f"LithoSeed pipeline: mode={cfg.mode}, mesh={cfg.mesh_type}, "
          f"cell_size={cfg.base_cell_size}m")
    print("Full pipeline (SimPEG inversion) is not yet implemented.")
    print("Use 'extract' subcommand with existing inversion results.")
