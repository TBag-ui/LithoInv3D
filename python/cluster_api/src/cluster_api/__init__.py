"""
Lithology classification and surface extraction from joint inversion results.

Usage
-----
    from api import load_inversion, cluster_lithology, extract_surfaces

    # Load SimPEG joint inversion output
    mesh, density, susc = load_inversion(
        mesh_path="forrestania_joint.msh",
        density_path="forrestania_density.xyz",
        susceptibility_path="forrestania_susceptibility.xyz",
    )

    # Cluster into distinct lithology groups
    labels, summary = cluster_lithology(
        density, susc, n_clusters=4, random_state=42,
    )

    # Extract closed triangulated surfaces
    surfaces = extract_surfaces(
        mesh, labels, summary,
        output_dir="surfaces",
        formats=("obj", "ts"),
    )
"""

from ._io import load_inversion, MeshData
from ._cluster import cluster_lithology, LithologySummary
from ._surface import extract_surfaces, SurfaceResult

__all__ = [
    "load_inversion",
    "MeshData",
    "cluster_lithology",
    "LithologySummary",
    "extract_surfaces",
    "SurfaceResult",
]
