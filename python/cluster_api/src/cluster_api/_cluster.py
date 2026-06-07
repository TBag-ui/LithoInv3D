"""
Clustering of inversion physical-property models into distinct lithology groups.

Uses Gaussian Mixture Models on (density, susceptibility) pairs.  A sentinel
value (default -9999) marks cells outside the region of interest; these are
excluded from clustering and assigned label -1.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

try:
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler
    _HAS_SKLEARN = True
except ImportError:
    _HAS_SKLEARN = False


@dataclass
class LithologySummary:
    """Per-group physical-property statistics."""

    labels: List[int]
    counts: List[int]
    density_mean: List[float]
    density_std: List[float]
    susc_mean: List[float]
    susc_std: List[float]

    def to_csv(self, path: str) -> None:
        import csv
        with open(path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(["LithoGroup", "NCells", "DensityMean_gcc",
                         "DensityStd_gcc", "SuscMean_SI", "SuscStd_SI"])
            for i in range(len(self.labels)):
                w.writerow([
                    self.labels[i], self.counts[i],
                    f"{self.density_mean[i]:.6f}", f"{self.density_std[i]:.6f}",
                    f"{self.susc_mean[i]:.6f}", f"{self.susc_std[i]:.6f}",
                ])

    @property
    def n_groups(self) -> int:
        return len(self.labels)


# ---------------------------------------------------------------------------
# Public
# ---------------------------------------------------------------------------


def cluster_lithology(
    density: np.ndarray,
    susceptibility: np.ndarray,
    n_clusters: int = 4,
    sentinel: float = -9999.0,
    random_state: int = 42,
) -> Tuple[np.ndarray, LithologySummary]:
    """
    Cluster (density, susceptibility) into *n_clusters* lithology groups.

    Cells whose density is below *sentinel* + 1 are excluded and assigned
    label -1.  The remaining cells are standardised and clustered with a
    Gaussian Mixture Model.  Clusters are ordered by increasing mean density.

    Parameters
    ----------
    density : (N,) array
    susceptibility : (N,) array
    n_clusters : int
    sentinel : float
    random_state : int

    Returns
    -------
    labels : (N,) int array       -1=sentinel, 1..n_clusters=litho group
    summary : LithologySummary
    """
    if not _HAS_SKLEARN:
        raise ImportError("scikit-learn is required for clustering")

    sentinel_mask = density < (sentinel + 1.0)
    valid = ~sentinel_mask

    d_val = density[valid].copy()
    s_val = susceptibility[valid].copy()

    # Standardise
    features = np.column_stack([d_val, s_val])
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)

    # Cluster
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type='full',
        init_params='random_from_data',
        random_state=random_state,
        n_init=10,
        max_iter=500,
    )
    labels_valid = gmm.fit_predict(features_scaled)

    # Order clusters by mean density (felsic → mafic)
    means = []
    for k in range(n_clusters):
        mask = labels_valid == k
        means.append((k, d_val[mask].mean()))
    means.sort(key=lambda x: x[1])

    remap = {old: new + 1 for new, (old, _) in enumerate(means)}
    labels_valid_ordered = np.array([remap[l] for l in labels_valid])

    # Full label array
    labels = np.full(len(density), -1, dtype=int)
    labels[valid] = labels_valid_ordered

    # Build summary
    summary = LithologySummary(
        labels=[], counts=[], density_mean=[], density_std=[],
        susc_mean=[], susc_std=[],
    )
    for k in range(1, n_clusters + 1):
        m = labels == k
        summary.labels.append(k)
        summary.counts.append(int(m.sum()))
        summary.density_mean.append(float(density[m].mean()))
        summary.density_std.append(float(density[m].std()))
        summary.susc_mean.append(float(susceptibility[m].mean()))
        summary.susc_std.append(float(susceptibility[m].std()))

    return labels, summary


def cluster_intersection(
    density: np.ndarray,
    susceptibility: np.ndarray,
    n_density: int = 4,
    n_susc: int = 3,
    sentinel: float = -9999.0,
    random_state: int = 42,
) -> Tuple[np.ndarray, LithologySummary]:
    """
    Intersection clustering: 1D GMM on density × 1D GMM on susceptibility.

    1. Run 1D GMM on density alone → labels_dens (0..n_density-1)
    2. Run 1D GMM on susceptibility alone → labels_susc (0..n_susc-1)
    3. Combined label = labels_dens * n_susc + labels_susc
    4. Map non-empty combined labels to sequential group IDs (1..N)
    5. Order groups by (density_mean, susc_mean) within each density band

    This preserves the "layered earth" from density while subdividing by
    magnetic character — cells in the same density band but different
    susceptibility clusters get different group IDs.

    Returns (labels, summary) where labels uses sequential group IDs
    starting at 1 (sentinel / excluded cells = -1).
    """
    if not _HAS_SKLEARN:
        raise ImportError("scikit-learn is required for clustering")

    sentinel_mask = density < (sentinel + 1.0)
    valid = ~sentinel_mask

    d_val = density[valid].copy()
    s_val = susceptibility[valid].copy()
    n_valid = len(d_val)

    # --- 1D GMM on density ---
    scaler_d = StandardScaler()
    d_scaled = scaler_d.fit_transform(d_val.reshape(-1, 1))
    gmm_d = GaussianMixture(
        n_components=n_density, covariance_type="full",
        init_params="random_from_data", random_state=random_state,
        n_init=10, max_iter=500,
    )
    labels_dens = gmm_d.fit_predict(d_scaled)

    # Order density clusters by mean density
    d_means_raw = []
    for k in range(n_density):
        mask = labels_dens == k
        d_means_raw.append((k, float(d_val[mask].mean()) if mask.any() else 0.0))
    d_means_raw.sort(key=lambda x: x[1])
    remap_d = {old: new for new, (old, _) in enumerate(d_means_raw)}
    labels_dens = np.array([remap_d[l] for l in labels_dens])

    # --- 1D GMM on susceptibility ---
    scaler_s = StandardScaler()
    s_scaled = scaler_s.fit_transform(s_val.reshape(-1, 1))
    gmm_s = GaussianMixture(
        n_components=n_susc, covariance_type="full",
        init_params="random_from_data", random_state=random_state,
        n_init=10, max_iter=500,
    )
    labels_susc = gmm_s.fit_predict(s_scaled)

    # Order susceptibility clusters by mean susceptibility
    s_means_raw = []
    for k in range(n_susc):
        mask = labels_susc == k
        s_means_raw.append((k, float(s_val[mask].mean()) if mask.any() else 0.0))
    s_means_raw.sort(key=lambda x: x[1])
    remap_s = {old: new for new, (old, _) in enumerate(s_means_raw)}
    labels_susc = np.array([remap_s[l] for l in labels_susc])

    # --- Intersection: combined label = labels_dens * n_susc + labels_susc ---
    combined_raw = labels_dens * n_susc + labels_susc

    # Map non-empty combinations to sequential group IDs
    # Order: by density/depth layer first, then by susc within each layer.
    # combined_raw already encodes this: lower layer labels have lower
    # combined_raw values, and within a layer lower susc means lower raw.
    present_raw = sorted(set(np.unique(combined_raw)))
    raw_to_seq = {raw: i + 1 for i, raw in enumerate(present_raw)}

    labels_valid = np.array([raw_to_seq[c] for c in combined_raw])

    # Full label array
    labels = np.full(len(density), -1, dtype=int)
    labels[valid] = labels_valid
    n_groups = len(present_raw)

    # Build per-group summary (compute mean density/susc per sequential group)
    summary = LithologySummary(
        labels=[], counts=[], density_mean=[], density_std=[],
        susc_mean=[], susc_std=[],
    )
    for k in range(1, n_groups + 1):
        m = labels == k
        summary.labels.append(k)
        summary.counts.append(int(m.sum()))
        summary.density_mean.append(float(density[m].mean()))
        summary.density_std.append(float(density[m].std()))
        summary.susc_mean.append(float(susceptibility[m].mean()))
        summary.susc_std.append(float(susceptibility[m].std()))

    return labels, summary
