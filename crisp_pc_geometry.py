#!/usr/bin/env python3
# ##############################################################
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
# ##############################################################
#
# Script : scripts/crisp_pc_geometry.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Sub-step 6j — PC geometry fingerprint for COMPASS-AI.
#
#   Every cohort produces a cloud of points in 20-dimensional PC
#   space. The shape, size, density, position, and inter-group
#   relationships of that cloud are a geometric fingerprint of the
#   cohort's genetic architecture. This script captures that
#   fingerprint in full and writes it to the Step 6 JSON and .xaux
#   file for COMPASS-AI training.
#
#   We never discard metrics on the grounds that they seem
#   unnecessary. COMPASS-AI decides what is useful during training
#   across thousands of cohorts — our job is to capture everything
#   faithfully.
#
# Metrics computed:
#
#   COHORT-LEVEL (full PC cloud):
#     pc_mean              — mean of each PC (20 values)
#     pc_sd                — SD of each PC (20 values)
#     covariance_matrix    — full 20×20 covariance matrix
#     correlation_matrix   — 20×20 PC correlation matrix
#     eigenvalues          — eigenvalues of the covariance matrix
#                            (principal directions of the cloud)
#     eigenvectors         — eigenvectors of the covariance matrix
#     cloud_volume         — log determinant of covariance matrix
#                            (log used for numerical stability)
#     elongation_ratio     — largest / smallest eigenvalue
#                            (sphere = 1.0, sausage >> 1.0)
#     sphericity           — geometric mean / arithmetic mean of
#                            eigenvalues (1.0 = perfect sphere)
#     mean_nn_distance     — mean nearest-neighbour distance in PC
#                            space (density proxy — small = compact)
#     effective_dimensionality — number of PCs needed to explain
#                            90% of within-cohort variance
#     kurtosis_per_pc      — excess kurtosis of each PC distribution
#                            (detects heavy tails / multi-modality)
#     skewness_per_pc      — skewness of each PC distribution
#
#   PER-GROUP (one block per GMM cluster):
#     group_centroid       — mean PC vector per cluster
#     group_covariance     — covariance matrix per cluster
#     group_volume         — log determinant per cluster
#     group_elongation     — eigenvalue ratio per cluster
#     group_n              — sample count per cluster
#     group_density        — mean intra-cluster NN distance
#
#   INTER-GROUP (pairwise, all group pairs):
#     bhattacharyya_coeff  — overlap between two Gaussian clusters
#                            (0 = no overlap, 1 = identical)
#     mahalanobis_dist     — scale-normalised centroid distance
#                            (accounts for cluster shape)
#     euclidean_dist       — raw centroid distance
#
#   PER-SAMPLE:
#     mahalanobis_score    — continuous distance from assigned
#                            cluster centroid (replaces binary
#                            OUTLIER flag for COMPASS-AI)
#
# Outputs:
#   {project}_pc_geometry.json       Full fingerprint (COMPASS-AI)
#   {project}_pc_geometry_sample.txt FID IID MAHAL_SCORE (per-sample)
#
# Usage:
#   python3 crisp_pc_geometry.py \
#       --eigenvec   sim_project_pca.eigenvec \
#       --satellite  sim_project_ancestry.txt \
#       --out-dir    step6_pca_pass1/ \
#       --project    sim_project \
#       --n-pcs      20
#
##########################################################################

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import kurtosis, skew


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6j — PC geometry fingerprint for COMPASS-AI"
    )
    p.add_argument("--eigenvec",  required=True,
                   help="PLINK .eigenvec file with PC scores")
    p.add_argument("--satellite", required=True,
                   help="Satellite ancestry file from Step 6f")
    p.add_argument("--out-dir",   required=True,
                   help="Output directory")
    p.add_argument("--project",   required=True,
                   help="Project name prefix")
    p.add_argument("--n-pcs",     type=int, default=20,
                   help="Number of PCs to include (default: 20)")
    return p.parse_args()


##########################################################################
# FILE READER
##########################################################################

def read_eigenvec(path: str, n_pcs: int) -> pd.DataFrame:
    """
    Read PLINK .eigenvec (PLINK 1.9 or PLINK 2 format).
    Returns FID, IID, PC1..PCn.
    """
    with open(path) as f:
        first = f.readline().strip()
    if first.startswith("#"):
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
    else:
        df = pd.read_csv(path, sep=r"\s+", header=None)
        n_cols = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i+1}" for i in range(n_cols)]
    pc_cols = [f"PC{i+1}" for i in range(n_pcs) if f"PC{i+1}" in df.columns]
    return df[["FID", "IID"] + pc_cols]


##########################################################################
# SAFE MATRIX OPERATIONS
##########################################################################

def safe_logdet(cov: np.ndarray) -> float:
    """
    Compute log determinant of a covariance matrix safely.
    WHY: det() of a high-dimensional covariance matrix underflows
    to zero for compact clusters. log-det avoids this and is
    numerically stable. We use the sign-logdet decomposition and
    return -inf for singular matrices rather than crashing.
    """
    try:
        sign, logdet = np.linalg.slogdet(cov)
        if sign <= 0:
            return float("-inf")
        return round(float(logdet), 6)
    except Exception:
        return float("-inf")


def safe_eigendecomp(cov: np.ndarray) -> tuple:
    """
    Eigendecomposition of a covariance matrix.
    Returns (eigenvalues, eigenvectors) sorted descending.
    WHY: eigh guarantees real eigenvalues for symmetric matrices.
    Small negative eigenvalues from floating point are clipped to 0.
    """
    try:
        vals, vecs = np.linalg.eigh(cov)
        vals = np.clip(vals, 0, None)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]
    except Exception:
        n = cov.shape[0]
        return np.zeros(n), np.eye(n)


def safe_inv(cov: np.ndarray) -> np.ndarray | None:
    """
    Safe matrix inverse — returns None if singular.
    WHY: Mahalanobis distance requires the inverse covariance.
    Singular matrices arise from degenerate clusters (n < n_pcs).
    We use pseudo-inverse as fallback.
    """
    try:
        return np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        try:
            return np.linalg.pinv(cov)
        except Exception:
            return None


##########################################################################
# COHORT-LEVEL METRICS
##########################################################################

def compute_cohort_metrics(X: np.ndarray,
                           n_pcs: int) -> dict:
    """
    Compute full cohort-level PC cloud geometry metrics.
    X: (n_samples, n_pcs) array of PC scores.
    """
    print(f"[CRISP] Computing cohort-level PC geometry "
          f"(n={X.shape[0]}, pcs={X.shape[1]})")

    metrics = {}

    # ── Mean and SD per PC ────────────────────────────────────────────
    pc_mean = X.mean(axis=0)
    pc_sd   = X.std(axis=0)
    metrics["pc_mean"] = [round(float(v), 8) for v in pc_mean]
    metrics["pc_sd"]   = [round(float(v), 8) for v in pc_sd]

    # ── Covariance matrix ─────────────────────────────────────────────
    # WHY: The full 20×20 covariance matrix is the richest single
    # descriptor of the cloud geometry. It encodes variance per axis,
    # correlations between axes, and the overall shape of the cloud.
    # COMPASS-AI can use it directly as input features (400 numbers
    # per cohort) or compute derived metrics from it.
    cov = np.cov(X, rowvar=False)
    metrics["covariance_matrix"] = [
        [round(float(v), 8) for v in row] for row in cov
    ]

    # ── Correlation matrix ────────────────────────────────────────────
    # WHY: The correlation matrix normalises the covariance by SD,
    # making it comparable across cohorts with different overall
    # variance scales. Off-diagonal elements reveal which PC axes
    # are correlated — unexpected correlations signal structure.
    std = np.sqrt(np.diag(cov))
    std[std == 0] = 1.0  # avoid division by zero
    corr = cov / np.outer(std, std)
    np.fill_diagonal(corr, 1.0)
    metrics["correlation_matrix"] = [
        [round(float(v), 6) for v in row] for row in corr
    ]

    # ── Eigendecomposition of covariance ─────────────────────────────
    # WHY: Eigenvalues of the covariance matrix reveal the principal
    # directions of variance in the cloud. Unlike the PLINK eigenvalues
    # (which reflect population-level allele frequency variance),
    # these eigenvalues reflect the within-cohort PC cloud shape.
    # A flat spectrum (all eigenvalues similar) = spherical cloud.
    # A steep spectrum (one dominant eigenvalue) = elongated cloud.
    vals, vecs = safe_eigendecomp(cov)
    metrics["covariance_eigenvalues"] = [round(float(v), 8) for v in vals]
    metrics["covariance_eigenvectors"] = [
        [round(float(v), 8) for v in row] for row in vecs
    ]

    # ── Cloud volume ──────────────────────────────────────────────────
    # WHY: Log-determinant of the covariance matrix = log of the
    # generalised variance = log of the volume of the ellipsoid
    # that best fits the cloud. Large volume = diverse cohort.
    # Small volume = compact, homogeneous cohort.
    metrics["cloud_log_volume"] = safe_logdet(cov)

    # ── Elongation ratio ──────────────────────────────────────────────
    # WHY: Ratio of largest to smallest eigenvalue of the covariance.
    # 1.0 = perfect sphere. Large ratio = sausage-shaped cloud
    # (one ancestry gradient dominates). Useful for COMPASS-AI to
    # detect cline-structured populations vs admixed clusters.
    vals_nonzero = vals[vals > 1e-12]
    if len(vals_nonzero) >= 2:
        elongation = float(vals_nonzero[0] / vals_nonzero[-1])
    else:
        elongation = 1.0
    metrics["elongation_ratio"] = round(elongation, 4)

    # ── Sphericity ────────────────────────────────────────────────────
    # WHY: Sphericity = geometric mean / arithmetic mean of eigenvalues.
    # 1.0 = perfect sphere (all axes equal variance). Low sphericity
    # = highly anisotropic cloud. Complements elongation_ratio which
    # only captures the extremes — sphericity reflects the full spectrum.
    if len(vals_nonzero) > 0:
        geo_mean   = float(np.exp(np.log(vals_nonzero).mean()))
        arith_mean = float(vals_nonzero.mean())
        sphericity = geo_mean / arith_mean if arith_mean > 0 else 0.0
    else:
        sphericity = 0.0
    metrics["sphericity"] = round(sphericity, 6)

    # ── Effective dimensionality ──────────────────────────────────────
    # WHY: How many PCs are needed to explain 90% of within-cohort
    # variance? Low effective dimensionality = most variation captured
    # by a few axes = dominant ancestry gradient. High effective
    # dimensionality = variance spread across many axes = complex
    # or admixed structure. A key COMPASS-AI signal for deciding
    # how many PCs to use as GWAS covariates.
    total_var = vals.sum()
    if total_var > 0:
        cumvar = np.cumsum(vals) / total_var
        eff_dim_90 = int(np.searchsorted(cumvar, 0.90)) + 1
        eff_dim_99 = int(np.searchsorted(cumvar, 0.99)) + 1
    else:
        eff_dim_90 = n_pcs
        eff_dim_99 = n_pcs
    metrics["effective_dimensionality_90pct"] = eff_dim_90
    metrics["effective_dimensionality_99pct"] = eff_dim_99
    metrics["variance_explained_cumulative"]  = [
        round(float(v), 6) for v in (np.cumsum(vals) / total_var if total_var > 0
                                      else np.zeros(len(vals)))
    ]

    # ── Mean nearest-neighbour distance ──────────────────────────────
    # WHY: Average distance from each sample to its nearest neighbour
    # in PC space. Dense clouds (small mean NN distance) indicate
    # a compact, homogeneous cohort. Sparse clouds indicate either
    # a diverse or under-sampled cohort. We subsample for efficiency
    # on large N — 2000 samples is sufficient for a stable estimate.
    n_nn_subsample = min(2000, X.shape[0])
    rng = np.random.default_rng(42)
    idx = rng.choice(X.shape[0], size=n_nn_subsample, replace=False)
    X_sub = X[idx]
    dists = cdist(X_sub, X_sub, metric="euclidean")
    np.fill_diagonal(dists, np.inf)
    nn_dists = dists.min(axis=1)
    metrics["mean_nn_distance"] = round(float(nn_dists.mean()), 8)
    metrics["sd_nn_distance"]   = round(float(nn_dists.std()),  8)

    # ── Kurtosis and skewness per PC ─────────────────────────────────
    # WHY: Kurtosis measures tail heaviness. High excess kurtosis on
    # a PC axis suggests outliers or bimodal structure on that axis —
    # a COMPASS-AI signal for potential sub-structure or contamination.
    # Skewness detects asymmetric distributions — expected when a
    # cohort contains one large and one small ancestry group.
    kurt = [round(float(kurtosis(X[:, i], fisher=True)), 4)
            for i in range(X.shape[1])]
    skewness = [round(float(skew(X[:, i])), 4)
                for i in range(X.shape[1])]
    metrics["kurtosis_per_pc"]  = kurt
    metrics["skewness_per_pc"]  = skewness

    print(f"[CRISP] Cloud log-volume    : {metrics['cloud_log_volume']}")
    print(f"[CRISP] Elongation ratio    : {metrics['elongation_ratio']}")
    print(f"[CRISP] Sphericity          : {metrics['sphericity']}")
    print(f"[CRISP] Effective dim (90%) : {eff_dim_90}")
    print(f"[CRISP] Mean NN distance    : {metrics['mean_nn_distance']}")

    return metrics


##########################################################################
# PER-GROUP METRICS
##########################################################################

def compute_group_metrics(X: np.ndarray,
                          groups: np.ndarray,
                          all_groups: list) -> dict:
    """
    Compute per-GMM-cluster geometry metrics.
    WHY: Cohort-level metrics describe the full cloud. Per-group
    metrics describe each ancestry cluster independently. COMPASS-AI
    needs both — two cohorts can have identical overall cloud volumes
    but very different within-group compactness.
    """
    print(f"[CRISP] Computing per-group PC geometry ({len(all_groups)} groups)")
    group_metrics = {}

    for grp in all_groups:
        mask  = groups == grp
        X_grp = X[mask]
        n_grp = int(mask.sum())

        gm = {"n": n_grp}

        if n_grp < 3:
            # Too few samples for meaningful geometry
            gm["centroid"]   = None
            gm["covariance"] = None
            gm["log_volume"] = None
            gm["elongation"] = None
            gm["density"]    = None
            group_metrics[grp] = gm
            continue

        centroid = X_grp.mean(axis=0)
        gm["centroid"] = [round(float(v), 8) for v in centroid]

        cov = np.cov(X_grp, rowvar=False) if n_grp > 1 else np.zeros((X.shape[1], X.shape[1]))
        gm["covariance"] = [[round(float(v), 8) for v in row] for row in cov]
        gm["log_volume"] = safe_logdet(cov)

        vals, _ = safe_eigendecomp(cov)
        vals_nz = vals[vals > 1e-12]
        gm["elongation"] = round(float(vals_nz[0] / vals_nz[-1]), 4) if len(vals_nz) >= 2 else 1.0

        # Within-group density
        if n_grp >= 2:
            d = cdist(X_grp, X_grp, metric="euclidean")
            np.fill_diagonal(d, np.inf)
            gm["density"] = round(float(d.min(axis=1).mean()), 8)
        else:
            gm["density"] = None

        group_metrics[grp] = gm
        print(f"[CRISP]   {grp}: n={n_grp}  logvol={gm['log_volume']}  elong={gm['elongation']}")

    return group_metrics


##########################################################################
# INTER-GROUP METRICS
##########################################################################

def bhattacharyya_coefficient(mu1: np.ndarray, cov1: np.ndarray,
                               mu2: np.ndarray, cov2: np.ndarray) -> float:
    """
    Bhattacharyya coefficient between two multivariate Gaussians.
    Measures overlap: 0 = no overlap, 1 = identical distributions.

    WHY: Standard pairwise distance metrics (Euclidean, Mahalanobis)
    measure centroid separation but ignore cluster shape. Two clusters
    can be far apart in centroid distance but heavily overlapping if
    they have large covariance. The Bhattacharyya coefficient accounts
    for both centroid distance and covariance structure simultaneously.
    It is the natural measure of cluster overlap for COMPASS-AI — a
    high coefficient means the two groups are genetically similar and
    may not be meaningfully distinct for stratification purposes.

    Formula:
      BC = exp(-D_B)
      D_B = (1/8)(mu1-mu2)^T * S^-1 * (mu1-mu2)
            + (1/2) * ln(det(S) / sqrt(det(cov1)*det(cov2)))
      where S = (cov1 + cov2) / 2
    """
    try:
        diff = mu1 - mu2
        S    = (cov1 + cov2) / 2.0
        S_inv = safe_inv(S)
        if S_inv is None:
            return float("nan")

        # Mahalanobis term
        maha_term = 0.125 * float(diff @ S_inv @ diff)

        # Log-determinant term
        _, ld_S    = np.linalg.slogdet(S)
        _, ld_cov1 = np.linalg.slogdet(cov1)
        _, ld_cov2 = np.linalg.slogdet(cov2)
        det_term = 0.5 * (ld_S - 0.5 * (ld_cov1 + ld_cov2))

        D_B = maha_term + det_term
        BC  = float(np.exp(-D_B))
        return round(max(0.0, min(1.0, BC)), 6)

    except Exception:
        return float("nan")


def mahalanobis_centroid_dist(mu1: np.ndarray, cov1: np.ndarray,
                               mu2: np.ndarray, cov2: np.ndarray) -> float:
    """
    Mahalanobis distance between two cluster centroids.
    Uses the pooled covariance matrix as the metric.

    WHY: Raw Euclidean distance between centroids is scale-dependent
    — a distance of 0.01 in PC1 means something very different from
    0.01 in PC10 (which has much less variance). Mahalanobis distance
    normalises by the pooled covariance, giving a scale-free measure
    of centroid separation that is comparable across cohorts.
    """
    try:
        S_pooled = (cov1 + cov2) / 2.0
        S_inv = safe_inv(S_pooled)
        if S_inv is None:
            return float("nan")
        diff = mu1 - mu2
        maha = float(np.sqrt(max(0.0, diff @ S_inv @ diff)))
        return round(maha, 6)
    except Exception:
        return float("nan")


def compute_intergroup_metrics(group_metrics: dict,
                               all_groups: list) -> dict:
    """
    Compute pairwise inter-group geometry metrics for all group pairs.
    """
    print(f"[CRISP] Computing inter-group metrics")
    inter = {}

    for i, g1 in enumerate(all_groups):
        for g2 in all_groups[i+1:]:
            gm1 = group_metrics[g1]
            gm2 = group_metrics[g2]

            pair_key = f"{g1}_vs_{g2}"
            pair = {}

            # Skip if either group has insufficient data
            if gm1["centroid"] is None or gm2["centroid"] is None:
                pair["bhattacharyya_coeff"] = None
                pair["mahalanobis_dist"]    = None
                pair["euclidean_dist"]      = None
                inter[pair_key] = pair
                continue

            mu1  = np.array(gm1["centroid"])
            mu2  = np.array(gm2["centroid"])
            cov1 = np.array(gm1["covariance"])
            cov2 = np.array(gm2["covariance"])

            # Euclidean centroid distance
            pair["euclidean_dist"] = round(float(np.linalg.norm(mu1 - mu2)), 6)

            # Bhattacharyya coefficient
            pair["bhattacharyya_coeff"] = bhattacharyya_coefficient(mu1, cov1, mu2, cov2)

            # Mahalanobis centroid distance
            pair["mahalanobis_dist"] = mahalanobis_centroid_dist(mu1, cov1, mu2, cov2)

            inter[pair_key] = pair
            print(f"[CRISP]   {pair_key}: "
                  f"BC={pair['bhattacharyya_coeff']}  "
                  f"Maha={pair['mahalanobis_dist']}  "
                  f"Euc={pair['euclidean_dist']}")

    return inter


##########################################################################
# PER-SAMPLE MAHALANOBIS SCORE
##########################################################################

def compute_sample_mahalanobis(X: np.ndarray,
                                groups: np.ndarray,
                                group_metrics: dict,
                                iids: np.ndarray,
                                fids: np.ndarray) -> pd.DataFrame:
    """
    Compute per-sample Mahalanobis distance from assigned cluster centroid.

    WHY: The binary OUTLIER flag from Step 6d is a hard threshold.
    COMPASS-AI needs a continuous score — a sample with Mahalanobis
    score of 2.1 is very different from one with a score of 8.7,
    even if both are below the OUTLIER_SD threshold. This continuous
    score is richer training signal and enables COMPASS-AI to learn
    non-linear relationships between cluster position and downstream
    QC outcomes.

    Also: samples from small or admixed clusters may have inflated
    Mahalanobis scores not because they are genuine outliers but
    because their cluster's covariance estimate is unstable (small n).
    We flag this with a reliability indicator.
    """
    print(f"[CRISP] Computing per-sample Mahalanobis scores")

    all_groups = list(group_metrics.keys())
    scores = []
    reliable = []

    for i, (iid, fid, grp) in enumerate(zip(iids, fids, groups)):
        gm = group_metrics.get(grp, {})

        if gm.get("centroid") is None or gm.get("covariance") is None:
            scores.append(float("nan"))
            reliable.append(False)
            continue

        mu  = np.array(gm["centroid"])
        cov = np.array(gm["covariance"])
        n   = gm.get("n", 0)

        cov_inv = safe_inv(cov)
        if cov_inv is None:
            scores.append(float("nan"))
            reliable.append(False)
            continue

        diff  = X[i] - mu
        maha  = float(np.sqrt(max(0.0, diff @ cov_inv @ diff)))
        scores.append(round(maha, 6))

        # WHY reliability flag: Mahalanobis distance from a cluster
        # with n < 2 * n_pcs samples has an unreliable covariance
        # estimate — the matrix may be near-singular. We flag these
        # so COMPASS-AI can down-weight them during training.
        reliable.append(n >= 2 * X.shape[1])

    result = pd.DataFrame({
        "FID"          : fids,
        "IID"          : iids,
        "GROUP"        : groups,
        "MAHAL_SCORE"  : scores,
        "MAHAL_RELIABLE": ["YES" if r else "NO" for r in reliable],
    })

    print(f"[CRISP] Mahalanobis scores: "
          f"mean={np.nanmean(scores):.3f}  "
          f"max={np.nanmax(scores):.3f}  "
          f"reliable={sum(reliable)}/{len(reliable)}")

    return result


##########################################################################
# MAIN
##########################################################################

def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP] crisp_pc_geometry.py — {args.project}")
    print(f"[CRISP] N PCs: {args.n_pcs}")

    # ── Load data ─────────────────────────────────────────────────────
    eigenvec  = read_eigenvec(args.eigenvec, args.n_pcs)
    satellite = pd.read_csv(args.satellite, sep="\t")

    pc_cols = [f"PC{i+1}" for i in range(args.n_pcs)
               if f"PC{i+1}" in eigenvec.columns]
    n_pcs_actual = len(pc_cols)

    # Merge satellite GROUP labels with PC scores
    df = satellite[["FID", "IID", "SUPERPOP", "OUTLIER"]].merge(
        eigenvec[["IID"] + pc_cols], on="IID", how="inner"
    ).dropna(subset=pc_cols)

    print(f"[CRISP] Samples loaded: {len(df)}")
    print(f"[CRISP] PCs available : {n_pcs_actual}")

    X      = df[pc_cols].values.astype(float)
    groups = df["SUPERPOP"].values
    iids   = df["IID"].values
    fids   = df["FID"].values

    all_groups = sorted(set(groups))

    # ── Cohort-level metrics ──────────────────────────────────────────
    cohort_metrics = compute_cohort_metrics(X, n_pcs_actual)

    # ── Per-group metrics ─────────────────────────────────────────────
    group_metrics = compute_group_metrics(X, groups, all_groups)

    # ── Inter-group metrics ───────────────────────────────────────────
    inter_metrics = compute_intergroup_metrics(group_metrics, all_groups)

    # ── Per-sample Mahalanobis scores ─────────────────────────────────
    sample_scores = compute_sample_mahalanobis(
        X, groups, group_metrics, iids, fids
    )

    # Write per-sample scores
    sample_path = out_dir / f"{args.project}_pc_geometry_sample.txt"
    sample_scores.to_csv(sample_path, sep="\t", index=False)
    print(f"[CRISP] Per-sample scores: {sample_path}")

    # ── Assemble full fingerprint JSON ────────────────────────────────
    # WHY: Every metric goes into a single versioned JSON. No metric
    # is omitted on grounds of seeming unnecessary — COMPASS-AI
    # decides what is useful during training. The schema_version
    # field allows future versions to handle format changes gracefully.
    fingerprint = {
        "schema_version"     : "1.0.0",
        "sub_step"           : "6j_pc_geometry",
        "project"            : args.project,
        "n_samples"          : int(len(df)),
        "n_pcs"              : n_pcs_actual,
        "n_groups"           : len(all_groups),
        "groups"             : all_groups,

        # Full cohort cloud geometry
        "cohort"             : cohort_metrics,

        # Per-group geometry (one block per GMM cluster)
        "groups_geometry"    : group_metrics,

        # Pairwise inter-group metrics
        "inter_group"        : inter_metrics,

        # Summary derived signals for COMPASS-AI quick access
        "compass_signals"    : {
            "cloud_log_volume"         : cohort_metrics["cloud_log_volume"],
            "elongation_ratio"         : cohort_metrics["elongation_ratio"],
            "sphericity"               : cohort_metrics["sphericity"],
            "effective_dim_90pct"      : cohort_metrics["effective_dimensionality_90pct"],
            "mean_nn_distance"         : cohort_metrics["mean_nn_distance"],
            "n_groups"                 : len(all_groups),
            "most_isolated_group"      : (
                max(inter_metrics,
                    key=lambda k: inter_metrics[k].get("mahalanobis_dist") or 0)
                if inter_metrics else None
            ),
            "most_overlapping_pair"    : (
                max(inter_metrics,
                    key=lambda k: inter_metrics[k].get("bhattacharyya_coeff") or 0)
                if inter_metrics else None
            ),
            "max_group_elongation"     : max(
                (gm["elongation"] or 0 for gm in group_metrics.values()),
                default=None
            ),
        },
    }

    json_path = out_dir / f"{args.project}_pc_geometry.json"
    with open(json_path, "w") as f:
        json.dump(fingerprint, f, indent=2)
    print(f"[CRISP] PC geometry JSON: {json_path}")

    # ── Summary ───────────────────────────────────────────────────────
    print(f"[CRISP] ══════════════════════════════════════════════════")
    print(f"[CRISP] PC geometry fingerprint complete")
    print(f"[CRISP] Samples          : {len(df)}")
    print(f"[CRISP] Groups           : {all_groups}")
    print(f"[CRISP] Cloud log-volume : {cohort_metrics['cloud_log_volume']}")
    print(f"[CRISP] Elongation ratio : {cohort_metrics['elongation_ratio']}")
    print(f"[CRISP] Effective dim    : {cohort_metrics['effective_dimensionality_90pct']} PCs → 90% variance")
    print(f"[CRISP] crisp_pc_geometry.py complete.")


if __name__ == "__main__":
    main()
