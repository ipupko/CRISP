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
# Script : scripts/crisp_gmm_cluster.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   GMM (Gaussian Mixture Model) clustering of PCA or MDS output,
#   followed by an SD-based outlier check within each cluster.
#
#   Two-method design rationale:
#     GMM  — answers "which group does this sample belong to?"
#             Returns a soft assignment with posterior probability
#             (MEMBERSHIP_P). Handles elliptical clusters. Admixed
#             samples that sit between clusters get low MEMBERSHIP_P,
#             making them visible without hard exclusion.
#     SD   — answers "is this sample a clean member of its group?"
#             Flags samples > PCA_OUTLIER_SD standard deviations from
#             their cluster centroid as OUTLIER=YES. Distinct from
#             the GMM assignment — a sample can be assigned to a group
#             but still sit at its outer edge.
#
#   AUTO cluster detection:
#     When PCA_N_CLUSTERS=AUTO, we fit GMM for K=1..8 and select K
#     via BIC (Bayesian Information Criterion). BIC penalises complexity
#     so it naturally avoids over-splitting. The BIC curve and selected
#     K are written to the summary JSON for analyst review.
#
#   1KG projection mode:
#     When PCA_MODE=1KG_PROJECTION, the cohort was merged with 1KG
#     before PCA. 1KG samples are used to anchor cluster centroids to
#     named superpopulations (EUR, AFR, EAS, SAS, AMR). Cohort samples
#     are then assigned the superpopulation label of their nearest
#     centroid. Without 1KG, clusters are labelled GROUP_1, GROUP_2 etc.
#
# Outputs (written to --out-dir):
#   {project}_gmm_assignments.txt   FID IID GROUP MEMBERSHIP_P OUTLIER
#   {project}_cluster_{N}.txt       FID IID per-cluster ID lists
#   {project}_gmm_summary.json      K, BIC scores, cluster sizes
#
# Usage:
#   python3 crisp_gmm_cluster.py \
#       --eigenvec  project.eigenvec \
#       --n-clusters AUTO \
#       --outlier-sd 6 \
#       --pca-mode  COHORT_ONLY \
#       --kg-prefix /path/to/1kg \
#       --out-dir   step6_pca_pass1/ \
#       --project   my_project
#
##########################################################################

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

# Suppress sklearn convergence warnings in logs — we handle them gracefully
warnings.filterwarnings("ignore", category=UserWarning)


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6 — GMM clustering + SD outlier detection"
    )
    p.add_argument("--eigenvec",    required=True,
                   help="PLINK .eigenvec file (PCA) or .mds file (MDS)")
    p.add_argument("--n-clusters",  default="AUTO",
                   help="Number of clusters (AUTO or integer, default: AUTO)")
    p.add_argument("--outlier-sd",  type=float, default=6.0,
                   help="SD threshold for outlier flagging within clusters (default: 6)")
    p.add_argument("--pca-mode",    default="COHORT_ONLY",
                   choices=["COHORT_ONLY", "1KG_PROJECTION", "MDS"],
                   help="PCA computation mode")
    p.add_argument("--kg-prefix",   default="",
                   help="Path to 1KG PLINK prefix (needed for 1KG_PROJECTION mode)")
    p.add_argument("--out-dir",     required=True,
                   help="Output directory")
    p.add_argument("--project",     required=True,
                   help="Project name prefix for output files")
    return p.parse_args()


##########################################################################
# FILE READERS
##########################################################################

def read_eigenvec(path: str) -> pd.DataFrame:
    """
    Read a PLINK .eigenvec file into a DataFrame.
    PLINK 1.9 format: FID IID PC1 PC2 ... PCn (space-delimited, no header)
    PLINK 2 format  : #FID IID PC1 PC2 ...    (space-delimited, header with #FID)
    WHY: We auto-detect format by checking whether the first line starts
    with '#' (PLINK 2) or not (PLINK 1.9), so the script works with both
    PLINK versions without manual configuration.
    """
    with open(path) as f:
        first = f.readline().strip()

    if first.startswith("#"):
        # PLINK 2 — strip the # from #FID
        df = pd.read_csv(path, sep=r"\s+", comment=None)
        df.columns = [c.lstrip("#") for c in df.columns]
    else:
        # PLINK 1.9 — no header, columns are FID IID PC1 PC2 ...
        df = pd.read_csv(path, sep=r"\s+", header=None)
        n_pcs = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i+1}" for i in range(n_pcs)]

    return df


def read_mds(path: str) -> pd.DataFrame:
    """
    Read a PLINK .mds file into a DataFrame.
    Format: FID IID SOL C1 C2 ... Cn (space-delimited, header present)
    WHY: MDS output has an extra SOL column (solution number) that we
    drop, then rename C1..Cn to PC1..PCn so downstream code is unified
    between PCA and MDS processing paths.
    """
    df = pd.read_csv(path, sep=r"\s+")
    # Drop SOL column if present
    if "SOL" in df.columns:
        df = df.drop(columns=["SOL"])
    # Rename C1..Cn to PC1..PCn for unified handling
    rename = {c: c.replace("C", "PC") for c in df.columns if c.startswith("C")}
    df = df.rename(columns=rename)
    return df


def read_kg_populations(kg_prefix: str) -> pd.DataFrame | None:
    """
    Read 1KG population annotation file alongside the PLINK prefix.
    Looks for: {prefix}.pop (FID IID SUPERPOP) or {prefix}.psam
    (PLINK2 format with SuperPop or Population column).
    WHY: When 1KG projection is used, we need to know which 1KG samples
    belong to which superpopulation so we can name cluster centroids.
    Without this file, clusters remain as GROUP_1, GROUP_2 etc. even in
    projection mode.
    """
    if not kg_prefix:
        return None

    # Try .pop first (simple three-column format)
    pop_path = kg_prefix + ".pop"
    if os.path.exists(pop_path):
        df = pd.read_csv(pop_path, sep=r"\s+", header=None,
                         names=["FID", "IID", "SUPERPOP"])
        return df[["IID", "SUPERPOP"]]

    # Try .psam (PLINK 2 sample annotation)
    psam_path = kg_prefix + ".psam"
    if os.path.exists(psam_path):
        df = pd.read_csv(psam_path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
        for col in ["SuperPop", "Population", "SUPERPOP", "pop"]:
            if col in df.columns:
                df = df.rename(columns={col: "SUPERPOP"})
                break
        if "IID" in df.columns and "SUPERPOP" in df.columns:
            return df[["IID", "SUPERPOP"]]

    return None


##########################################################################
# AUTO K SELECTION VIA BIC
##########################################################################

def select_k_bic(X: np.ndarray,
                 k_range: range = range(1, 9),
                 n_init: int = 10,
                 random_state: int = 42) -> tuple[int, dict]:
    """
    Fit GMM for each K in k_range, return the K with lowest BIC.
    WHY: BIC (Bayesian Information Criterion) penalises model complexity.
    For population structure, it naturally selects the smallest K that
    explains the variance — avoiding over-splitting small cohorts into
    spurious sub-groups while correctly detecting real structure in
    larger multi-ancestry cohorts.
    Returns (best_k, {k: bic_score}) for analyst review.
    """
    bic_scores = {}
    best_k = 1
    best_bic = np.inf

    print(f"[CRISP] AUTO K selection via BIC (K=1..{max(k_range)})")

    for k in k_range:
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                n_init=n_init,
                random_state=random_state,
                max_iter=500,
            )
            gmm.fit(X)
            bic = gmm.bic(X)
            bic_scores[k] = round(float(bic), 2)
            print(f"[CRISP]   K={k}  BIC={bic:.1f}")

            if bic < best_bic:
                best_bic = bic
                best_k = k

        except Exception as e:
            print(f"[CRISP]   K={k}  FAILED: {e}", file=sys.stderr)
            bic_scores[k] = None

    print(f"[CRISP] Selected K={best_k} (BIC={best_bic:.1f})")
    return best_k, bic_scores


##########################################################################
# GMM FITTING
##########################################################################

def fit_gmm(X: np.ndarray,
            k: int,
            n_init: int = 10,
            random_state: int = 42) -> GaussianMixture:
    """
    Fit final GMM with the chosen K.
    WHY: n_init=10 runs the EM algorithm from 10 random initialisations
    and keeps the best result (highest log-likelihood). This avoids
    getting stuck in local optima, which is a known problem for GMM with
    elliptical clusters in high-dimensional PC space.
    """
    gmm = GaussianMixture(
        n_components=k,
        covariance_type="full",
        n_init=n_init,
        random_state=random_state,
        max_iter=500,
    )
    gmm.fit(X)
    return gmm


##########################################################################
# SD OUTLIER CHECK WITHIN CLUSTERS
##########################################################################

def flag_outliers(X: np.ndarray,
                  labels: np.ndarray,
                  sd_threshold: float) -> np.ndarray:
    """
    For each cluster, compute the centroid and per-sample Euclidean
    distance from centroid in units of SD. Flag samples beyond
    sd_threshold as outliers.
    WHY: GMM assigns every sample to some cluster — even extreme outliers
    get assigned to the nearest cluster. The SD check catches samples that
    are assigned to a cluster but sit far from its centre. These are
    either genuine ancestry outliers (not fitting any cohort group cleanly)
    or admixed samples at the boundary between two groups.
    Returns boolean array: True = outlier.
    """
    n_samples = X.shape[0]
    outlier_flags = np.zeros(n_samples, dtype=bool)

    for cluster_id in np.unique(labels):
        mask = labels == cluster_id
        X_cluster = X[mask]

        if X_cluster.shape[0] < 3:
            # WHY: Cannot compute meaningful SD with fewer than 3 samples.
            # Mark all as non-outliers to avoid false positives in
            # tiny clusters that may be legitimate rare-ancestry samples.
            continue

        centroid = X_cluster.mean(axis=0)
        distances = np.linalg.norm(X_cluster - centroid, axis=1)

        # SD-based threshold
        sd = distances.std()
        if sd == 0:
            continue

        threshold = sd_threshold * sd
        cluster_outliers = distances > threshold

        # Map back to full array indices
        cluster_indices = np.where(mask)[0]
        outlier_flags[cluster_indices] = cluster_outliers

    return outlier_flags


##########################################################################
# CLUSTER LABEL ASSIGNMENT
##########################################################################

def assign_labels(labels: np.ndarray,
                  df: pd.DataFrame,
                  pc_cols: list[str],
                  kg_pops: pd.DataFrame | None,
                  pca_mode: str) -> list[str]:
    """
    Assign human-readable group labels to GMM cluster IDs.
    WHY:
      COHORT_ONLY / MDS mode — no reference anchor, so clusters get
        generic labels: GROUP_1, GROUP_2 etc. The analyst can visually
        inspect the PCA/MDS plots and map these to known populations.
      1KG_PROJECTION mode — 1KG samples are present in the dataset with
        known superpopulation labels. We find the modal superpopulation
        within each GMM cluster (using only 1KG samples) and assign that
        name to the cluster. Cohort samples then inherit the superpop
        label of their cluster.
    """
    n_clusters = len(np.unique(labels))

    if pca_mode != "1KG_PROJECTION" or kg_pops is None:
        # Generic labels
        label_map = {i: f"GROUP_{i+1}" for i in range(n_clusters)}
        return [label_map[l] for l in labels]

    # 1KG projection — find modal superpop per cluster
    df_tmp = df.copy()
    df_tmp["_CLUSTER"] = labels
    merged = df_tmp.merge(kg_pops, on="IID", how="left")

    label_map = {}
    used_pops = set()

    for cluster_id in sorted(np.unique(labels)):
        cluster_mask = merged["_CLUSTER"] == cluster_id
        kg_in_cluster = merged.loc[cluster_mask, "SUPERPOP"].dropna()

        if len(kg_in_cluster) == 0:
            # No 1KG samples in this cluster — fall back to generic label
            label = f"GROUP_{cluster_id + 1}"
        else:
            modal_pop = kg_in_cluster.value_counts().idxmax()
            # Avoid duplicate superpop labels if two clusters have same modal pop
            if modal_pop in used_pops:
                label = f"{modal_pop}_2"
            else:
                label = modal_pop
                used_pops.add(modal_pop)

        label_map[cluster_id] = label
        print(f"[CRISP]   Cluster {cluster_id} → {label} "
              f"(n_1KG={len(kg_in_cluster)})")

    return [label_map[l] for l in labels]


##########################################################################
# MAIN
##########################################################################

def main():
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP] crisp_gmm_cluster.py — {args.project}")
    print(f"[CRISP] Input   : {args.eigenvec}")
    print(f"[CRISP] Mode    : {args.pca_mode}")
    print(f"[CRISP] Clusters: {args.n_clusters}")
    print(f"[CRISP] OutlierSD: {args.outlier_sd}")

    # ── Read PC / MDS scores ──────────────────────────────────────────
    is_mds = args.pca_mode == "MDS" or args.eigenvec.endswith(".mds")

    if is_mds:
        df = read_mds(args.eigenvec)
    else:
        df = read_eigenvec(args.eigenvec)

    pc_cols = [c for c in df.columns if c.startswith("PC")]
    n_pcs_available = len(pc_cols)

    # WHY: For clustering we use the first 10 PCs by default (or all if
    # fewer are available). PCs beyond 10 typically capture noise rather
    # than population structure and can degrade GMM performance.
    n_pcs_for_clustering = min(10, n_pcs_available)
    cluster_cols = pc_cols[:n_pcs_for_clustering]

    print(f"[CRISP] Samples   : {len(df)}")
    print(f"[CRISP] PCs avail : {n_pcs_available}  (using first {n_pcs_for_clustering} for clustering)")

    # ── Separate 1KG samples from cohort samples ──────────────────────
    # WHY: In 1KG projection mode, the merged dataset contains both cohort
    # and 1KG samples. We need to cluster all samples together (so cohort
    # samples are placed in the 1KG PC space) but only report assignments
    # for the cohort samples. 1KG samples are identified by their IID
    # appearing in the 1KG population annotation file.

    kg_pops = read_kg_populations(args.kg_prefix)

    if args.pca_mode == "1KG_PROJECTION" and kg_pops is not None:
        is_kg = df["IID"].isin(kg_pops["IID"])
        df_cohort = df[~is_kg].copy()
        print(f"[CRISP] 1KG samples  : {is_kg.sum()}  (used for centroid anchoring)")
        print(f"[CRISP] Cohort samples: {(~is_kg).sum()}")
    else:
        df_cohort = df.copy()
        is_kg = pd.Series(False, index=df.index)

    # ── Scale features ────────────────────────────────────────────────
    # WHY: GMM is sensitive to the scale of input dimensions. PC variance
    # decreases with PC number (PC1 explains the most variance, PC10 the
    # least). Without scaling, PC1 dominates the clustering and later PCs
    # are effectively ignored. StandardScaler brings all PCs to unit
    # variance so they contribute equally to the cluster geometry.
    X_all = df[cluster_cols].values
    scaler = StandardScaler()
    X_all_scaled = scaler.fit_transform(X_all)

    # Cohort-only indices for later extraction
    cohort_idx = np.where(~is_kg)[0] if args.pca_mode == "1KG_PROJECTION" else np.arange(len(df))
    X_cohort_scaled = X_all_scaled[cohort_idx]

    # ── Select K ──────────────────────────────────────────────────────
    bic_scores = {}

    if args.n_clusters.upper() == "AUTO":
        k, bic_scores = select_k_bic(X_all_scaled)
    else:
        try:
            k = int(args.n_clusters)
            print(f"[CRISP] Using fixed K={k}")
        except ValueError:
            print(f"[CRISP] ERROR: --n-clusters must be AUTO or integer, got: {args.n_clusters}",
                  file=sys.stderr)
            sys.exit(1)

    # ── Fit GMM ───────────────────────────────────────────────────────
    print(f"[CRISP] Fitting GMM (K={k})...")
    gmm = fit_gmm(X_all_scaled, k)

    # Get labels and posterior probabilities for all samples
    all_labels = gmm.predict(X_all_scaled)
    all_probs  = gmm.predict_proba(X_all_scaled)
    # MEMBERSHIP_P = probability of the assigned cluster (max posterior)
    all_membership_p = all_probs.max(axis=1)

    # Cohort-only assignments
    cohort_labels       = all_labels[cohort_idx]
    cohort_membership_p = all_membership_p[cohort_idx]

    # ── SD outlier check within clusters ──────────────────────────────
    print(f"[CRISP] SD outlier check (threshold={args.outlier_sd} SD within each cluster)...")
    # WHY: SD check is performed on UNSCALED PC coordinates so the SD
    # threshold is interpretable in terms of the original PC space.
    # Scaling would compress distances and make threshold comparison
    # to the unscaled eigenvalue-weighted PC space misleading.
    X_cohort_raw = df_cohort[cluster_cols].values
    outlier_flags = flag_outliers(X_cohort_raw, cohort_labels, args.outlier_sd)

    n_outliers = outlier_flags.sum()
    print(f"[CRISP] Outliers flagged: {n_outliers} / {len(df_cohort)}")

    # ── Assign human-readable labels ──────────────────────────────────
    group_labels = assign_labels(
        cohort_labels, df_cohort, cluster_cols, kg_pops, args.pca_mode
    )

    # ── Build output DataFrame ────────────────────────────────────────
    result = df_cohort[["FID", "IID"]].copy().reset_index(drop=True)
    result["GROUP"]        = group_labels
    result["MEMBERSHIP_P"] = np.round(cohort_membership_p, 4)
    result["OUTLIER"]      = ["YES" if o else "NO" for o in outlier_flags]

    # ── Write assignment file ─────────────────────────────────────────
    assign_path = out_dir / f"{args.project}_gmm_assignments.txt"
    result.to_csv(assign_path, sep="\t", index=False)
    print(f"[CRISP] Assignments written: {assign_path}")

    # ── Write per-cluster sample ID files ─────────────────────────────
    # WHY: These are consumed by crisp_pca_pass1.sh to run plink --keep
    # for each group — mirroring the DIVERGE pattern where separate
    # per-province datasets are extracted for downstream analyses.
    unique_groups = sorted(set(group_labels))
    cluster_sizes = {}

    for group in unique_groups:
        mask = result["GROUP"] == group
        cluster_df = result.loc[mask, ["FID", "IID"]]
        # Extract cluster number for filename
        cluster_n = group.replace("GROUP_", "").replace("_", "")
        cluster_path = out_dir / f"{args.project}_cluster_{cluster_n}.txt"
        cluster_df.to_csv(cluster_path, sep="\t", index=False, header=False)
        cluster_sizes[group] = int(mask.sum())
        print(f"[CRISP]   {group}: {mask.sum()} samples → {cluster_path.name}")

    # ── Write JSON summary ─────────────────────────────────────────────
    summary = {
        "schema_version"   : "1.0.0",
        "project"          : args.project,
        "pca_mode"         : args.pca_mode,
        "n_clusters"       : k,
        "k_selection"      : "AUTO_BIC" if args.n_clusters.upper() == "AUTO" else "FIXED",
        "bic_scores"       : bic_scores,
        "n_pcs_used"       : n_pcs_for_clustering,
        "n_samples_total"  : len(df),
        "n_samples_cohort" : len(df_cohort),
        "n_outliers"       : int(n_outliers),
        "outlier_sd"       : args.outlier_sd,
        "cluster_sizes"    : cluster_sizes,
        "groups"           : unique_groups,
        "gmm_converged"    : bool(gmm.converged_),
        "gmm_n_iter"       : int(gmm.n_iter_),
        "gmm_log_likelihood": round(float(gmm.lower_bound_), 4),
    }

    json_path = out_dir / f"{args.project}_gmm_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[CRISP] Summary JSON written: {json_path}")
    print(f"[CRISP] crisp_gmm_cluster.py complete.")


if __name__ == "__main__":
    main()
