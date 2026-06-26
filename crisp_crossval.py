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
# Script : scripts/crisp_crossval.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Cross-validate cluster assignments between PCA-GMM and MDS-GMM.
#
#   PCA and MDS approach population structure from different angles:
#     PCA — eigendecomposition of the allele frequency covariance matrix
#     MDS — multidimensional scaling of pairwise IBS distances
#   For a well-structured cohort the two methods should produce very
#   similar cluster assignments. High disagreement signals either:
#     a) Complex admixture that one method captures better than the other
#     b) An unstable K choice — the structure is ambiguous
#     c) A technical artefact in one of the two analyses
#
#   Agreement is computed using the Adjusted Rand Index (ARI), which
#   measures cluster assignment concordance adjusted for chance. ARI=1
#   means perfect agreement, ARI=0 means no better than random.
#   We also report raw per-sample agreement (YES/NO) for the satellite
#   file and a percentage agreement rate for the report.
#
#   WHY Hungarian matching:
#     GMM cluster labels are arbitrary integers — cluster 0 in PCA may
#     correspond to cluster 2 in MDS. Before computing agreement, we
#     must find the optimal label mapping between the two assignments.
#     The Hungarian algorithm (scipy.optimize.linear_sum_assignment)
#     finds this mapping in O(n^3) time by maximising the overlap
#     between cluster pairs.
#
# Outputs:
#   {project}_mds_agree.txt          FID IID MDS_AGREE
#   {project}_crossval_summary.json  agreement_rate_pct, ARI, mapping
#
# Usage:
#   python3 crisp_crossval.py \
#       --pca-assign  sim_project_gmm_assignments.txt \
#       --mds-assign  sim_project_mds_gmm_assignments.txt \
#       --out         sim_project_mds_agree.txt \
#       --project     sim_project
#
##########################################################################

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6 — PCA/MDS cluster assignment cross-validation"
    )
    p.add_argument("--pca-assign", required=True,
                   help="PCA GMM assignment file ({project}_gmm_assignments.txt)")
    p.add_argument("--mds-assign", required=True,
                   help="MDS GMM assignment file ({project}_mds_gmm_assignments.txt)")
    p.add_argument("--out",        required=True,
                   help="Output per-sample agreement file ({project}_mds_agree.txt)")
    p.add_argument("--project",    required=True,
                   help="Project name prefix for JSON summary")
    return p.parse_args()


##########################################################################
# HUNGARIAN LABEL MATCHING
##########################################################################

def match_labels(labels_a: np.ndarray,
                 labels_b: np.ndarray) -> dict:
    """
    Find the optimal mapping from labels_b to labels_a that maximises
    cluster overlap, using the Hungarian algorithm.

    WHY: GMM assigns cluster IDs as arbitrary integers. The same
    population cluster may be labelled GROUP_1 in the PCA run and
    GROUP_3 in the MDS run. Without matching, comparing labels directly
    would severely underestimate agreement. The Hungarian algorithm finds
    the bijective mapping that maximises the number of concordant
    assignments in O(n^3) time.

    Returns a dict mapping each label in labels_b to its best match
    in labels_a.
    """
    unique_a = sorted(set(labels_a))
    unique_b = sorted(set(labels_b))

    # Build confusion matrix: rows = labels_a, cols = labels_b
    # Cell (i,j) = number of samples with label unique_a[i] in A
    # and unique_b[j] in B
    n_a = len(unique_a)
    n_b = len(unique_b)
    confusion = np.zeros((n_a, n_b), dtype=int)

    a_idx = {v: i for i, v in enumerate(unique_a)}
    b_idx = {v: i for i, v in enumerate(unique_b)}

    for la, lb in zip(labels_a, labels_b):
        confusion[a_idx[la], b_idx[lb]] += 1

    # Hungarian algorithm minimises cost — negate to maximise overlap
    row_ind, col_ind = linear_sum_assignment(-confusion)

    # Build mapping: label_b → label_a
    mapping = {}
    for r, c in zip(row_ind, col_ind):
        mapping[unique_b[c]] = unique_a[r]

    # Any labels in B not covered by the assignment (when |B| > |A|)
    # get mapped to a fallback
    for lb in unique_b:
        if lb not in mapping:
            mapping[lb] = lb  # keep as-is

    return mapping


##########################################################################
# MAIN
##########################################################################

def main():
    args = parse_args()
    out_path = Path(args.out)
    out_dir  = out_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP] crisp_crossval.py — {args.project}")
    print(f"[CRISP] PCA assignments: {args.pca_assign}")
    print(f"[CRISP] MDS assignments: {args.mds_assign}")

    # ── Read assignment files ─────────────────────────────────────────
    pca = pd.read_csv(args.pca_assign, sep="\t")
    mds = pd.read_csv(args.mds_assign, sep="\t")

    # Standardise column names — both files have FID IID GROUP MEMBERSHIP_P OUTLIER
    for df, name in [(pca, "PCA"), (mds, "MDS")]:
        missing = [c for c in ["FID", "IID", "GROUP"] if c not in df.columns]
        if missing:
            print(f"[CRISP] ERROR: {name} assignment file missing columns: {missing}",
                  file=sys.stderr)
            sys.exit(1)

    print(f"[CRISP] PCA samples: {len(pca)}")
    print(f"[CRISP] MDS samples: {len(mds)}")

    # ── Merge on IID ──────────────────────────────────────────────────
    # WHY: We merge on IID (not FID+IID) because family IDs may differ
    # between the PCA and MDS runs if any sample had FID reassignment.
    # IID is the stable per-sample identifier in PLINK.
    merged = pca[["FID", "IID", "GROUP"]].merge(
        mds[["IID", "GROUP"]].rename(columns={"GROUP": "GROUP_MDS"}),
        on="IID",
        how="inner"
    )

    n_merged = len(merged)
    n_pca_only = len(pca) - n_merged
    n_mds_only = len(mds) - n_merged

    if n_merged == 0:
        print("[CRISP] ERROR: No samples in common between PCA and MDS assignments.",
              file=sys.stderr)
        sys.exit(1)

    if n_pca_only > 0 or n_mds_only > 0:
        print(f"[CRISP] WARNING: {n_pca_only} samples in PCA only, "
              f"{n_mds_only} in MDS only — excluded from cross-validation.")

    print(f"[CRISP] Samples in cross-validation: {n_merged}")

    # ── Match cluster labels via Hungarian algorithm ───────────────────
    pca_labels = merged["GROUP"].values
    mds_labels = merged["GROUP_MDS"].values

    print("[CRISP] Matching cluster labels (Hungarian algorithm)...")
    label_mapping = match_labels(pca_labels, mds_labels)
    print(f"[CRISP] Label mapping (MDS → PCA): {label_mapping}")

    # Apply mapping to MDS labels
    mds_labels_matched = np.array([label_mapping[l] for l in mds_labels])

    # ── Per-sample agreement ───────────────────────────────────────────
    # WHY: MDS_AGREE is YES if the Hungarian-matched MDS cluster label
    # equals the PCA cluster label for that sample. This is written to
    # the satellite file so analysts can filter on agreement per sample.
    agree = pca_labels == mds_labels_matched
    merged["MDS_AGREE"] = ["YES" if a else "NO" for a in agree]

    n_agree    = agree.sum()
    n_disagree = (~agree).sum()
    agree_pct  = round(float(n_agree / n_merged * 100), 2)

    print(f"[CRISP] Agreement: {n_agree}/{n_merged} ({agree_pct}%)")
    print(f"[CRISP] Disagreement: {n_disagree} samples")

    # ── Adjusted Rand Index ───────────────────────────────────────────
    # WHY: Raw % agreement is affected by the number of clusters —
    # with many clusters, random chance gives lower agreement. ARI
    # corrects for this, giving a more interpretable concordance metric
    # that is comparable across cohorts with different numbers of clusters.
    ari = adjusted_rand_score(pca_labels, mds_labels)
    ari_rounded = round(float(ari), 4)
    print(f"[CRISP] Adjusted Rand Index (ARI): {ari_rounded}")

    # Interpret ARI
    if ari >= 0.90:
        ari_interpretation = "EXCELLENT — structures are highly concordant"
    elif ari >= 0.70:
        ari_interpretation = "GOOD — minor labelling differences"
    elif ari >= 0.50:
        ari_interpretation = "MODERATE — some structural differences; inspect plots"
    else:
        ari_interpretation = "POOR — significant disagreement; complex structure likely"

    print(f"[CRISP] ARI interpretation: {ari_interpretation}")

    # Warn if below threshold
    if agree_pct < 80:
        print(f"[CRISP] WARNING: Agreement below 80% ({agree_pct}%). "
              "Population structure may be complex — inspect PCA and MDS plots.",
              file=sys.stderr)

    # ── Disagreement breakdown by cluster ─────────────────────────────
    # WHY: Knowing which clusters disagree helps localise the problem.
    # A single cluster with low agreement (e.g. an admixed group that
    # sits differently in PCA vs MDS space) is different from global
    # disagreement across all clusters.
    cluster_agreement = {}
    for group in sorted(set(pca_labels)):
        mask = pca_labels == group
        n_group = mask.sum()
        n_group_agree = (agree & mask).sum()
        pct = round(float(n_group_agree / n_group * 100), 1) if n_group > 0 else 0.0
        cluster_agreement[group] = {
            "n_samples"  : int(n_group),
            "n_agree"    : int(n_group_agree),
            "agree_pct"  : pct,
        }
        print(f"[CRISP]   {group}: {n_group_agree}/{n_group} agree ({pct}%)")

    # ── Write per-sample agreement file ───────────────────────────────
    out_df = merged[["FID", "IID", "MDS_AGREE"]]
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[CRISP] Per-sample agreement written: {out_path}")

    # ── Write JSON summary ─────────────────────────────────────────────
    summary = {
        "schema_version"      : "1.0.0",
        "project"             : args.project,
        "n_samples"           : n_merged,
        "n_pca_only"          : n_pca_only,
        "n_mds_only"          : n_mds_only,
        "n_agree"             : int(n_agree),
        "n_disagree"          : int(n_disagree),
        "agreement_rate_pct"  : agree_pct,
        "adjusted_rand_index" : ari_rounded,
        "ari_interpretation"  : ari_interpretation,
        "label_mapping"       : label_mapping,
        "cluster_agreement"   : cluster_agreement,
        "warning_low_agreement": agree_pct < 80,
    }

    json_path = out_path.parent / f"{args.project}_crossval_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[CRISP] Cross-val summary JSON: {json_path}")
    print(f"[CRISP] crisp_crossval.py complete.")


if __name__ == "__main__":
    main()
