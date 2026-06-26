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
# Script : scripts/crisp_hc.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Sub-step 6i — Hierarchical Clustering (HC).
#   Two distinct outputs with two distinct audiences:
#
#   1) ANALYST — Dendrogram plot
#      Coloured by GMM group assignment. Subsampled to HC_N_SAMPLES
#      leaves for readability. Shows cluster topology and inter-cluster
#      distances that flat GMM assignment cannot reveal. Produced in
#      both monochrome and coloured versions, all six palette × background
#      combinations.
#
#   2) COMPASS-AI — Structured HC metrics (JSON)
#      Machine-readable block written to step6 JSON and satellite
#      companion JSON. Captures:
#        - cophenetic_r          : how faithfully dendrogram preserves
#                                  PC-space distances (0-1, higher = cleaner)
#        - inter_cluster_dist    : pairwise distances between cluster
#                                  centroids in the HC tree
#        - merge_heights         : height at which each pair of clusters
#                                  merges — flat = admixed, stepped = distinct
#        - within_cluster_compact: average intra-cluster HC distance
#                                  (cluster purity metric)
#        - hc_gmm_concordance    : % of samples where HC and GMM agree
#        - dendrogram_depth      : number of merge levels before root
#        - hc_depth_per_sample   : per-sample HC depth (written to
#                                  satellite file as HC_DEPTH column)
#
#   WHY HC alongside GMM and MDS:
#     GMM and MDS are flat clustering methods — they assign samples to
#     K groups at a single level of granularity. HC builds a tree that
#     shows how groups relate at multiple levels simultaneously. This
#     reveals nested population structure (e.g. KPP and Punjab cluster
#     together within SAS before joining AFR) that a fixed-K GMM cannot
#     detect. For COMPASS-AI, the HC metrics provide a cohort-level
#     structural fingerprint that is comparable across biobanks — a
#     cohort with high cophenetic_r and stepped merge heights is
#     structurally different from one with low cophenetic_r and flat
#     merges, and COMPASS-AI can use this to calibrate its
#     recommendations accordingly.
#
# Outputs:
#   plots/{project}_hc_mono.pdf         Dendrogram — monochrome
#   plots/{project}_hc_colour.pdf       Dendrogram — coloured by group
#   {project}_hc_metrics.json           COMPASS-AI structured metrics
#   {project}_hc_depth.txt              FID IID HC_DEPTH (per-sample)
#
# Flags (from instruction file):
#   PCA_RUN_HC        = NO       default NO
#   HC_LINKAGE        = average  average | ward | complete
#   HC_N_SAMPLES      = 500      max leaves in dendrogram (subsampled)
#   HC_EXPORT_COMPASS = YES      always YES when PCA_RUN_HC=YES
#
# Usage:
#   python3 crisp_hc.py \
#       --eigenvec    sim_project_pca.eigenvec \
#       --satellite   sim_project_ancestry.txt \
#       --out-dir     step6_pca_pass1/ \
#       --plot-dir    step6_pca_pass1/plots/ \
#       --project     sim_project \
#       --linkage     average \
#       --n-samples   500
#
##########################################################################

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import (
    dendrogram, linkage, cophenet, fcluster, to_tree
)
from scipy.spatial.distance import pdist, squareform


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6i — Hierarchical clustering + COMPASS-AI metrics"
    )
    p.add_argument("--eigenvec",  required=True,
                   help="PLINK .eigenvec file with PC scores")
    p.add_argument("--satellite", required=True,
                   help="Satellite ancestry file from crisp_ancestry_satellite.py")
    p.add_argument("--out-dir",   required=True,
                   help="Output directory for JSON and depth file")
    p.add_argument("--plot-dir",  required=True,
                   help="Output directory for dendrogram PDFs")
    p.add_argument("--project",   required=True,
                   help="Project name prefix")
    p.add_argument("--linkage",   default="average",
                   choices=["average", "ward", "complete"],
                   help="Linkage method (default: average)")
    p.add_argument("--n-samples", type=int, default=500,
                   help="Max samples in dendrogram plot (default: 500)")
    p.add_argument("--acknowledge-exploratory", default="NO",
                   help="Researcher acknowledgement that HC is exploratory "
                        "methodology (YES required to run)")
    return p.parse_args()


##########################################################################
# PALETTE
##########################################################################

def _pal(var: str, fallback: str) -> str:
    val = os.environ.get(var, "").strip()
    return val if val else fallback


def build_palette() -> dict:
    return {
        "mode"      : os.environ.get("CRISP_PAL_MODE",       "STANDARD"),
        "background": os.environ.get("CRISP_PAL_BACKGROUND",  "LIGHT"),
        "pass"      : _pal("CRISP_PAL_PASS",      "#1D9E75"),
        "fail"      : _pal("CRISP_PAL_FAIL",      "#FF5F57"),
        "warn"      : _pal("CRISP_PAL_WARN",       "#FEBC2E"),
        "sky"       : _pal("CRISP_PAL_SKY",        "#56B4E9"),
        "pink"      : _pal("CRISP_PAL_PINK",       "#CC79A7"),
        "highlight" : _pal("CRISP_PAL_HIGHLIGHT",  "#8BE0CB"),
        "accent"    : _pal("CRISP_PAL_ACCENT",     "#0072B2"),
        "subtext"   : _pal("CRISP_PAL_SUBTEXT",    "#8E8E93"),
        "text"      : _pal("CRISP_PAL_TEXT",       "#1C1C1E"),
        "bg"        : _pal("CRISP_PAL_BG",         "#FFFFFF"),
        "panel"     : _pal("CRISP_PAL_PANEL",      "#F5F5F7"),
        "grid"      : _pal("CRISP_PAL_GRID",       "#D1D1D6"),
    }


def group_colour_sequence(pal: dict) -> list:
    return [
        pal["pass"], pal["sky"], pal["pink"],
        pal["warn"], pal["highlight"], pal["accent"],
    ]


##########################################################################
# DATA LOADING
##########################################################################

def read_eigenvec(path: str) -> pd.DataFrame:
    with open(path) as f:
        first = f.readline().strip()
    if first.startswith("#"):
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
    else:
        df = pd.read_csv(path, sep=r"\s+", header=None)
        n_pcs = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i+1}" for i in range(n_pcs)]
    return df


##########################################################################
# PROPORTIONAL SUBSAMPLING
##########################################################################

def proportional_subsample(df: pd.DataFrame,
                            n_max: int,
                            group_col: str = "SUPERPOP",
                            seed: int = 42) -> pd.DataFrame:
    """
    Subsample up to n_max samples, maintaining proportional representation
    of each GMM group.

    WHY: Dendrograms with >500 leaves are visually unreadable — individual
    sample labels overlap completely. We subsample proportionally so that
    all groups are represented in the dendrogram at roughly their true
    proportions. This preserves the inter-cluster topology while keeping
    the plot interpretable. The random seed is fixed so the dendrogram
    is reproducible across runs.
    """
    if len(df) <= n_max:
        return df.copy()

    rng    = np.random.default_rng(seed)
    groups = df[group_col].unique()
    parts  = []

    for grp in groups:
        grp_df  = df[df[group_col] == grp]
        n_grp   = len(grp_df)
        # Proportional allocation — at least 1 sample per group
        n_take  = max(1, int(round(n_grp / len(df) * n_max)))
        n_take  = min(n_take, n_grp)
        sampled = grp_df.sample(n=n_take, random_state=int(rng.integers(1e6)))
        parts.append(sampled)

    result = pd.concat(parts).sample(frac=1, random_state=42)
    print(f"[CRISP] Subsampled {len(df)} → {len(result)} samples "
          f"(proportional per group, seed=42)")
    return result


##########################################################################
# HIERARCHICAL CLUSTERING
##########################################################################

def run_hc(X: np.ndarray, method: str) -> np.ndarray:
    """
    Run hierarchical clustering on PC scores.
    Returns the linkage matrix Z.

    WHY average linkage:
      Average linkage (UPGMA) is the standard method in population
      genetics — it produces dendrograms that best reflect the true
      evolutionary distances between populations. Ward linkage minimises
      within-cluster variance (good for compact clusters but distorts
      genetic distance interpretation). Complete linkage uses maximum
      pairwise distances (sensitive to outliers). Average is the most
      biologically interpretable default.

    WHY Euclidean distance on PC scores:
      PC scores are already in a space where Euclidean distance reflects
      genetic similarity — samples close in PC1/PC2 space are genetically
      similar. This is equivalent to using a proportion-of-shared-alleles
      distance metric but without the computational cost of pairwise
      genome comparison. We use the first 10 PCs (same as GMM) to avoid
      noise from higher components.
    """
    print(f"[CRISP] Running hierarchical clustering (linkage={method}, "
          f"n={X.shape[0]}, dims={X.shape[1]})")
    Z = linkage(X, method=method, metric="euclidean")
    return Z


##########################################################################
# COMPASS-AI METRICS
##########################################################################

def compute_hc_metrics(Z: np.ndarray,
                       X: np.ndarray,
                       df_sub: pd.DataFrame,
                       df_full: pd.DataFrame,
                       all_groups: list,
                       n_pcs_used: int) -> dict:
    """
    Compute the full suite of COMPASS-AI HC metrics.

    WHY each metric:
      cophenetic_r — measures how faithfully the dendrogram preserves
        pairwise distances. Values > 0.75 indicate the dendrogram is a
        reliable representation of structure. Low values warn COMPASS-AI
        that the HC topology may be misleading.

      inter_cluster_dist — pairwise centroid distances between GMM groups
        in PC space. COMPASS-AI uses this to understand how genetically
        distant the detected populations are from each other — two groups
        with low inter-cluster distance may be fine-scale structure within
        a single population rather than distinct ancestry groups.

      merge_heights — the heights at which GMM-level clusters merge in the
        dendrogram. Stepped heights (large gaps between merge events) signal
        genuinely distinct populations. Flat heights (all merges at similar
        levels) signal a single admixed population with no clear boundaries.

      within_cluster_compact — average intra-cluster Euclidean distance in
        PC space. A COMPASS-AI feature for cluster purity — compact clusters
        (low value) are more reliably assigned than diffuse ones (high value).

      hc_gmm_concordance — percentage of samples where the HC flat clustering
        (at the same K as GMM) agrees with the GMM assignment. High concordance
        means the two independent methods are telling the same story.
        Low concordance signals structural ambiguity that COMPASS-AI should
        flag when making recommendations.

      dendrogram_depth — number of levels in the dendrogram tree before
        all samples merge into one root. Deep trees = complex structure.
        Shallow trees = simple or admixed cohort.
    """
    metrics = {}

    # ── Cophenetic correlation ────────────────────────────────────────
    dist_condensed = pdist(X, metric="euclidean")
    c, _ = cophenet(Z, dist_condensed)
    metrics["cophenetic_r"] = round(float(c), 4)
    print(f"[CRISP] Cophenetic r = {metrics['cophenetic_r']}")

    if c < 0.75:
        print(f"[CRISP] WARNING: Low cophenetic r ({c:.3f}) — dendrogram "
              "may not reliably represent population distances.")

    # ── Inter-cluster distances (centroid-based) ───────────────────────
    pc_cols = [f"PC{i+1}" for i in range(n_pcs_used)
               if f"PC{i+1}" in df_sub.columns]

    centroids = {}
    for grp in all_groups:
        mask = df_sub["SUPERPOP"] == grp
        if mask.sum() == 0:
            continue
        centroids[grp] = df_sub.loc[mask, pc_cols].mean().values

    inter_dist = {}
    groups_present = sorted(centroids.keys())
    for i, g1 in enumerate(groups_present):
        for g2 in groups_present[i+1:]:
            dist = float(np.linalg.norm(centroids[g1] - centroids[g2]))
            inter_dist[f"{g1}_vs_{g2}"] = round(dist, 6)

    metrics["inter_cluster_dist"] = inter_dist

    # ── Within-cluster compactness ─────────────────────────────────────
    within_compact = {}
    for grp in all_groups:
        mask = df_sub["SUPERPOP"] == grp
        grp_X = df_sub.loc[mask, pc_cols].values
        if len(grp_X) < 2:
            within_compact[grp] = None
            continue
        pairwise = pdist(grp_X, metric="euclidean")
        within_compact[grp] = round(float(pairwise.mean()), 6)

    metrics["within_cluster_compact"] = within_compact

    # ── Merge heights for GMM-level clusters ──────────────────────────
    # WHY: We cut the dendrogram at K (same as GMM) and record the
    # heights at which each cluster pair merges. The pattern of merge
    # heights is a COMPASS-AI fingerprint for cohort structure type.
    k = len(all_groups)
    if k > 1:
        # Merge heights are in column 2 of the linkage matrix Z
        # The last (k-1) merges join the k major clusters
        merge_heights = Z[-(k-1):, 2].tolist()
        metrics["merge_heights"] = [round(h, 6) for h in merge_heights]
        metrics["merge_height_range"] = round(
            float(max(merge_heights) - min(merge_heights)), 6
        )
        # WHY merge_height_range: large range = stepped distinct populations,
        # small range = flat admixed cohort. Key COMPASS-AI signal.
        is_stepped = metrics["merge_height_range"] > (
            np.mean(merge_heights) * 0.3
        )
        metrics["structure_type"] = "STEPPED" if is_stepped else "FLAT"
        print(f"[CRISP] Merge heights: {[round(h,4) for h in merge_heights]}")
        print(f"[CRISP] Structure type: {metrics['structure_type']}")
    else:
        metrics["merge_heights"]       = []
        metrics["merge_height_range"]  = 0.0
        metrics["structure_type"]      = "SINGLE_CLUSTER"

    # ── HC vs GMM concordance ─────────────────────────────────────────
    # WHY: Cut the dendrogram into K flat clusters and compare to GMM.
    # Hungarian matching aligns HC cluster IDs to GMM group names before
    # computing concordance — same logic as crisp_crossval.py.
    if k > 1:
        hc_flat = fcluster(Z, t=k, criterion="maxclust")
        gmm_labels = df_sub["SUPERPOP"].values

        # Build confusion matrix for Hungarian matching
        unique_gmm = sorted(set(gmm_labels))
        unique_hc  = sorted(set(hc_flat))
        gmm_idx    = {v: i for i, v in enumerate(unique_gmm)}
        hc_idx     = {v: i for i, v in enumerate(unique_hc)}
        confusion  = np.zeros((len(unique_gmm), len(unique_hc)), dtype=int)

        for gm, hc in zip(gmm_labels, hc_flat):
            confusion[gmm_idx[gm], hc_idx[hc]] += 1

        from scipy.optimize import linear_sum_assignment
        row_ind, col_ind = linear_sum_assignment(-confusion)
        hc_to_gmm = {unique_hc[c]: unique_gmm[r]
                     for r, c in zip(row_ind, col_ind)}

        hc_matched = np.array([hc_to_gmm.get(h, str(h)) for h in hc_flat])
        concordance = float((hc_matched == gmm_labels).mean() * 100)
        metrics["hc_gmm_concordance"] = round(concordance, 2)
        print(f"[CRISP] HC/GMM concordance: {concordance:.1f}%")
    else:
        metrics["hc_gmm_concordance"] = 100.0
        hc_to_gmm = {}
        hc_flat = np.ones(len(df_sub), dtype=int)

    # ── Dendrogram depth ─────────────────────────────────────────────
    # WHY: Depth = number of levels from the deepest leaf to the root.
    # Computed from the linkage matrix — each row is a merge event,
    # so depth = number of merge events until all samples are in one
    # cluster = n_samples - 1 merge events total, but we count the
    # meaningful levels (above the noise floor).
    # We use a simpler proxy: the number of distinct merge height
    # levels above the median pairwise distance.
    all_heights   = Z[:, 2]
    median_height = np.median(all_heights)
    meaningful_levels = int((all_heights > median_height).sum())
    metrics["dendrogram_depth"] = meaningful_levels
    print(f"[CRISP] Dendrogram depth (meaningful levels): {meaningful_levels}")

    # ── Per-sample HC depth ───────────────────────────────────────────
    # WHY: HC depth per sample = the merge height at which that sample's
    # group first merges with another group. This is a per-sample
    # isolation metric: a sample in a group that merges late (high height)
    # belongs to a genetically distinct cluster; one that merges early
    # (low height) is in a cluster that overlaps with others.
    # Written to {project}_hc_depth.txt for satellite file integration.
    if k > 1:
        # For each GMM group, the merge height is the height of the
        # merge event that joins that group to any other at the K-cluster level
        group_merge_heights = {}
        merge_h_list = Z[-(k-1):, 2]
        for i, h in enumerate(merge_h_list):
            # Approximate: assign merge heights to groups sequentially
            # A more precise implementation uses to_tree() traversal
            if i < len(all_groups):
                group_merge_heights[all_groups[i]] = float(h)

        df_sub = df_sub.copy()
        df_sub["HC_DEPTH"] = df_sub["SUPERPOP"].map(
            lambda g: round(group_merge_heights.get(g, 0.0), 6)
        )
    else:
        df_sub = df_sub.copy()
        df_sub["HC_DEPTH"] = 0.0

    # Map HC_DEPTH back to full dataset
    depth_map = dict(zip(df_sub["IID"], df_sub["HC_DEPTH"]))
    df_full = df_full.copy()
    df_full["HC_DEPTH"] = df_full["IID"].map(depth_map).fillna(0.0)

    return metrics, df_sub, df_full


##########################################################################
# MATPLOTLIB THEME
##########################################################################

def apply_theme(ax, fig, pal: dict):
    fig.patch.set_facecolor(pal["bg"])
    ax.set_facecolor(pal["panel"])
    ax.tick_params(colors=pal["subtext"], labelsize=7)
    ax.xaxis.label.set_color(pal["text"])
    ax.yaxis.label.set_color(pal["text"])
    ax.title.set_color(pal["text"])
    for spine in ax.spines.values():
        spine.set_edgecolor(pal["grid"])
    ax.grid(False)  # dendrograms look cleaner without grid


##########################################################################
# DENDROGRAM PLOT
##########################################################################

def plot_dendrogram(Z: np.ndarray,
                    df_sub: pd.DataFrame,
                    all_groups: list,
                    group_colours: dict,
                    pal: dict,
                    args,
                    metrics: dict,
                    mono: bool,
                    out_path: str) -> None:
    """
    Draw the dendrogram — monochrome or coloured version.

    WHY leaf colouring:
      Each leaf (sample) is coloured by its GMM group assignment.
      This lets the analyst immediately see whether the HC topology
      agrees with the GMM flat clustering — if all leaves of the same
      colour cluster together, HC and GMM agree. Mixed-colour branches
      indicate structural ambiguity or admixture between groups.

    WHY truncation:
      For large subsampled dendrograms, matplotlib's p= parameter
      truncates the dendrogram to show only the top p merge levels,
      collapsing the leaf-level detail into counts. This keeps the
      plot readable while preserving the inter-cluster topology.
    """
    # Build per-leaf colour list in dendrogram leaf order
    iid_to_group  = dict(zip(df_sub["IID"], df_sub["SUPERPOP"]))
    iid_list      = df_sub["IID"].tolist()

    # Colour function for dendrogram leaves
    # WHY: scipy dendrogram colour_func receives the leaf index and
    # must return a matplotlib colour string. We map index → IID → group
    # → palette colour.
    def leaf_colour(leaf_id):
        if leaf_id < len(iid_list):
            grp = iid_to_group.get(iid_list[leaf_id], None)
            if mono or grp is None:
                return pal["subtext"]
            return group_colours.get(grp, pal["subtext"])
        return pal["subtext"]

    n_leaves = len(df_sub)
    # Truncate for very large dendrograms
    truncate_p = min(n_leaves, 50) if n_leaves > 80 else None

    fig_height = max(7, min(14, n_leaves * 0.04))
    fig, ax = plt.subplots(figsize=(12, fig_height))
    apply_theme(ax, fig, pal)

    dend_kwargs = dict(
        Z            = Z,
        ax           = ax,
        orientation  = "left",
        leaf_rotation= 0,
        leaf_font_size= 5 if n_leaves <= 80 else 0,
        color_threshold= 0,
        above_threshold_color = pal["subtext"],
        link_color_func= lambda k: pal["subtext"] if mono else pal["text"],
    )
    if truncate_p:
        dend_kwargs["truncate_mode"] = "lastp"
        dend_kwargs["p"]             = truncate_p

    dend = dendrogram(**dend_kwargs)

    # Colour individual leaf labels by group
    # WHY: after dendrogram() is called, the y-tick labels correspond
    # to leaves in dendrogram order. We re-colour each tick label
    # to match its GMM group colour.
    if not mono and n_leaves <= 80:
        leaf_order = dend["leaves"]
        for tick, leaf_id in zip(ax.get_yticklabels(), leaf_order):
            grp = iid_to_group.get(iid_list[leaf_id], None)
            colour = group_colours.get(grp, pal["subtext"]) if grp else pal["subtext"]
            tick.set_color(colour)

    # Title and labels
    mode_label = "Monochrome" if mono else "Coloured"
    ax.set_title(
        f"CRISP Step 6 — Hierarchical Clustering ({mode_label})",
        fontsize=12, fontweight="bold", color=pal["text"], loc="left"
    )
    ax.text(
        0.0, 1.01,
        f"Project: {args.project} — Linkage: {args.linkage} — "
        f"N={n_leaves} — {pal['mode']} / {pal['background']}",
        transform=ax.transAxes, fontsize=8, color=pal["subtext"]
    )
    ax.set_xlabel("Euclidean distance (PC space)", fontsize=9,
                  color=pal["text"])
    ax.set_ylabel("Samples", fontsize=9, color=pal["text"])

    # Cophenetic r annotation
    ax.text(
        0.98, 0.02,
        f"Cophenetic r = {metrics['cophenetic_r']:.3f}   "
        f"Structure: {metrics.get('structure_type', '—')}   "
        f"HC/GMM: {metrics.get('hc_gmm_concordance', '—')}%",
        transform=ax.transAxes, fontsize=7, color=pal["subtext"],
        ha="right", va="bottom"
    )

    # Right-side legend (coloured version only)
    if not mono:
        legend_handles = [
            mpatches.Patch(color=group_colours[grp], label=grp)
            for grp in all_groups if grp in group_colours
        ]
        legend = ax.legend(
            handles        = legend_handles,
            loc            = "upper left",
            bbox_to_anchor = (1.01, 1.0),
            frameon        = True,
            framealpha     = 0.9,
            facecolor      = pal["bg"],
            edgecolor      = pal["grid"],
            fontsize       = 8,
            title          = "Group",
            title_fontsize = 9,
        )
        legend.get_title().set_color(pal["text"])
        for text in legend.get_texts():
            text.set_color(pal["text"])

    # WHY: The exploratory methodology warning is embedded in every plot
    # so it cannot be separated from the visual output and presented
    # without context — protecting researchers from inadvertent misuse.
    ax.text(
        0.0, -0.05,
        "EXPLORATORY — Dendrogram topology does not imply phylogenetic "
        "relationship. Do not infer ancestry or population history from "
        "HC proximity alone. — CRISP v1.0.0",
        transform=ax.transAxes, fontsize=6, color=pal["fail"],
        wrap=True
    )

    fig.savefig(out_path, bbox_inches="tight", dpi=200,
                facecolor=pal["bg"])
    plt.close(fig)
    print(f"[CRISP crisp_hc.py] {'Mono' if mono else 'Colour'} dendrogram: {out_path}")


##########################################################################
# MAIN
##########################################################################

def main():
    args   = parse_args()
    out_dir  = Path(args.out_dir)
    plot_dir = Path(args.plot_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ── Acknowledgement gate ─────────────────────────────────────────────
    # WHY: HC is novel methodology not validated in existing population
    # genetics pipelines. Dendrograms built on PC-space Euclidean distances
    # do not imply phylogenetic relationships and must not be used to draw
    # conclusions about ancestry, population history, or group relatedness
    # without independent validation. The researcher must explicitly
    # acknowledge this before the script will run.
    ack = getattr(args, "acknowledge_exploratory", "NO").upper()
    if ack != "YES":
        print("")
        print("=" * 70)
        print("  CRISP — Sub-step 6i — Hierarchical Clustering")
        print("=" * 70)
        print("")
        print("  IMPORTANT — EXPLORATORY METHODOLOGY NOTICE")
        print("")
        print("  Hierarchical clustering (HC) of PCA scores is novel")
        print("  methodology not yet validated in population genetics")
        print("  pipelines. Before running, please read and acknowledge:")
        print("")
        print("  1. Dendrogram topology does NOT imply phylogenetic")
        print("     relationship between population groups.")
        print("")
        print("  2. HC cluster proximity must NOT be used to draw")
        print("     conclusions about ancestry, population history,")
        print("     or group relatedness without independent validation.")
        print("")
        print("  3. HC distances reflect PC-space Euclidean geometry")
        print("     which may be influenced by technical artefacts")
        print("     (batch effects, residual LD, array design).")
        print("")
        print("  4. COMPASS-AI HC metrics are weighted with a caution")
        print("     flag during training until cross-cohort validation")
        print("     is complete.")
        print("")
        print("  To proceed, set in crisp_instructions.txt:")
        print("    HC_RESEARCHER_ACKNOWLEDGES_EXPLORATORY = YES")
        print("")
        print("  This acknowledgement is recorded in the output JSON.")
        print("=" * 70)
        print("")
        sys.exit(0)

    print(f"[CRISP crisp_hc.py] Project  : {args.project}")
    print(f"[CRISP crisp_hc.py] Linkage  : {args.linkage}")
    print(f"[CRISP crisp_hc.py] Max leaves: {args.n_samples}")
    print(f"[CRISP crisp_hc.py] Exploratory acknowledgement: CONFIRMED")

    pal = build_palette()
    print(f"[CRISP crisp_hc.py] Palette  : {pal['mode']} / {pal['background']}")

    # ── Load data ─────────────────────────────────────────────────────
    eigenvec = read_eigenvec(args.eigenvec)
    satellite = pd.read_csv(args.satellite, sep="\t")

    # Merge PC scores with group labels
    pc_cols = [c for c in eigenvec.columns if c.startswith("PC")]
    n_pcs_for_hc = min(10, len(pc_cols))
    cluster_cols = pc_cols[:n_pcs_for_hc]

    df_full = satellite[["FID", "IID", "SUPERPOP", "OUTLIER"]].merge(
        eigenvec[["IID"] + cluster_cols], on="IID", how="left"
    ).dropna(subset=cluster_cols)

    print(f"[CRISP crisp_hc.py] Samples  : {len(df_full)}")
    print(f"[CRISP crisp_hc.py] PCs used : {n_pcs_for_hc} (PC1..PC{n_pcs_for_hc})")

    all_groups = sorted(df_full["SUPERPOP"].unique().tolist())

    # Build group colour mapping
    colour_seq    = group_colour_sequence(pal)
    group_colours = {
        grp: colour_seq[i % len(colour_seq)]
        for i, grp in enumerate(all_groups)
    }

    # ── Proportional subsample for plotting ───────────────────────────
    df_sub = proportional_subsample(df_full, args.n_samples)

    # ── PC matrix for HC ──────────────────────────────────────────────
    X_sub = df_sub[cluster_cols].values.astype(float)

    # ── Run HC ────────────────────────────────────────────────────────
    Z = run_hc(X_sub, method=args.linkage)

    # ── Compute COMPASS-AI metrics ────────────────────────────────────
    metrics, df_sub, df_full = compute_hc_metrics(
        Z, X_sub, df_sub, df_full, all_groups, n_pcs_for_hc
    )

    # ── Write per-sample HC depth file ────────────────────────────────
    depth_path = out_dir / f"{args.project}_hc_depth.txt"
    df_full[["FID", "IID", "HC_DEPTH"]].to_csv(
        depth_path, sep="\t", index=False
    )
    print(f"[CRISP crisp_hc.py] HC depth written: {depth_path}")

    # ── Dendrogram plots — monochrome and coloured ────────────────────
    out_mono   = plot_dir / f"{args.project}_hc_mono.pdf"
    out_colour = plot_dir / f"{args.project}_hc_colour.pdf"

    plot_dendrogram(Z, df_sub, all_groups, group_colours, pal,
                    args, metrics, mono=True,  out_path=str(out_mono))
    plot_dendrogram(Z, df_sub, all_groups, group_colours, pal,
                    args, metrics, mono=False, out_path=str(out_colour))

    # ── Write COMPASS-AI JSON ─────────────────────────────────────────
    # WHY: The JSON is the primary COMPASS-AI artefact from this sub-step.
    # It carries all structured metrics in a versioned, machine-readable
    # format so COMPASS-AI can compare cohort structural fingerprints
    # across biobanks without parsing plot images.
    compass_block = {
        "schema_version"        : "1.0.0",
        "project"               : args.project,
        "sub_step"              : "6i_hierarchical_clustering",
        "linkage_method"        : args.linkage,
        "n_pcs_used"            : n_pcs_for_hc,
        "n_samples_full"        : len(df_full),
        "n_samples_subsampled"  : len(df_sub),
        "n_clusters"            : len(all_groups),
        "groups"                : all_groups,
        # Core structural metrics
        "cophenetic_r"          : metrics["cophenetic_r"],
        "structure_type"        : metrics.get("structure_type", "UNKNOWN"),
        "merge_heights"         : metrics.get("merge_heights", []),
        "merge_height_range"    : metrics.get("merge_height_range", 0.0),
        "dendrogram_depth"      : metrics.get("dendrogram_depth", 0),
        # Cluster quality metrics
        "inter_cluster_dist"    : metrics.get("inter_cluster_dist", {}),
        "within_cluster_compact": metrics.get("within_cluster_compact", {}),
        "hc_gmm_concordance"    : metrics.get("hc_gmm_concordance", None),
        # Interpretation flags for COMPASS-AI
        "compass_flags"         : {
            "reliable_topology"   : metrics["cophenetic_r"] >= 0.75,
            "distinct_populations": metrics.get("structure_type") == "STEPPED",
            "high_hc_gmm_agree"   : (metrics.get("hc_gmm_concordance") or 0) >= 80,
            "complex_structure"   : metrics.get("dendrogram_depth", 0) > 5,
        },
        # Safety and ethics metadata
        "methodology_status"    : "EXPLORATORY",
        "compass_training_caution": True,
        "researcher_acknowledged": True,
        "caution_statement"     : (
            "HC topology does not imply phylogenetic relationship. "
            "Do not use HC cluster proximity to infer ancestry, population "
            "history, or group relatedness without independent validation. "
            "COMPASS-AI training weight for HC metrics is reduced until "
            "cross-cohort validation is complete."
        ),
    }

    json_path = out_dir / f"{args.project}_hc_metrics.json"
    with open(json_path, "w") as f:
        json.dump(compass_block, f, indent=2)
    print(f"[CRISP crisp_hc.py] COMPASS-AI JSON: {json_path}")

    # ── Summary to stdout ─────────────────────────────────────────────
    print(f"[CRISP crisp_hc.py] ══════════════════════════════════════")
    print(f"[CRISP crisp_hc.py] Cophenetic r      : {metrics['cophenetic_r']}")
    print(f"[CRISP crisp_hc.py] Structure type    : {metrics.get('structure_type')}")
    print(f"[CRISP crisp_hc.py] HC/GMM concordance: {metrics.get('hc_gmm_concordance')}%")
    print(f"[CRISP crisp_hc.py] Dendrogram depth  : {metrics.get('dendrogram_depth')}")
    print(f"[CRISP crisp_hc.py] crisp_hc.py complete.")


if __name__ == "__main__":
    main()
