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
# Script : scripts/plot_pca_pass1.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Python/matplotlib parity of plot_pca_pass1.R.
#   Produces four plots for CRISP Step 6 — PCA Pass 1:
#
#   Plot 1 — Monochrome PCA scatter + ellipses
#   Plot 2 — Coloured PCA scatter + ellipses
#   Plot 3 — Monochrome MDS scatter + ellipses  (if MDS was run)
#   Plot 4 — Coloured MDS scatter + ellipses    (if MDS was run)
#   Plot 5 — Scree plot (PCA only)
#
#   All six palette × background combinations are fully supported.
#   Colours consumed exclusively from CRISP_PAL_* environment
#   variables — zero hardcoded hex values in this script.
#
#   Right-side external legend — non-negotiable CRISP convention.
#   Outliers always CRISP_PAL_FAIL, shape = × (marker='x').
#
# Args:
#   --input       path to .eigenvec or .mds file
#   --satellite   path to {project}_ancestry.txt satellite file
#   --plot-type   PCA or MDS
#   --out-dir     output directory for PDFs
#   --project     project name prefix
#   --n-clusters  number of detected clusters (integer)
#   --ellipse-sd  SD radius for confidence ellipse (default 3)
#   --admixed     COHORT_ADMIXED flag (YES/NO)
#
##########################################################################

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # non-interactive backend for server/pipeline use
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.stats import chi2


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6 — PCA Pass 1 matplotlib plots"
    )
    p.add_argument("--input",      required=True)
    p.add_argument("--satellite",  required=True)
    p.add_argument("--plot-type",  default="PCA", choices=["PCA", "MDS"])
    p.add_argument("--out-dir",    required=True)
    p.add_argument("--project",    required=True)
    p.add_argument("--n-clusters", type=int, default=None)
    p.add_argument("--ellipse-sd", type=float, default=3.0)
    p.add_argument("--admixed",    default="NO")
    return p.parse_args()


##########################################################################
# PALETTE
##########################################################################

def _pal(var: str, fallback: str) -> str:
    """
    Read a colour from a CRISP_PAL_* environment variable.
    Returns fallback if the variable is unset or empty.
    WHY: Zero hardcoded hex values in this script. All colours flow
    from the environment set by _crisp_palette() in crisp_pca_pass1.sh.
    Fallbacks are STANDARD/LIGHT values so the script degrades
    gracefully in manual/test runs without the shell wrapper.
    """
    val = os.environ.get(var, "").strip()
    return val if val else fallback


def build_palette() -> dict:
    """
    Build the full colour dictionary from CRISP_PAL_* env vars.
    WHY: Centralising palette construction means every plot function
    receives the same colour object — no risk of individual functions
    reading env vars inconsistently.
    """
    return {
        "mode"       : os.environ.get("CRISP_PAL_MODE",       "STANDARD"),
        "background" : os.environ.get("CRISP_PAL_BACKGROUND",  "LIGHT"),
        "pass"       : _pal("CRISP_PAL_PASS",      "#1D9E75"),
        "fail"       : _pal("CRISP_PAL_FAIL",      "#FF5F57"),
        "warn"       : _pal("CRISP_PAL_WARN",       "#FEBC2E"),
        "sky"        : _pal("CRISP_PAL_SKY",        "#56B4E9"),
        "pink"       : _pal("CRISP_PAL_PINK",       "#CC79A7"),
        "highlight"  : _pal("CRISP_PAL_HIGHLIGHT",  "#8BE0CB"),
        "accent"     : _pal("CRISP_PAL_ACCENT",     "#0072B2"),
        "subtext"    : _pal("CRISP_PAL_SUBTEXT",    "#8E8E93"),
        "text"       : _pal("CRISP_PAL_TEXT",       "#1C1C1E"),
        "bg"         : _pal("CRISP_PAL_BG",         "#FFFFFF"),
        "panel"      : _pal("CRISP_PAL_PANEL",      "#F5F5F7"),
        "grid"       : _pal("CRISP_PAL_GRID",       "#D1D1D6"),
    }


def group_colour_sequence(pal: dict) -> list:
    """
    Return the ordered sequence of group colours.
    WHY: Group colours cycle through this sequence — GROUP_1 gets
    pal['pass'], GROUP_2 gets pal['sky'], etc. The sequence is
    defined once so R and Python plots are always consistent.
    """
    return [
        pal["pass"],
        pal["sky"],
        pal["pink"],
        pal["warn"],
        pal["highlight"],
        pal["accent"],
    ]


##########################################################################
# DATA READERS
##########################################################################

def read_eigenvec(path: str) -> pd.DataFrame:
    """
    Read PLINK .eigenvec (PLINK 1.9 or PLINK 2 format).
    """
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


def read_mds(path: str) -> pd.DataFrame:
    """
    Read PLINK .mds file, renaming C1..Cn → PC1..PCn.
    """
    df = pd.read_csv(path, sep=r"\s+")
    if "SOL" in df.columns:
        df = df.drop(columns=["SOL"])
    df.columns = [c.replace("C", "PC", 1) if c.startswith("C") and c[1:].isdigit()
                  else c for c in df.columns]
    return df


def read_satellite(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_plot_data(args) -> tuple[pd.DataFrame, str, str]:
    """
    Load and merge the data needed for plotting.
    Returns (plot_df, x_label, y_label).
    """
    sat = read_satellite(args.satellite)

    if args.plot_type == "MDS":
        mds = read_mds(args.input)
        plot_df = mds.merge(
            sat[["IID", "SUPERPOP", "OUTLIER", "MEMBERSHIP_P"]],
            on="IID", how="left"
        )
        x_label, y_label = "MDS Dimension 1", "MDS Dimension 2"
    else:
        plot_df = sat.copy()
        x_label, y_label = "PC1", "PC2"

    plot_df["SUPERPOP"] = plot_df["SUPERPOP"].astype(str)
    plot_df["OUTLIER"]  = plot_df["OUTLIER"].astype(str)
    plot_df["PC1"]      = pd.to_numeric(plot_df["PC1"], errors="coerce")
    plot_df["PC2"]      = pd.to_numeric(plot_df["PC2"], errors="coerce")
    plot_df = plot_df.dropna(subset=["PC1", "PC2"])

    return plot_df, x_label, y_label


##########################################################################
# MATPLOTLIB THEME
##########################################################################

def apply_theme(ax, fig, pal: dict):
    """
    Apply CRISP theme to a matplotlib axes and figure.
    WHY: All background, grid, text, and spine colours are driven by
    CRISP_PAL_* variables so NIGHT/DARK mode works without code changes.
    """
    fig.patch.set_facecolor(pal["bg"])
    ax.set_facecolor(pal["panel"])

    ax.tick_params(colors=pal["subtext"], labelsize=8)
    ax.xaxis.label.set_color(pal["text"])
    ax.yaxis.label.set_color(pal["text"])
    ax.title.set_color(pal["text"])

    for spine in ax.spines.values():
        spine.set_edgecolor(pal["grid"])

    ax.grid(True, color=pal["grid"], linewidth=0.3, alpha=0.7)
    ax.set_axisbelow(True)


##########################################################################
# CONFIDENCE ELLIPSE
##########################################################################

def confidence_ellipse(x: np.ndarray,
                       y: np.ndarray,
                       ax,
                       sd: float,
                       colour: str,
                       fill: bool = False,
                       fill_alpha: float = 0.08,
                       linestyle: str = "-",
                       linewidth: float = 0.9) -> None:
    """
    Draw a confidence ellipse at sd standard deviations in 2D.

    WHY: We use the eigendecomposition of the covariance matrix to find
    the major and minor axes of the ellipse and their orientations.
    The chi-squared threshold at df=2 converts the SD radius to the
    correct 2D confidence boundary:
        threshold = chi2.ppf(pchisq(sd^2, df=2), df=2) = sd^2
    The ellipse width and height are scaled by sqrt(eigenvalue * threshold).

    This is the matplotlib equivalent of ggplot2's stat_ellipse(type='norm').
    Results are numerically identical to the R version.
    """
    if len(x) < 5:
        return

    # Covariance matrix of the cluster
    cov = np.cov(x, y)
    # WHY: eigh guarantees real eigenvalues for symmetric matrices,
    # unlike eig which may return complex values due to floating point.
    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    # Sort by descending eigenvalue so width >= height
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    # Chi-squared threshold for sd standard deviations in 2D
    # WHY: chi2.ppf(level, df=2) where level = pchisq(sd^2, df=2)
    # = chi2.cdf(sd^2, df=2). This is the same as sd^2 by definition
    # of the chi-squared distribution, so the threshold = sd^2.
    chi2_threshold = sd ** 2

    # Ellipse width and height (2 * semi-axis length)
    width  = 2.0 * np.sqrt(eigenvalues[0] * chi2_threshold)
    height = 2.0 * np.sqrt(eigenvalues[1] * chi2_threshold)

    # Rotation angle of the major axis
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))

    # Centroid
    cx, cy = np.mean(x), np.mean(y)

    ellipse = Ellipse(
        xy        = (cx, cy),
        width     = width,
        height    = height,
        angle     = angle,
        edgecolor = colour,
        facecolor = colour if fill else "none",
        alpha     = fill_alpha if fill else 0.9,
        linewidth = linewidth,
        linestyle = linestyle,
        zorder    = 2,
    )
    ax.add_patch(ellipse)


##########################################################################
# CLUSTER CENTROIDS
##########################################################################

def compute_centroids(plot_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-cluster centroid from non-outlier samples.
    WHY: Centroid labels are placed at the geometric centre of each
    cluster (excluding outliers) so the label sits inside the main
    body of the cluster rather than being dragged by outliers.
    """
    non_out = plot_df[plot_df["OUTLIER"] != "YES"]
    centroids = (
        non_out.groupby("SUPERPOP")[["PC1", "PC2"]]
        .mean()
        .reset_index()
        .rename(columns={"PC1": "cx", "PC2": "cy"})
    )
    return centroids


##########################################################################
# PLOT 1 — MONOCHROME
##########################################################################

def plot_mono(plot_df: pd.DataFrame,
              pal: dict,
              args,
              x_label: str,
              y_label: str,
              all_groups: list,
              out_path: str) -> None:
    """
    Monochrome scatter — all samples in PAL_SUBTEXT, outliers in PAL_FAIL.
    Dashed ellipses in PAL_TEXT. Centroid labels in bold.
    WHY: The analyst sees cluster geometry before colour is applied.
    """
    fig, ax = plt.subplots(figsize=(9, 7))
    apply_theme(ax, fig, pal)

    non_out = plot_df[plot_df["OUTLIER"] != "YES"]
    out     = plot_df[plot_df["OUTLIER"] == "YES"]

    # Non-outlier samples
    ax.scatter(
        non_out["PC1"], non_out["PC2"],
        c     = pal["subtext"],
        s     = 6,
        alpha = 0.55,
        marker= "o",
        zorder= 3,
        label = f"Sample (n={len(non_out)})",
    )

    # Outlier samples — × marker
    if len(out) > 0:
        ax.scatter(
            out["PC1"], out["PC2"],
            c     = pal["fail"],
            s     = 20,
            alpha = 0.85,
            marker= "x",
            linewidths=1.2,
            zorder= 4,
            label = f"Outlier (n={len(out)})",
        )

    # Dashed ellipses per cluster (non-outliers only)
    for grp in all_groups:
        grp_data = plot_df[(plot_df["SUPERPOP"] == grp) &
                           (plot_df["OUTLIER"] != "YES")]
        if len(grp_data) < 5:
            continue
        confidence_ellipse(
            grp_data["PC1"].values,
            grp_data["PC2"].values,
            ax,
            sd        = args.ellipse_sd,
            colour    = pal["text"],
            fill      = False,
            linestyle = "--",
            linewidth = 0.7,
        )

    # Centroid labels
    centroids = compute_centroids(plot_df)
    for _, row in centroids.iterrows():
        ax.text(
            row["cx"], row["cy"], row["SUPERPOP"],
            fontsize   = 8,
            fontweight = "bold",
            color      = pal["text"],
            ha         = "center",
            va         = "center",
            zorder     = 5,
        )

    n_grp = len(all_groups)
    n_out = len(out)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_title(
        f"CRISP Step 6 — {args.plot_type} Pass 1 (Monochrome)",
        fontsize=13, fontweight="bold", color=pal["text"], loc="left"
    )
    ax.text(
        0.0, 1.01,
        f"N={len(plot_df)} samples — {n_grp} cluster(s) — "
        f"{n_out} outlier(s) — {pal['mode']} / {pal['background']}",
        transform=ax.transAxes, fontsize=8, color=pal["subtext"]
    )

    # Right-side legend — non-negotiable CRISP convention
    legend = ax.legend(
        loc            = "upper left",
        bbox_to_anchor = (1.01, 1.0),
        frameon        = True,
        framealpha     = 0.9,
        facecolor      = pal["bg"],
        edgecolor      = pal["grid"],
        fontsize       = 8,
        title          = "Status",
        title_fontsize = 9,
    )
    legend.get_title().set_color(pal["text"])
    for text in legend.get_texts():
        text.set_color(pal["text"])

    ax.text(
        0.0, -0.07,
        f"Confidence: PROBABLE — ellipse SD={args.ellipse_sd:.0f} — CRISP v1.0.0",
        transform=ax.transAxes, fontsize=7, color=pal["subtext"]
    )

    fig.savefig(out_path, bbox_inches="tight", dpi=200,
                facecolor=pal["bg"])
    plt.close(fig)
    print(f"[CRISP plot_pca_pass1.py] Plot 1 written: {out_path}")


##########################################################################
# PLOT 2 — COLOURED
##########################################################################

def plot_colour(plot_df: pd.DataFrame,
                pal: dict,
                args,
                x_label: str,
                y_label: str,
                all_groups: list,
                group_colours: dict,
                out_path: str) -> None:
    """
    Coloured scatter — one CRISP palette colour per group.
    Outliers always PAL_FAIL regardless of group.
    Coloured solid ellipses. Filled alpha overlay if COHORT_ADMIXED=YES.
    """
    fig, ax = plt.subplots(figsize=(9, 7))
    apply_theme(ax, fig, pal)

    # Plot non-outlier samples group by group
    for grp in all_groups:
        grp_data = plot_df[(plot_df["SUPERPOP"] == grp) &
                           (plot_df["OUTLIER"] != "YES")]
        if len(grp_data) == 0:
            continue
        ax.scatter(
            grp_data["PC1"], grp_data["PC2"],
            c     = group_colours[grp],
            s     = 6,
            alpha = 0.60,
            marker= "o",
            zorder= 3,
            label = f"{grp} (n={len(grp_data)})",
        )

    # Outliers — × marker, always PAL_FAIL
    out = plot_df[plot_df["OUTLIER"] == "YES"]
    if len(out) > 0:
        ax.scatter(
            out["PC1"], out["PC2"],
            c     = pal["fail"],
            s     = 22,
            alpha = 0.90,
            marker= "x",
            linewidths=1.3,
            zorder= 5,
            label = f"Outlier (n={len(out)})",
        )

    # Ellipses per group
    for grp in all_groups:
        grp_data = plot_df[(plot_df["SUPERPOP"] == grp) &
                           (plot_df["OUTLIER"] != "YES")]
        if len(grp_data) < 5:
            continue
        gcol = group_colours[grp]

        # Solid border ellipse
        confidence_ellipse(
            grp_data["PC1"].values,
            grp_data["PC2"].values,
            ax,
            sd        = args.ellipse_sd,
            colour    = gcol,
            fill      = False,
            linestyle = "-",
            linewidth = 0.9,
        )

        # WHY: When COHORT_ADMIXED=YES we additionally draw a filled
        # ellipse at low alpha so the cluster boundary is emphasised
        # visually — the analyst is looking at a single admixed group
        # and the filled ellipse makes the core vs periphery distinction
        # immediately obvious.
        if args.admixed.upper() == "YES":
            confidence_ellipse(
                grp_data["PC1"].values,
                grp_data["PC2"].values,
                ax,
                sd         = args.ellipse_sd,
                colour     = gcol,
                fill       = True,
                fill_alpha = 0.08,
                linestyle  = "-",
                linewidth  = 0.0,
            )

    # Centroid labels
    centroids = compute_centroids(plot_df)
    for _, row in centroids.iterrows():
        ax.text(
            row["cx"], row["cy"], row["SUPERPOP"],
            fontsize   = 8,
            fontweight = "bold",
            color      = pal["text"],
            ha         = "center",
            va         = "center",
            zorder     = 6,
        )

    n_grp = len(all_groups)
    n_out = len(out)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.set_title(
        f"CRISP Step 6 — {args.plot_type} Pass 1 (Coloured)",
        fontsize=13, fontweight="bold", color=pal["text"], loc="left"
    )
    ax.text(
        0.0, 1.01,
        f"N={len(plot_df)} samples — {n_grp} cluster(s) — "
        f"{n_out} outlier(s) — {pal['mode']} / {pal['background']}",
        transform=ax.transAxes, fontsize=8, color=pal["subtext"]
    )

    # Right-side legend
    legend = ax.legend(
        loc            = "upper left",
        bbox_to_anchor = (1.01, 1.0),
        frameon        = True,
        framealpha     = 0.9,
        facecolor      = pal["bg"],
        edgecolor      = pal["grid"],
        fontsize       = 8,
        title          = "Group",
        title_fontsize = 9,
        markerscale    = 2.0,
    )
    legend.get_title().set_color(pal["text"])
    for text in legend.get_texts():
        text.set_color(pal["text"])

    ax.text(
        0.0, -0.07,
        f"Confidence: PROBABLE — ellipse SD={args.ellipse_sd:.0f} — "
        f"Palette: {pal['mode']} — CRISP v1.0.0",
        transform=ax.transAxes, fontsize=7, color=pal["subtext"]
    )

    fig.savefig(out_path, bbox_inches="tight", dpi=200,
                facecolor=pal["bg"])
    plt.close(fig)
    print(f"[CRISP plot_pca_pass1.py] Plot 2 written: {out_path}")


##########################################################################
# SCREE PLOT
##########################################################################

def plot_scree(eigenval_path: str,
               pal: dict,
               args,
               out_path: str) -> None:
    """
    Scree plot — bars per PC (% variance) + cumulative line.
    WHY: Helps analyst decide how many PCs to include as GWAS covariates.
    Bars = PAL_PASS, cumulative line = PAL_WARN — consistent with R version.
    """
    if not os.path.exists(eigenval_path):
        print(f"[CRISP plot_pca_pass1.py] Eigenval not found, scree skipped: "
              f"{eigenval_path}")
        return

    eigenvals = np.array([float(l.strip()) for l in open(eigenval_path)
                          if l.strip()])
    pct_var  = eigenvals / eigenvals.sum() * 100
    cum_var  = np.cumsum(pct_var)
    n_shown  = min(20, len(eigenvals))
    pcs      = np.arange(1, n_shown + 1)

    fig, ax = plt.subplots(figsize=(9, 5))
    apply_theme(ax, fig, pal)

    ax.bar(pcs, pct_var[:n_shown],
           color=pal["pass"], alpha=0.80, width=0.65, zorder=3)

    ax2 = ax.twinx()
    ax2.plot(pcs, cum_var[:n_shown],
             color=pal["warn"], linewidth=1.2, marker="o",
             markersize=4, zorder=4)
    ax2.set_ylabel("Cumulative % Variance", color=pal["warn"], fontsize=9)
    ax2.tick_params(colors=pal["warn"], labelsize=8)
    ax2.set_facecolor("none")

    ax.set_xlabel("Principal Component", fontsize=10)
    ax.set_ylabel("% Variance Explained", fontsize=10)
    ax.set_xticks(pcs)
    ax.set_title(
        "CRISP Step 6 — PCA Scree Plot",
        fontsize=13, fontweight="bold", color=pal["text"], loc="left"
    )
    ax.text(
        0.0, 1.01,
        f"Project: {args.project} — {pal['mode']} / {pal['background']}",
        transform=ax.transAxes, fontsize=8, color=pal["subtext"]
    )
    ax.text(
        0.0, -0.10,
        "Bars = per-PC variance — Line = cumulative variance — CRISP v1.0.0",
        transform=ax.transAxes, fontsize=7, color=pal["subtext"]
    )

    fig.savefig(out_path, bbox_inches="tight", dpi=200,
                facecolor=pal["bg"])
    plt.close(fig)
    print(f"[CRISP plot_pca_pass1.py] Scree plot written: {out_path}")


##########################################################################
# MAIN
##########################################################################

def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP plot_pca_pass1.py] Project  : {args.project}")
    print(f"[CRISP plot_pca_pass1.py] Type     : {args.plot_type}")

    # Build palette
    pal = build_palette()
    print(f"[CRISP plot_pca_pass1.py] Palette  : {pal['mode']} / {pal['background']}")

    # Load data
    plot_df, x_label, y_label = load_plot_data(args)
    print(f"[CRISP plot_pca_pass1.py] Samples  : {len(plot_df)}")

    # Group setup
    all_groups    = sorted(plot_df["SUPERPOP"].unique().tolist())
    colour_seq    = group_colour_sequence(pal)
    group_colours = {
        grp: colour_seq[i % len(colour_seq)]
        for i, grp in enumerate(all_groups)
    }

    print(f"[CRISP plot_pca_pass1.py] Groups   : {all_groups}")

    plot_type_lower = args.plot_type.lower()

    # ── Plot 1 — Monochrome ───────────────────────────────────────────
    out_mono = out_dir / f"{args.project}_{plot_type_lower}_mono_py.pdf"
    plot_mono(plot_df, pal, args, x_label, y_label,
              all_groups, str(out_mono))

    # ── Plot 2 — Coloured ─────────────────────────────────────────────
    out_colour = out_dir / f"{args.project}_{plot_type_lower}_colour_py.pdf"
    plot_colour(plot_df, pal, args, x_label, y_label,
                all_groups, group_colours, str(out_colour))

    # ── Scree plot (PCA only) ─────────────────────────────────────────
    if args.plot_type == "PCA":
        eigenval_path = args.input.replace(".eigenvec", ".eigenval")
        out_scree = out_dir / f"{args.project}_pca_scree_py.pdf"
        plot_scree(eigenval_path, pal, args, str(out_scree))

    print(f"[CRISP plot_pca_pass1.py] All plots complete — "
          f"{pal['mode']} / {pal['background']}")


if __name__ == "__main__":
    main()
