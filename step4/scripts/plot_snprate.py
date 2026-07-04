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
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 4: Variant Call Rate — Plotting Script
# Version: 0.3.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# Python alternative to plot_snprate.R for users who prefer
# matplotlib over ggplot2. Selected via PLOT_ENGINE = PYTHON.
#
# SIMPLE mode:
#   Single histogram of per-variant missingness rates
#   Vertical line marking the GENO threshold
#   Variants failing threshold highlighted in red
#
# CASCADE mode:
#   Four individual histograms, one per pass
#   Single faceted 2x2 plot showing all passes
#
# CUSTOM mode:
#   One individual histogram per tier
#   Single faceted plot, 2 columns, rows scale to tier count
#   Last cell hidden if odd number of tiers
#
# Usage:
#   python3 plot_snprate.py <mode> <geno> <report_dir> <lmiss_files...>
##########################################################################

import sys
import os
import re
import math
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker


# ─────────────────────────────────────────────
# COLOUR PALETTE
# Reads CRISP_PAL_* environment variables exported by _crisp_flavour.sh's
# _crisp_palette() function in crisp_snprate.sh.
#
# WHY ENV VARS: Python cannot source a bash script. Environment variable
# passing is the only clean cross-language mechanism. The _crisp_palette()
# call in crisp_snprate.sh runs before any python3 call, so all vars
# are available in the Python process environment.
#
# WHY FALLBACK VALUES: Hardcoded fallbacks match the previous hardcoded
# values exactly, so STANDARD mode produces byte-identical output to the
# pre-v0.3.0 script. Zero-risk upgrade for users not using COLOURBLIND
# or NIGHT mode.
#
# NOTE: CRISP_PAL_SUBTITLE and CRISP_PAL_CAPTION are not yet exported
# by _crisp_flavour.sh. Fallbacks are used until they are added.
# ─────────────────────────────────────────────

_COLOUR_MODE = os.environ.get("CRISP_PAL_MODE", "STANDARD").upper()

if _COLOUR_MODE == "COLOURBLIND":
    PAL_PASS     = os.environ.get("CRISP_PAL_PASS",     "#0072B2")
    PAL_FAIL     = os.environ.get("CRISP_PAL_FAIL",     "#D55E00")
    PAL_SUBTITLE = os.environ.get("CRISP_PAL_SUBTITLE", "#555555")
    PAL_CAPTION  = os.environ.get("CRISP_PAL_CAPTION",  "#888888")
    BG_COLOUR    = "#FFFFFF"
    FG_COLOUR    = "#000000"
    AX_COLOUR    = "#f8f9fa"
elif _COLOUR_MODE == "NIGHT":
    # Dark-background variant, Okabe-Ito-derived.
    # Provisional pending final NIGHT palette in _crisp_flavour.sh.
    PAL_PASS     = os.environ.get("CRISP_PAL_PASS",     "#56B4E9")
    PAL_FAIL     = os.environ.get("CRISP_PAL_FAIL",     "#E69F00")
    PAL_SUBTITLE = os.environ.get("CRISP_PAL_SUBTITLE", "#AAAAAA")
    PAL_CAPTION  = os.environ.get("CRISP_PAL_CAPTION",  "#777777")
    BG_COLOUR    = "#1E1E1E"
    FG_COLOUR    = "#E0E0E0"
    AX_COLOUR    = "#2A2A2A"
else:
    # STANDARD — fallbacks match pre-v0.3.0 hardcoded values exactly
    PAL_PASS     = os.environ.get("CRISP_PAL_PASS",     "#1D9E75")
    PAL_FAIL     = os.environ.get("CRISP_PAL_FAIL",     "#ff5f57")
    PAL_SUBTITLE = os.environ.get("CRISP_PAL_SUBTITLE", "#555555")
    PAL_CAPTION  = os.environ.get("CRISP_PAL_CAPTION",  "#888888")
    BG_COLOUR    = "#FFFFFF"
    FG_COLOUR    = "#000000"
    AX_COLOUR    = "#f8f9fa"

COLORS = {
    'Passing': PAL_PASS,
    'Failing': PAL_FAIL,
}

CAPTION = "CRISP | Comprehensive Robust Integrated SNP Processing"


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
def log(msg):
    print(f"[PLOT] {msg}", flush=True)

def abort(msg):
    print(f"[PLOT] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# LOAD LMISS FILE
# ─────────────────────────────────────────────
def load_lmiss(filepath: str) -> pd.DataFrame:
    """Load a PLINK .lmiss file. Columns: CHR SNP N_MISS N_GENO F_MISS"""
    if not os.path.isfile(filepath):
        abort(f".lmiss file not found: {filepath}")
    df = pd.read_csv(filepath, sep=r'\s+')
    if 'F_MISS' not in df.columns:
        abort(f"F_MISS column not found in: {filepath}")
    return df


# ─────────────────────────────────────────────
# SHARED STYLE
# ─────────────────────────────────────────────
def apply_style(ax):
    # WHY AX_COLOUR: axis background adapts to NIGHT mode (#2A2A2A) while
    # remaining the standard light grey (#f8f9fa) in STANDARD/COLOURBLIND.
    ax.set_facecolor(AX_COLOUR)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.spines[['left', 'bottom']].set_color(FG_COLOUR)
    ax.tick_params(labelsize=9, colors=FG_COLOUR)
    ax.xaxis.label.set_color(FG_COLOUR)
    ax.yaxis.label.set_color(FG_COLOUR)


def add_legend(ax):
    # WHY loc='right': CRISP plotting standard. Right-side legends do not
    # obscure the histogram body and scale gracefully with label length.
    patches = [mpatches.Patch(color=COLORS[s], label=s) for s in COLORS]
    ax.legend(handles=patches, loc='upper right', fontsize=9,
              framealpha=0.9, edgecolor=FG_COLOUR,
              labelcolor=FG_COLOUR, facecolor=BG_COLOUR)


def format_count(n: int) -> str:
    return f"{n:,}"


# ─────────────────────────────────────────────
# SINGLE HISTOGRAM
# ─────────────────────────────────────────────
def draw_histogram(ax, df: pd.DataFrame, threshold: float,
                   title: str, subtitle: str):
    """Draw a single missingness histogram on ax."""
    df = df.copy()
    df['status'] = df['F_MISS'].apply(
        lambda x: 'Failing' if x > threshold else 'Passing'
    )

    for status, grp in df.groupby('status'):
        ax.hist(grp['F_MISS'], bins=80, color=COLORS[status],
                alpha=0.85, edgecolor='white', linewidth=0.2,
                label=status)

    ax.axvline(threshold, color=PAL_FAIL, linestyle='--', linewidth=1.5)
    ax.text(threshold + 0.002, ax.get_ylim()[1] * 0.94,
            f'GENO = {threshold}',
            color=PAL_FAIL, fontsize=9, fontweight='bold')

    ax.set_xlabel('Per-variant missingness rate (F_MISS)',
                  fontweight='bold', fontsize=10, color=FG_COLOUR)
    ax.set_ylabel('Number of variants', fontweight='bold', fontsize=10,
                  color=FG_COLOUR)
    ax.set_title(title, fontweight='bold', fontsize=12, pad=8, color=FG_COLOUR)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'{x*100:.1f}%'))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: format_count(int(x))))
    apply_style(ax)
    add_legend(ax)

    n_total   = len(df)
    n_failing = (df['F_MISS'] > threshold).sum()

    return n_total, n_failing


# ─────────────────────────────────────────────
# SIMPLE MODE
# ─────────────────────────────────────────────
def plot_simple(lmiss_file: str, geno: float, report_dir: str):

    log(f"Loading: {lmiss_file}")
    df = load_lmiss(lmiss_file)

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor(BG_COLOUR)

    n_total   = len(df)
    n_failing = int((df['F_MISS'] > geno).sum())

    draw_histogram(
        ax, df, geno,
        title    = "CRISP — Step 4: Variant Call Rate",
        subtitle = f"SIMPLE mode  |  GENO={geno:.4g}  |  "
                   f"{format_count(n_failing)}/{format_count(n_total)} failing"
    )

    fig.suptitle(
        f"SIMPLE mode  |  GENO threshold: {geno:.4g}  |  "
        f"{format_count(n_failing)}/{format_count(n_total)} variants failing",
        fontsize=9, color=PAL_SUBTITLE, y=0.96
    )
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_CAPTION)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    out_file = os.path.join(report_dir, "step4_snprate_simple.pdf")
    fig.savefig(out_file, bbox_inches='tight', dpi=150)
    plt.close(fig)

    log(f"Simple histogram saved: {out_file}")
    log(f"  Variants total   : {format_count(n_total)}")
    log(f"  Variants failing : {format_count(n_failing)} (F_MISS > {geno})")
    log(f"  Variants passing : {format_count(n_total - n_failing)}")

    return out_file


# ─────────────────────────────────────────────
# MULTI-PASS MODE (CASCADE and CUSTOM)
# ─────────────────────────────────────────────
def plot_multipass(lmiss_files: list, thresholds: list,
                   mode: str, report_dir: str):

    n         = len(lmiss_files)
    labels    = [f"Pass {i+1} (geno={t:.4g})" for i, t in enumerate(thresholds)]
    pdf_files = []

    # ── Individual plots ─────────────────────────────────────
    for i, (filepath, threshold, label) in enumerate(
            zip(lmiss_files, thresholds, labels)):

        log(f"\nPass {i+1} — Loading: {filepath}")
        df = load_lmiss(filepath)

        fig, ax = plt.subplots(figsize=(10, 6))
        fig.patch.set_facecolor(BG_COLOUR)
        n_total   = len(df)
        n_failing = int((df['F_MISS'] > threshold).sum())

        draw_histogram(
            ax, df, threshold,
            title    = f"CRISP — Step 4: Variant Call Rate — {label}",
            subtitle = f"{mode} mode  |  "
                       f"{format_count(n_failing)}/{format_count(n_total)} variants failing"
        )
        fig.suptitle(
            f"{mode} mode  |  "
            f"{format_count(n_failing)}/{format_count(n_total)} variants failing at this threshold",
            fontsize=9, color=PAL_SUBTITLE, y=0.96
        )
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_CAPTION)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        out_file = os.path.join(
            report_dir,
            f"step4_snprate_{mode.lower()}_pass{i+1}.pdf"
        )
        fig.savefig(out_file, bbox_inches='tight', dpi=150)
        plt.close(fig)

        pdf_files.append(out_file)
        log(f"  Individual plot saved: {out_file}")
        log(f"  Threshold : {threshold}")
        log(f"  Variants  : {format_count(n_total)} total, "
            f"{format_count(n_failing)} failing")

    # ── Faceted plot ──────────────────────────────────────────
    log(f"\nGenerating faceted {mode} comparison plot ({n} passes)...")

    n_cols  = 2
    n_rows  = math.ceil(n / n_cols)
    fig_w   = 14
    fig_h   = max(6, n_rows * 4.5)

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(fig_w, fig_h),
                              squeeze=False)
    fig.patch.set_facecolor(BG_COLOUR)
    axes_flat = axes.flatten()

    for i, (filepath, threshold, label) in enumerate(
            zip(lmiss_files, thresholds, labels)):

        df = load_lmiss(filepath)
        df['status'] = df['F_MISS'].apply(
            lambda x: 'Failing' if x > threshold else 'Passing'
        )

        ax = axes_flat[i]
        for status, grp in df.groupby('status'):
            ax.hist(grp['F_MISS'], bins=60, color=COLORS[status],
                    alpha=0.85, edgecolor='white', linewidth=0.15)

        ax.axvline(threshold, color=PAL_FAIL,
                   linestyle='--', linewidth=1.0)
        ax.set_title(label, fontweight='bold', fontsize=10, color=FG_COLOUR)
        ax.set_xlabel('F_MISS', fontsize=8, color=FG_COLOUR)
        ax.set_ylabel('Variants', fontsize=8, color=FG_COLOUR)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{x*100:.0f}%'))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: format_count(int(x))))
        apply_style(ax)

        n_fail = (df['F_MISS'] > threshold).sum()
        ax.text(0.98, 0.95, f'{format_count(n_fail)} failing',
                transform=ax.transAxes, ha='right', va='top',
                color=PAL_FAIL, fontsize=8, fontweight='bold')

    # Hide unused panels
    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    # Shared legend — placed below plot area, clear of caption text
    patches = [mpatches.Patch(color=COLORS[s], label=s) for s in COLORS]
    fig.legend(handles=patches, loc='lower center', ncol=2,
               fontsize=9, bbox_to_anchor=(0.5, 0.0),
               framealpha=0.9, edgecolor=FG_COLOUR,
               labelcolor=FG_COLOUR, facecolor=BG_COLOUR)

    fig.suptitle(
        f"CRISP — Step 4: Variant Call Rate — {mode} Mode\n"
        f"Per-variant missingness distribution across {n} threshold passes",
        fontweight='bold', fontsize=12, color=FG_COLOUR
    )
    fig.text(0.5, -0.01, CAPTION, ha='center', fontsize=8, color=PAL_CAPTION)
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])

    facet_file = os.path.join(
        report_dir,
        f"step4_snprate_{mode.lower()}_faceted.pdf"
    )
    fig.savefig(facet_file, bbox_inches='tight', dpi=150)
    plt.close(fig)

    log(f"Faceted plot saved: {facet_file}")

    return pdf_files + [facet_file]


# ─────────────────────────────────────────────
# PARSE THRESHOLD FROM FILENAME
# ─────────────────────────────────────────────
def parse_threshold_from_filename(filepath: str, fallback: float) -> float:
    """
    Extract threshold from filename pattern geno{value}.lmiss
    e.g. step4_custom_pass1_geno0.30.lmiss -> 0.30
    """
    fname = os.path.basename(filepath)
    match = re.search(r'geno([0-9]+\.?[0-9]*)', fname)
    if match:
        return float(match.group(1))
    return fallback


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRISP — Step 4: Variant Call Rate Plotting (Python)"
    )
    parser.add_argument("mode",        help="SIMPLE, CASCADE, or CUSTOM")
    parser.add_argument("geno",        help="Final GENO threshold", type=float)
    parser.add_argument("report_dir",  help="Output directory for plots")
    parser.add_argument("lmiss_files", nargs='+',
                        help="One or more .lmiss files")
    args = parser.parse_args()

    mode        = args.mode.upper()
    geno        = args.geno
    report_dir  = args.report_dir
    lmiss_files = args.lmiss_files
    n_passes    = len(lmiss_files)

    os.makedirs(report_dir, exist_ok=True)

    log("CRISP — Step 4: Variant Call Rate Plotting (Python)")
    log(f"Mode       : {mode}")
    log(f"GENO       : {geno}")
    log(f"Report dir : {report_dir}")
    log(f"Files      : {n_passes}")

    if mode == "SIMPLE":
        plot_simple(lmiss_files[0], geno, report_dir)

    elif mode == "CASCADE":
        thresholds = [0.25, 0.20, 0.10, 0.05]
        if n_passes != 4:
            abort(f"CASCADE expects 4 .lmiss files, got {n_passes}")
        plot_multipass(lmiss_files, thresholds, "CASCADE", report_dir)

    elif mode == "CUSTOM":
        # Infer thresholds from filenames
        thresholds = []
        for i, f in enumerate(lmiss_files):
            fallback = geno + (n_passes - i - 1) * 0.05
            thresholds.append(parse_threshold_from_filename(f, fallback))

        if len(thresholds) < 2:
            abort("CUSTOM mode requires at least 2 .lmiss files")

        plot_multipass(lmiss_files, thresholds, "CUSTOM", report_dir)

    else:
        abort(f"Unknown mode '{mode}'. Valid options: SIMPLE, CASCADE, CUSTOM")

    log("\nPlotting complete.")


if __name__ == "__main__":
    main()
