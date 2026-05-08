#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 3: Sample Call Rate, Plotting Script
# Version: 0.2.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
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
# CONSTANTS
# ─────────────────────────────────────────────
COLORS = {
    'Passing' : '#1D9E75',
    'Failing' : '#ff5f57',
}

CAPTION = "CRISP | Comprehensive Robust Integrated SNP Processing"


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
def log(msg: str):
    print(f"[PLOT] {msg}", flush=True)

def abort(msg: str):
    print(f"[PLOT] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# LOAD IMISS FILE
# ─────────────────────────────────────────────
def load_imiss(filepath: str) -> pd.DataFrame:
    """
    Load a PLINK .imiss file.
    Expected columns: FID IID MISS_PHENO N_MISS N_GENO F_MISS
    """
    if not os.path.isfile(filepath):
        abort(f".imiss file not found: {filepath}")
    df = pd.read_csv(filepath, sep=r'\s+')
    if 'F_MISS' not in df.columns:
        abort(f"F_MISS column not found in: {filepath}")
    return df


# ─────────────────────────────────────────────
# SHARED STYLE
# ─────────────────────────────────────────────
def apply_style(ax):
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.tick_params(labelsize=9)


def add_legend(ax):
    patches = [mpatches.Patch(color=COLORS[s], label=s) for s in COLORS]
    ax.legend(handles=patches, loc='upper right', fontsize=9,
              framealpha=0.9, edgecolor='#cccccc')


def fmt(n: int) -> str:
    return f"{n:,}"


# ─────────────────────────────────────────────
# DRAW SINGLE HISTOGRAM
# ─────────────────────────────────────────────
def draw_histogram(ax, df: pd.DataFrame, threshold: float):
    """Draw per-sample missingness histogram on ax."""
    df = df.copy()
    df['status'] = df['F_MISS'].apply(
        lambda x: 'Failing' if x > threshold else 'Passing'
    )
    for status, grp in df.groupby('status'):
        ax.hist(grp['F_MISS'], bins=50, color=COLORS[status],
                alpha=0.85, edgecolor='white', linewidth=0.2,
                label=status)
    ax.axvline(threshold, color='#ff5f57', linestyle='--', linewidth=1.5)
    ax.text(threshold + 0.001, ax.get_ylim()[1] * 0.93,
            f'MIND = {threshold}',
            color='#ff5f57', fontsize=9, fontweight='bold')
    ax.set_xlabel('Per-sample missingness rate (F_MISS)',
                  fontweight='bold', fontsize=10)
    ax.set_ylabel('Number of samples', fontweight='bold', fontsize=10)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'{x*100:.1f}%'))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: fmt(int(x))))
    apply_style(ax)
    add_legend(ax)
    n_total   = len(df)
    n_failing = int((df['F_MISS'] > threshold).sum())
    return n_total, n_failing


# ─────────────────────────────────────────────
# SIMPLE MODE
# ─────────────────────────────────────────────
def plot_simple(imiss_file: str, mind: float, report_dir: str):

    log(f"Loading: {imiss_file}")
    df        = load_imiss(imiss_file)
    n_total   = len(df)
    n_failing = int((df['F_MISS'] > mind).sum())

    fig, ax = plt.subplots(figsize=(10, 6))
    draw_histogram(ax, df, mind)
    ax.set_title("CRISP: Step 3 Sample Call Rate",
                 fontweight='bold', fontsize=12)
    fig.suptitle(
        f"SIMPLE mode  |  MIND threshold: {mind}  |  "
        f"{fmt(n_failing)}/{fmt(n_total)} samples failing",
        fontsize=9, color='#555555', y=0.96
    )
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    out = os.path.join(report_dir, "step3_callrate_simple.pdf")
    fig.savefig(out, bbox_inches='tight', dpi=150)
    plt.close(fig)

    log(f"Simple histogram saved: {out}")
    log(f"  Samples total   : {fmt(n_total)}")
    log(f"  Samples failing : {fmt(n_failing)} (F_MISS > {mind})")
    log(f"  Samples passing : {fmt(n_total - n_failing)}")
    return out


# ─────────────────────────────────────────────
# MULTI-PASS MODE (CASCADE and CUSTOM)
# ─────────────────────────────────────────────
def plot_multipass(imiss_files: list, thresholds: list,
                   mode: str, report_dir: str):

    n         = len(imiss_files)
    labels    = [f"Pass {i+1} (mind={t})" for i, t in enumerate(thresholds)]
    pdf_files = []

    # Individual plots
    for i, (filepath, threshold, label) in enumerate(
            zip(imiss_files, thresholds, labels)):

        log(f"\nPass {i+1} , Loading: {filepath}")
        df        = load_imiss(filepath)
        n_total   = len(df)
        n_failing = int((df['F_MISS'] > threshold).sum())

        fig, ax = plt.subplots(figsize=(10, 6))
        draw_histogram(ax, df, threshold)
        ax.set_title(f"CRISP: Step 3 Sample Call Rate , {label}",
                     fontweight='bold', fontsize=12)
        fig.suptitle(
            f"{mode} mode  |  "
            f"{fmt(n_failing)}/{fmt(n_total)} samples failing at this threshold",
            fontsize=9, color='#555555', y=0.96
        )
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        out = os.path.join(
            report_dir,
            f"step3_callrate_{mode.lower()}_pass{i+1}.pdf"
        )
        fig.savefig(out, bbox_inches='tight', dpi=150)
        plt.close(fig)

        pdf_files.append(out)
        log(f"  Individual plot saved: {out}")
        log(f"  Threshold : {threshold}")
        log(f"  Samples   : {fmt(n_total)} total, {fmt(n_failing)} failing")

    # Faceted plot
    log(f"\nGenerating faceted {mode} comparison plot ({n} passes)...")

    n_cols = 2
    n_rows = math.ceil(n / n_cols)
    fig_w  = 14
    fig_h  = max(6, n_rows * 4.5)

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(fig_w, fig_h),
                              squeeze=False)
    axes_flat = axes.flatten()

    for i, (filepath, threshold, label) in enumerate(
            zip(imiss_files, thresholds, labels)):

        df = load_imiss(filepath)
        df['status'] = df['F_MISS'].apply(
            lambda x: 'Failing' if x > threshold else 'Passing'
        )
        ax = axes_flat[i]
        for status, grp in df.groupby('status'):
            ax.hist(grp['F_MISS'], bins=40, color=COLORS[status],
                    alpha=0.85, edgecolor='white', linewidth=0.15)
        ax.axvline(threshold, color='#ff5f57',
                   linestyle='--', linewidth=1.0)
        ax.set_title(label, fontweight='bold', fontsize=10)
        ax.set_xlabel('F_MISS', fontsize=8)
        ax.set_ylabel('Samples', fontsize=8)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{x*100:.0f}%'))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: fmt(int(x))))
        apply_style(ax)
        n_fail = int((df['F_MISS'] > threshold).sum())
        ax.text(0.98, 0.95, f'{fmt(n_fail)} failing',
                transform=ax.transAxes, ha='right', va='top',
                color='#ff5f57', fontsize=8, fontweight='bold')

    # Hide unused panels
    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    patches = [mpatches.Patch(color=COLORS[s], label=s) for s in COLORS]
    fig.legend(handles=patches, loc='lower center', ncol=2,
               fontsize=9, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(
        f"CRISP: Step 3 Sample Call Rate , {mode} Mode\n"
        f"Per-sample missingness distribution across {n} threshold passes",
        fontweight='bold', fontsize=12
    )
    fig.text(0.5, -0.01, CAPTION, ha='center', fontsize=8, color='#888888')
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])

    facet = os.path.join(
        report_dir,
        f"step3_callrate_{mode.lower()}_faceted.pdf"
    )
    fig.savefig(facet, bbox_inches='tight', dpi=150)
    plt.close(fig)

    log(f"Faceted plot saved: {facet}")
    return pdf_files + [facet]


# ─────────────────────────────────────────────
# PARSE THRESHOLD FROM FILENAME
# ─────────────────────────────────────────────
def parse_threshold(filepath: str, fallback: float) -> float:
    """
    Extract threshold from filename pattern mind{value}.imiss
    e.g. step3_custom_pass1_mind0.20.imiss -> 0.20
    """
    match = re.search(r'mind([0-9]+\.?[0-9]*)',
                      os.path.basename(filepath))
    if match:
        return float(match.group(1))
    return fallback


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRISP-py Step 3: Sample Call Rate Plotting"
    )
    parser.add_argument("mode",        help="SIMPLE, CASCADE, or CUSTOM")
    parser.add_argument("mind",        help="Final MIND threshold", type=float)
    parser.add_argument("report_dir",  help="Output directory for plots")
    parser.add_argument("imiss_files", nargs='+',
                        help="One or more .imiss files")
    args = parser.parse_args()

    mode        = args.mode.upper()
    mind        = args.mind
    report_dir  = args.report_dir
    imiss_files = args.imiss_files
    n_passes    = len(imiss_files)

    os.makedirs(report_dir, exist_ok=True)

    log("CRISP-py Step 3: Sample Call Rate Plotting")
    log(f"Mode       : {mode}")
    log(f"MIND       : {mind}")
    log(f"Report dir : {report_dir}")
    log(f"Files      : {n_passes}")

    if mode == "SIMPLE":
        plot_simple(imiss_files[0], mind, report_dir)

    elif mode == "CASCADE":
        thresholds = [0.25, 0.20, 0.10, 0.05]
        if n_passes != 4:
            abort(f"CASCADE expects 4 .imiss files, got {n_passes}")
        plot_multipass(imiss_files, thresholds, "CASCADE", report_dir)

    elif mode == "CUSTOM":
        thresholds = []
        for i, f in enumerate(imiss_files):
            fallback = mind + (n_passes - i - 1) * 0.05
            thresholds.append(parse_threshold(f, fallback))
        if len(thresholds) < 2:
            abort("CUSTOM mode requires at least 2 .imiss files")
        plot_multipass(imiss_files, thresholds, "CUSTOM", report_dir)

    else:
        abort(f"Unknown mode '{mode}'. Valid: SIMPLE, CASCADE, CUSTOM")

    log("\nPlotting complete.")


if __name__ == "__main__":
    main()
