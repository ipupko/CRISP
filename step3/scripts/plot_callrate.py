#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 3: Sample Call Rate, Plotting Script
# Version: 0.3.5
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
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
# PALETTE (from environment, with safe STANDARD/light fallbacks)
# ─────────────────────────────────────────────
def env_or(name: str, default: str) -> str:
    val = os.environ.get(name, '')
    return val if val else default

PAL_PASS    = env_or('CRISP_PAL_PASS',    '#1D9E75')
PAL_FAIL    = env_or('CRISP_PAL_FAIL',    '#ff5f57')
PAL_BG      = env_or('CRISP_PAL_BG',      '#FFFFFF')
PAL_PANEL   = env_or('CRISP_PAL_PANEL',   '#f8f9fa')
PAL_TEXT    = env_or('CRISP_PAL_TEXT',    '#1a1a1a')
PAL_SUBTEXT = env_or('CRISP_PAL_SUBTEXT', '#555555')
PAL_GRID    = env_or('CRISP_PAL_GRID',    '#e0e0e0')
PAL_MODE    = env_or('CRISP_PAL_MODE',    'STANDARD')
PAL_BACKGROUND = env_or('CRISP_PAL_BACKGROUND', 'LIGHT')

COLORS = {
    'Passing' : PAL_PASS,
    'Failing' : PAL_FAIL,
}

CAPTION = "CRISP | Comprehensive Robust Integrated SNP Processing"

# apply background/text colours globally via rcParams so every figure
# and axes created after this point picks them up automatically
plt.rcParams.update({
    'figure.facecolor': PAL_BG,
    'savefig.facecolor': PAL_BG,
    'axes.facecolor': PAL_PANEL,
    'axes.edgecolor': PAL_SUBTEXT,
    'axes.labelcolor': PAL_TEXT,
    'text.color': PAL_TEXT,
    'xtick.color': PAL_SUBTEXT,
    'ytick.color': PAL_SUBTEXT,
    'grid.color': PAL_GRID,
})


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
    ax.set_facecolor(PAL_PANEL)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.spines[['left', 'bottom']].set_color(PAL_SUBTEXT)
    ax.tick_params(labelsize=9, colors=PAL_SUBTEXT)


def add_legend(ax, loc='upper right'):
    patches = [mpatches.Patch(color=COLORS[s], label=s) for s in COLORS]
    legend = ax.legend(handles=patches, loc=loc, fontsize=9,
                        framealpha=0.9, edgecolor=PAL_GRID,
                        facecolor=PAL_BG)
    for text in legend.get_texts():
        text.set_color(PAL_TEXT)


def fmt(n: int) -> str:
    return f"{n:,}"


def x_axis_limit(df: pd.DataFrame, threshold: float) -> float:
    """Sensible x-axis upper limit: at least 2x the threshold, but never
    less than the max observed F_MISS * 1.05 (so points aren't clipped)."""
    return max(df['F_MISS'].max() * 1.05, threshold * 2)


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
                alpha=0.85, edgecolor=PAL_BG, linewidth=0.2,
                label=status)
    ax.axvline(threshold, color=PAL_FAIL, linestyle='--', linewidth=1.5)
    xlim_max = x_axis_limit(df, threshold)
    ax.set_xlim(0, xlim_max)
    # draw the histogram (and set xlim) BEFORE reading ylim, so the
    # annotation is placed relative to the actual plotted data
    ymax = ax.get_ylim()[1]
    label_offset = xlim_max * 0.01  # 1% of axis width, clears the bars
    ax.text(threshold + label_offset, ymax * 0.93,
            f'MIND = {threshold}',
            color=PAL_FAIL, fontsize=9, fontweight='bold')
    ax.set_xlabel('Per-sample missingness rate (F_MISS)',
                  fontweight='bold', fontsize=10)
    ax.set_ylabel('Number of samples', fontweight='bold', fontsize=10)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'{x*100:.1f}%'))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        lambda y, _: fmt(int(y))))
    apply_style(ax)
    add_legend(ax, loc='upper right')
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
                 fontweight='bold', fontsize=12, color=PAL_TEXT)
    fig.suptitle(
        f"SIMPLE mode  |  MIND threshold: {mind}  |  "
        f"{fmt(n_failing)}/{fmt(n_total)} samples failing",
        fontsize=9, color=PAL_SUBTEXT, y=0.96
    )
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_SUBTEXT)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    out = os.path.join(report_dir, "step3_callrate_simple.pdf")
    fig.savefig(out, bbox_inches='tight', dpi=150, facecolor=PAL_BG)
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

    # Load each .imiss file once and reuse for both the individual plot
    # and the faceted plot (v0.2.0 loaded each file twice)
    dataframes = []
    for filepath, threshold in zip(imiss_files, thresholds):
        df = load_imiss(filepath)
        df['status'] = df['F_MISS'].apply(
            lambda x: 'Failing' if x > threshold else 'Passing'
        )
        dataframes.append(df)

    # Individual plots
    for i, (filepath, threshold, label, df) in enumerate(
            zip(imiss_files, thresholds, labels, dataframes)):

        log(f"\nPass {i+1} , Loading: {filepath}")
        n_total   = len(df)
        n_failing = int((df['status'] == 'Failing').sum())

        fig, ax = plt.subplots(figsize=(10, 6))
        draw_histogram(ax, df, threshold)
        ax.set_title(f"CRISP: Step 3 Sample Call Rate , {label}",
                     fontweight='bold', fontsize=12, color=PAL_TEXT)
        fig.suptitle(
            f"{mode} mode  |  "
            f"{fmt(n_failing)}/{fmt(n_total)} samples failing at this threshold",
            fontsize=9, color=PAL_SUBTEXT, y=0.96
        )
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_SUBTEXT)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        out = os.path.join(
            report_dir,
            f"step3_callrate_{mode.lower()}_pass{i+1}.pdf"
        )
        fig.savefig(out, bbox_inches='tight', dpi=150, facecolor=PAL_BG)
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

    # shared formatter functions (avoids the per-iteration lambda
    # closure-over-loop-variable bug from v0.2.0)
    def pct_formatter(x, _pos):
        return f'{x*100:.0f}%'

    def count_formatter(y, _pos):
        return fmt(int(y))

    for i, (filepath, threshold, label, df) in enumerate(
            zip(imiss_files, thresholds, labels, dataframes)):

        ax = axes_flat[i]
        for status, grp in df.groupby('status'):
            ax.hist(grp['F_MISS'], bins=40, color=COLORS[status],
                    alpha=0.85, edgecolor=PAL_BG, linewidth=0.15)
        ax.axvline(threshold, color=PAL_FAIL,
                   linestyle='--', linewidth=1.0)
        ax.set_xlim(0, x_axis_limit(df, threshold))
        ax.set_title(label, fontweight='bold', fontsize=10, color=PAL_TEXT)
        ax.set_xlabel('F_MISS', fontsize=8)
        ax.set_ylabel('Samples', fontsize=8)
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(pct_formatter))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(count_formatter))
        apply_style(ax)
        n_fail = int((df['status'] == 'Failing').sum())
        ax.text(0.98, 0.95, f'{fmt(n_fail)} failing',
                transform=ax.transAxes, ha='right', va='top',
                color=PAL_FAIL, fontsize=8, fontweight='bold')

    # Hide unused panels
    for j in range(n, len(axes_flat)):
        axes_flat[j].set_visible(False)

    # legend placed upper-right of the figure, consistent with the
    # per-pass plots (v0.2.0 used lower-center here only)
    patches = [mpatches.Patch(color=COLORS[s], label=s) for s in COLORS]
    legend = fig.legend(handles=patches, loc='upper right', ncol=1,
                         fontsize=9, bbox_to_anchor=(0.995, 0.98),
                         facecolor=PAL_BG, edgecolor=PAL_GRID)
    for text in legend.get_texts():
        text.set_color(PAL_TEXT)

    fig.suptitle(
        f"CRISP: Step 3 Sample Call Rate , {mode} Mode\n"
        f"Per-sample missingness distribution across {n} threshold passes",
        fontweight='bold', fontsize=12, color=PAL_TEXT
    )
    fig.text(0.5, -0.01, CAPTION, ha='center', fontsize=8, color=PAL_SUBTEXT)
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])

    facet = os.path.join(
        report_dir,
        f"step3_callrate_{mode.lower()}_faceted.pdf"
    )
    fig.savefig(facet, bbox_inches='tight', dpi=150, facecolor=PAL_BG)
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
    log(f"Colour mode: {PAL_MODE}")
    log(f"Background : {PAL_BACKGROUND}")

    if mode == "SIMPLE":
        plot_simple(imiss_files[0], mind, report_dir)

    elif mode == "CASCADE":
        if n_passes != 4:
            abort(f"CASCADE expects 4 .imiss files, got {n_passes}")
        # thresholds parsed from filenames (mind{value} pattern), with
        # the documented CASCADE defaults as fallback if unparseable
        cascade_defaults = [0.25, 0.20, 0.10, 0.05]
        thresholds = [
            parse_threshold(f, fallback=cascade_defaults[i])
            for i, f in enumerate(imiss_files)
        ]
        plot_multipass(imiss_files, thresholds, "CASCADE", report_dir)

    elif mode == "CUSTOM":
        if n_passes < 2:
            abort(f"CUSTOM mode requires at least 2 .imiss files, got {n_passes}")
        # thresholds parsed from filenames; fallback descends from `mind`
        # in 0.05 steps for any pass whose filename doesn't match mind{value}
        thresholds = [
            parse_threshold(f, fallback=mind + (n_passes - i - 1) * 0.05)
            for i, f in enumerate(imiss_files)
        ]
        plot_multipass(imiss_files, thresholds, "CUSTOM", report_dir)

    else:
        abort(f"Unknown mode '{mode}'. Valid: SIMPLE, CASCADE, CUSTOM")

    log("\nPlotting complete.")


if __name__ == "__main__":
    main()
