#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 8: HWE Filtering, Plotting Script
# Version: 0.1.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# Generates HWE p-value distribution plots from PLINK --hardy output.
#
# Plots:
#   Plot 1: Histogram of -log10(HWE p-values) with threshold line
#   Plot 2: QQ plot of observed vs expected -log10(p-values)
#
# Usage:
#   python3 plot_hwe.py <output_dir> <hwe_file> <threshold> <mode>
##########################################################################

import sys
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


##########################################################################
### COLOUR PALETTE
### Reads PLOT_COLOUR_MODE from environment (set by crisp_*.sh)
### STANDARD: default CRISP palette
### COLOURBLIND: Okabe-Ito palette (Nature Methods recommended)
##########################################################################

import os as _os
_mode = _os.getenv('CRISP_PAL_MODE', 'STANDARD').upper()

if _mode == 'COLOURBLIND':
    PAL_PASS      = _os.getenv('CRISP_PAL_PASS',      '#0072B2')
    PAL_FAIL      = _os.getenv('CRISP_PAL_FAIL',      '#D55E00')
    PAL_WARN      = _os.getenv('CRISP_PAL_WARN',      '#E69F00')
    PAL_HIGHLIGHT = _os.getenv('CRISP_PAL_HIGHLIGHT', '#009E73')
    PAL_SKY       = _os.getenv('CRISP_PAL_SKY',       '#56B4E9')
    PAL_PINK      = _os.getenv('CRISP_PAL_PINK',      '#CC79A7')
    PAL_ORANGE    = _os.getenv('CRISP_PAL_YELLOW',    '#F0E442')
else:
    PAL_PASS      = '#1D9E75'
    PAL_FAIL      = '#ff5f57'
    PAL_WARN      = '#febc2e'
    PAL_HIGHLIGHT = '#5dcaa5'
    PAL_SKY       = '#7ec8e3'
    PAL_PINK      = '#afa9ec'
    PAL_ORANGE    = '#ff8c00'

CAPTION = 'CRISP | Comprehensive Robust Integrated SNP Processing | Compass Genomics'


def log(msg):
    print(f"[PLOT] {msg}", flush=True)


def apply_style(ax):
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=9)


##########################################################################
### PLOT 1: -LOG10(P) HISTOGRAM
##########################################################################

def plot_histogram(df, output_dir, threshold, mode):

    threshold_log = -np.log10(threshold)
    n_fail        = (df['P'] < threshold).sum()

    fig, ax = plt.subplots(figsize=(10, 6))

    pass_mask = df['log10P'] < threshold_log
    fail_mask = ~pass_mask

    ax.hist(df.loc[pass_mask, 'log10P'],
            bins=80, color='#1D9E75', alpha=0.85,
            edgecolor='white', linewidth=0.15,
            label=f'Passing HWE ({pass_mask.sum():,})')

    if fail_mask.sum() > 0:
        ax.hist(df.loc[fail_mask, 'log10P'],
                bins=80, color='#ff5f57', alpha=0.85,
                edgecolor='white', linewidth=0.15,
                label=f'Failing HWE (p < {threshold:.0e}) ({fail_mask.sum():,})')

    ax.axvline(threshold_log, color='#ff5f57',
               linestyle='--', linewidth=1.2)
    ax.text(threshold_log + 0.1, ax.get_ylim()[1] * 0.92,
            f'-log10(p) = {threshold_log:.1f}',
            color='#ff5f57', fontsize=8.5, fontweight='bold')

    ax.set_xlabel(r'$-\log_{10}(p)$', fontweight='bold', fontsize=10)
    ax.set_ylabel('Number of variants', fontweight='bold', fontsize=10)
    ax.set_title('CRISP: Step 8 HWE, P-value Distribution',
                 fontweight='bold', fontsize=13)
    apply_style(ax)
    fig.suptitle(
        f'Mode: {mode}  |  Threshold: p < {threshold:.0e}  |  '
        f'{n_fail:,} failing HWE',
        fontsize=9, color='#555', y=0.97)
    ax.legend(fontsize=9, loc='center right',
              bbox_to_anchor=(1.18, 0.5), framealpha=0.9)
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    out = os.path.join(output_dir, 'HWE_distribution.pdf')
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log(f"Histogram saved: {out}")


##########################################################################
### PLOT 2: QQ PLOT
##########################################################################

def plot_qq(df, output_dir, threshold, mode):

    threshold_log = -np.log10(threshold)
    n             = len(df)

    obs = np.sort(df['log10P'].values)[::-1]
    exp = -np.log10(np.arange(1, n + 1) / n)

    # sample for large datasets
    if n > 50000:
        rng = np.random.default_rng(42)
        idx = np.sort(rng.choice(n, 50000, replace=False))
        obs = obs[idx]
        exp = exp[idx]

    # genomic inflation factor lambda
    lambda_gc = np.median(df['log10P'].values) / (-np.log10(0.5))

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot([0, exp.max()], [0, exp.max()],
            color='#888', linestyle='--', linewidth=0.9, zorder=1)

    ax.scatter(exp, obs,
               c='#1D9E75', s=8, alpha=0.5, zorder=2)

    ax.axhline(threshold_log, color='#ff5f57',
               linestyle='--', linewidth=0.9)
    ax.text(0.02, threshold_log + 0.1,
            f'p = {threshold:.0e}',
            color='#ff5f57', fontsize=8, fontweight='bold')

    ax.set_xlabel(r'Expected $-\log_{10}(p)$',
                  fontweight='bold', fontsize=10)
    ax.set_ylabel(r'Observed $-\log_{10}(p)$',
                  fontweight='bold', fontsize=10)
    ax.set_title('CRISP: Step 8 HWE, QQ Plot',
                 fontweight='bold', fontsize=13)
    apply_style(ax)
    fig.suptitle(
        f'Mode: {mode}  |  lambda = {lambda_gc:.4f}  |  {n:,} variants',
        fontsize=9, color='#555', y=0.97)
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    out = os.path.join(output_dir, 'HWE_qq.pdf')
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log(f"QQ plot saved: {out}")


##########################################################################
### PLOT 3: META-ANALYSIS CLASSIFICATION
##########################################################################

def plot_meta_classification(meta_file, output_dir, threshold, mode):

    df_meta = pd.read_csv(meta_file, sep='\t')

    if len(df_meta) == 0 or 'CLASS' not in df_meta.columns or \
       'P_META' not in df_meta.columns:
        log("WARNING: Meta-analysis file missing expected columns. Skipping plot.")
        return

    df_meta = df_meta[df_meta['P_META'].notna() & (df_meta['P_META'] > 0)].copy()
    df_meta['log10P_meta'] = -np.log10(df_meta['P_META'])

    CLASS_COLORS = {
        'PASS'               : '#1D9E75',
        'FAIL_ALL'           : '#ff5f57',
        'FAIL_CASES_ONLY'    : '#febc2e',
        'FAIL_CONTROLS_ONLY' : '#ff8c00',
    }
    CLASS_LABELS = {
        'PASS'               : 'Pass',
        'FAIL_ALL'           : 'Fail all strata (likely error)',
        'FAIL_CASES_ONLY'    : 'Fail cases only (potential signal)',
        'FAIL_CONTROLS_ONLY' : 'Fail controls only',
    }

    threshold_log = -np.log10(threshold)
    fig, ax = plt.subplots(figsize=(11, 6))

    for cls in ['PASS', 'FAIL_ALL', 'FAIL_CASES_ONLY', 'FAIL_CONTROLS_ONLY']:
        grp = df_meta[df_meta['CLASS'] == cls]
        if len(grp) == 0:
            continue
        ax.hist(grp['log10P_meta'],
                bins=80, color=CLASS_COLORS[cls], alpha=0.85,
                edgecolor='white', linewidth=0.15,
                label=f"{CLASS_LABELS[cls]} ({len(grp):,})")

    ax.axvline(threshold_log, color='#ff5f57',
               linestyle='--', linewidth=1.2)
    ax.text(threshold_log + 0.1, ax.get_ylim()[1] * 0.92,
            f'p = {threshold:.0e}',
            color='#ff5f57', fontsize=8.5, fontweight='bold')

    ax.set_xlabel(r'$-\log_{10}(p_{meta})$', fontweight='bold', fontsize=10)
    ax.set_ylabel('Number of variants', fontweight='bold', fontsize=10)
    ax.set_title("CRISP: Step 8 HWE, Meta-analysis Classification",
                 fontweight='bold', fontsize=13)
    apply_style(ax)
    fig.suptitle(
        f"Fisher's combined p-value across strata  |  threshold p < {threshold:.0e}",
        fontsize=9, color='#555', y=0.97)
    ax.legend(fontsize=9, loc='center right',
              bbox_to_anchor=(1.18, 0.5), framealpha=0.9)
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    out = os.path.join(output_dir, 'HWE_meta_classification.pdf')
    plt.savefig(out, bbox_inches='tight')
    plt.close()
    log(f"Meta-analysis classification plot saved: {out}")


##########################################################################
### MAIN
##########################################################################

def main():
    parser = argparse.ArgumentParser(
        description='CRISP Step 8: HWE plotting'
    )
    parser.add_argument('output_dir')
    parser.add_argument('hwe_file')
    parser.add_argument('threshold', type=float)
    parser.add_argument('mode')
    parser.add_argument('hwe_meta',   nargs='?', default='NO')
    parser.add_argument('meta_file',  nargs='?', default='')
    args = parser.parse_args()

    mode = args.mode.upper()

    log(f"CRISP Step 8: HWE plotting")
    log(f"File      : {args.hwe_file}")
    log(f"Threshold : {args.threshold:.0e}")
    log(f"Mode      : {mode}")

    if not os.path.isfile(args.hwe_file):
        log(f"ERROR: HWE file not found: {args.hwe_file}")
        sys.exit(1)

    df = pd.read_csv(args.hwe_file, sep=r'\s+')

    if len(df) == 0:
        log("ERROR: HWE file is empty.")
        sys.exit(1)

    # keep ALL test rows, drop missing and zero p-values
    df = df[df['TEST'] == 'ALL'].copy()
    df = df[df['P'].notna() & (df['P'] > 0)].copy()
    df['log10P'] = -np.log10(df['P'])

    log(f"Variants with valid p-values: {len(df):,}")
    n_fail = (df['P'] < args.threshold).sum()
    log(f"Variants failing HWE (p < {args.threshold:.0e}): {n_fail:,}")

    os.makedirs(args.output_dir, exist_ok=True)

    plot_histogram(df, args.output_dir, args.threshold, mode)
    plot_qq(df, args.output_dir, args.threshold, mode)

    # meta-analysis classification plot
    if args.hwe_meta.upper() == 'YES' and args.meta_file and \
       os.path.isfile(args.meta_file):
        plot_meta_classification(args.meta_file, args.output_dir,
                                 args.threshold, mode)

    log("HWE plotting complete.")


if __name__ == '__main__':
    main()
