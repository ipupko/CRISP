#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 7: Relatedness, Plotting Script
# Version: 0.1.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# Generates relatedness plots from PLINK IBD or KING output.
#
# IBD mode:
#   Plot 1: IBD0 vs IBD1 scatter coloured by relationship class
#           Theoretical relationship positions overlaid
#   Plot 2: PI_HAT distribution histogram with threshold line
#
# KING mode:
#   Plot 1: Kinship coefficient distribution histogram
#   Plot 2: Kinship vs IBS0 scatter (if IBS0 column available)
#
# Usage:
#   python3 plot_relatedness.py <output_dir> <relatedness_file>
#                               <method> <pihat_cutoff> <king_cutoff>
##########################################################################

import sys
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

CAPTION = 'CRISP | Comprehensive Robust Integrated SNP Processing | Compass Genomics'

REL_COLORS = {
    'Duplicate / MZ twin' : '#ff5f57',
    'First degree'        : '#ff8c00',
    'Second degree'       : '#febc2e',
    'Distant'             : '#1D9E75',
}

REL_ORDER = ['Duplicate / MZ twin', 'First degree', 'Second degree', 'Distant']


def log(msg):
    print(f"[PLOT] {msg}", flush=True)


def apply_style(ax):
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=9)


def classify_ibd(pihat):
    if pihat >= 0.9:    return 'Duplicate / MZ twin'
    elif pihat >= 0.4:  return 'First degree'
    elif pihat >= 0.185:return 'Second degree'
    else:               return 'Distant'


def classify_king(kinship):
    if kinship >= 0.354:  return 'Duplicate / MZ twin'
    elif kinship >= 0.177:return 'First degree'
    elif kinship >= 0.0884:return 'Second degree'
    else:                  return 'Distant'


def legend_patches():
    return [mpatches.Patch(color=REL_COLORS[r], label=r) for r in REL_ORDER
            if r in REL_COLORS]


def place_labels(ax, df_grp, x_col, y_col, iid1_col, iid2_col,
                 rel_label, col):
    """Place IID1:IID2 (relationship) labels with boundary-aware
    collision-resolving curved arrows."""
    if len(df_grp) == 0:
        return

    ax.scatter(df_grp[x_col], df_grp[y_col],
               c=col, s=120, marker='*', zorder=6,
               edgecolors='#333', linewidths=0.6)

    xlim    = ax.get_xlim()
    ylim    = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]
    margin  = 0.01

    candidates = [
        ( x_range*0.10,  y_range*0.10),
        (-x_range*0.10,  y_range*0.10),
        ( x_range*0.10, -y_range*0.10),
        (-x_range*0.10, -y_range*0.10),
        ( x_range*0.16,  y_range*0.05),
        (-x_range*0.16,  y_range*0.05),
        ( x_range*0.05,  y_range*0.16),
        ( x_range*0.05, -y_range*0.16),
        ( x_range*0.20,  y_range*0.08),
    ]

    placed = []

    for _, row in df_grp.iterrows():
        lbl = f"{row[iid1_col]}:{row[iid2_col]} ({rel_label})"
        x0, y0 = row[x_col], row[y_col]
        best = None

        for dx, dy in candidates:
            tx = x0 + dx
            ty = y0 + dy
            if tx < xlim[0] + margin*x_range: continue
            if tx > xlim[1] - margin*x_range: continue
            if ty < ylim[0] + margin*y_range: continue
            if ty > ylim[1] - margin*y_range: continue
            overlap = any(
                abs(tx-px) < x_range*0.09 and abs(ty-py) < y_range*0.07
                for px, py in placed
            )
            if not overlap:
                best = (tx, ty)
                break

        if best is None:
            dx, dy = candidates[0]
            best = (
                float(np.clip(x0+dx, xlim[0]+margin*x_range, xlim[1]-margin*x_range)),
                float(np.clip(y0+dy, ylim[0]+margin*y_range, ylim[1]-margin*y_range))
            )

        placed.append(best)

        ax.annotate(
            lbl,
            xy         = (x0, y0),
            xytext     = best,
            fontsize   = 7.5,
            fontweight = 'bold',
            color      = col,
            bbox       = dict(boxstyle  = 'round,pad=0.3',
                              facecolor = 'white',
                              edgecolor = col,
                              alpha     = 0.90,
                              linewidth = 0.9),
            arrowprops = dict(
                arrowstyle      = '->',
                color           = col,
                lw              = 1.0,
                connectionstyle = 'arc3,rad=0.2'
            ),
            zorder = 7
        )


##########################################################################
### IBD PLOTS
##########################################################################

def plot_ibd(df, output_dir, pihat_cutoff):

    df['RELATIONSHIP'] = df['PI_HAT'].apply(classify_ibd)
    log(f"Pairs loaded : {len(df):,}")
    for rel, n in df['RELATIONSHIP'].value_counts().items():
        log(f"  {rel:<25}: {n:,}")

    # Plot 1: IBD0 vs IBD1 scatter
    fig, ax = plt.subplots(figsize=(10, 7))

    for rel in REL_ORDER:
        grp = df[df['RELATIONSHIP'] == rel]
        if len(grp) == 0:
            continue
        ax.scatter(grp['Z0'], grp['Z1'],
                   c=REL_COLORS[rel], s=18, alpha=0.65,
                   label=rel, zorder=3)

    # theoretical relationship positions
    theory = [
        (0.0,  1.0,   'Parent-offspring'),
        (0.25, 0.5,   'Full siblings'),
        (0.5,  0.5,   'Half siblings'),
        (1.0,  0.0,   'Unrelated'),
    ]
    for tx, ty, tlbl in theory:
        ax.plot(tx, ty, 'x', color='#333', markersize=8, markeredgewidth=1.5,
                zorder=5)
        ax.annotate(tlbl, xy=(tx, ty),
                    xytext=(tx + 0.02, ty + 0.03),
                    fontsize=8, color='#333')

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('IBD0 (proportion sharing 0 alleles)', fontweight='bold', fontsize=10)
    ax.set_ylabel('IBD1 (proportion sharing 1 allele)',  fontweight='bold', fontsize=10)
    ax.set_title('CRISP: Step 7 Relatedness, IBD0 vs IBD1',
                 fontweight='bold', fontsize=13)
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(
        f'PI_HAT cutoff: {pihat_cutoff:.3f}  |  {len(df):,} related pairs',
        fontsize=9, color='#555', y=0.97)
    ax.legend(handles=legend_patches(), loc='center right',
              bbox_to_anchor=(1.18, 0.5),
              fontsize=9, framealpha=0.9)
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    # label duplicates and MZ twins only
    grp_dup = df[df['RELATIONSHIP'] == 'Duplicate / MZ twin'].copy()
    if len(grp_dup) > 0:
        place_labels(ax, grp_dup, 'Z0', 'Z1', 'IID1', 'IID2',
                     'Dup/MZ', REL_COLORS['Duplicate / MZ twin'])

    out1 = os.path.join(output_dir, 'Relatedness_IBD_scatter.pdf')
    plt.savefig(out1, bbox_inches='tight')
    plt.close()
    log(f"IBD scatter saved: {out1}")

    # Plot 2: PI_HAT histogram
    fig, ax = plt.subplots(figsize=(10, 6))

    for rel in REL_ORDER:
        grp = df[df['RELATIONSHIP'] == rel]
        if len(grp) == 0:
            continue
        ax.hist(grp['PI_HAT'], bins=60,
                color=REL_COLORS[rel], alpha=0.82,
                edgecolor='white', linewidth=0.2, label=rel)

    ax.axvline(pihat_cutoff, color='#ff5f57', linestyle='--', linewidth=1.2)
    ax.axvline(0.4,          color='#ff8c00', linestyle=':', linewidth=0.9)
    ax.axvline(0.9,          color='#ff5f57', linestyle=':', linewidth=0.9)
    ax.text(pihat_cutoff + 0.005, ax.get_ylim()[1] * 0.92,
            f'PI_HAT = {pihat_cutoff}',
            color='#ff5f57', fontsize=8, fontweight='bold')

    ax.set_xlabel('PI_HAT (proportion of IBD)', fontweight='bold', fontsize=10)
    ax.set_ylabel('Number of pairs',            fontweight='bold', fontsize=10)
    ax.set_title('CRISP: Step 7 Relatedness, PI_HAT Distribution',
                 fontweight='bold', fontsize=13)
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(
        f'{len(df):,} pairs above PI_HAT >= {pihat_cutoff:.3f}',
        fontsize=9, color='#555', y=0.97)
    ax.legend(handles=legend_patches(), loc='center right',
              bbox_to_anchor=(1.18, 0.5),
              fontsize=9, framealpha=0.9)
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    out2 = os.path.join(output_dir, 'Relatedness_IBD_pihat.pdf')
    plt.savefig(out2, bbox_inches='tight')
    plt.close()
    log(f"PI_HAT histogram saved: {out2}")

    # text report
    report_file = os.path.join(output_dir, 'relatedness_report.txt')
    with open(report_file, 'w') as rep:
        rep.write("=" * 66 + "\n")
        rep.write("  CRISP: STEP 7 RELATEDNESS PAIRS REPORT\n")
        rep.write("  Comprehensive Robust Integrated SNP Processing\n")
        rep.write("=" * 66 + "\n")
        rep.write(f"  Method       : IBD\n")
        rep.write(f"  PI_HAT cutoff: {pihat_cutoff:.3f}\n")
        rep.write(f"  Total pairs  : {len(df):,}\n")
        rep.write("-" * 66 + "\n")
        for rel in ['Duplicate / MZ twin', 'First degree',
                    'Second degree', 'Distant']:
            grp = df[df['RELATIONSHIP'] == rel]
            rep.write(f"\n  {rel} ({len(grp):,} pairs)\n")
            rep.write("-" * 60 + "\n")
            for _, row in grp.iterrows():
                rep.write(
                    f"  FID1={str(row['FID1']):<12} IID1={str(row['IID1']):<12}"
                    f"  FID2={str(row['FID2']):<12} IID2={str(row['IID2']):<12}"
                    f"  PI_HAT={row['PI_HAT']:.4f}"
                    f"  Z0={row['Z0']:.3f}  Z1={row['Z1']:.3f}  Z2={row['Z2']:.3f}\n"
                )
        rep.write("\n" + "=" * 66 + "\n")
        rep.write("  END OF REPORT\n")
        rep.write("=" * 66 + "\n")
    log(f"Relatedness report saved: {report_file}")


##########################################################################
### KING PLOTS
##########################################################################

def plot_king(df, output_dir, king_cutoff):

    if '#FID1' in df.columns:
        df = df.rename(columns={'#FID1': 'FID1'})

    df['RELATIONSHIP'] = df['KINSHIP'].apply(classify_king)
    log(f"Pairs loaded : {len(df):,}")
    for rel, n in df['RELATIONSHIP'].value_counts().items():
        log(f"  {rel:<25}: {n:,}")

    # Plot 1: kinship histogram
    fig, ax = plt.subplots(figsize=(10, 6))

    for rel in REL_ORDER:
        grp = df[df['RELATIONSHIP'] == rel]
        if len(grp) == 0:
            continue
        ax.hist(grp['KINSHIP'], bins=60,
                color=REL_COLORS[rel], alpha=0.82,
                edgecolor='white', linewidth=0.2, label=rel)

    ax.axvline(king_cutoff, color='#ff5f57', linestyle='--', linewidth=1.2)
    ax.axvline(0.177,       color='#ff8c00', linestyle=':', linewidth=0.9)
    ax.axvline(0.354,       color='#ff5f57', linestyle=':', linewidth=0.9)
    ax.text(king_cutoff + 0.003, ax.get_ylim()[1] * 0.92,
            f'KING = {king_cutoff}',
            color='#ff5f57', fontsize=8, fontweight='bold')

    ax.set_xlabel('Kinship coefficient',  fontweight='bold', fontsize=10)
    ax.set_ylabel('Number of pairs',      fontweight='bold', fontsize=10)
    ax.set_title('CRISP: Step 7 Relatedness, KING Kinship Distribution',
                 fontweight='bold', fontsize=13)
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle(
        f'{len(df):,} pairs above KING >= {king_cutoff:.4f}',
        fontsize=9, color='#555', y=0.97)
    ax.legend(handles=legend_patches(), loc='center right',
              bbox_to_anchor=(1.18, 0.5),
              fontsize=9, framealpha=0.9)
    fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])

    out1 = os.path.join(output_dir, 'Relatedness_KING_kinship.pdf')
    plt.savefig(out1, bbox_inches='tight')
    plt.close()
    log(f"KING kinship histogram saved: {out1}")

    # Plot 2: Kinship vs IBS0 if available
    if 'IBS0' in df.columns:

        fig, ax = plt.subplots(figsize=(10, 7))

        for rel in REL_ORDER:
            grp = df[df['RELATIONSHIP'] == rel]
            if len(grp) == 0:
                continue
            ax.scatter(grp['IBS0'], grp['KINSHIP'],
                       c=REL_COLORS[rel], s=18, alpha=0.65,
                       label=rel, zorder=3)

        for thresh, lbl in [(0.0884, '2nd degree'), (0.177, '1st degree'), (0.354, 'Dup')]:
            ax.axhline(thresh, color='#888', linestyle=':', linewidth=0.7)
            ax.text(ax.get_xlim()[1] * 0.98, thresh + 0.005,
                    lbl, ha='right', fontsize=7, color='#888')

        ax.set_xlabel('IBS0 (proportion sharing 0 alleles)', fontweight='bold', fontsize=10)
        ax.set_ylabel('Kinship coefficient',                  fontweight='bold', fontsize=10)
        ax.set_title('CRISP: Step 7 Relatedness, KING Kinship vs IBS0',
                     fontweight='bold', fontsize=13)
        ax.set_facecolor('#f8f9fa')
        ax.spines[['top', 'right']].set_visible(False)
        fig.suptitle(
            f'{len(df):,} related pairs  |  KING cutoff: {king_cutoff:.4f}',
            fontsize=9, color='#555', y=0.97)
        ax.legend(handles=legend_patches(), loc='center right',
              bbox_to_anchor=(1.18, 0.5),
                  fontsize=9, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888')
        plt.tight_layout(rect=[0, 0.02, 1, 0.97])

        # label duplicates and MZ twins only
        grp_dup = df[df['RELATIONSHIP'] == 'Duplicate / MZ twin'].copy()
        if len(grp_dup) > 0:
            place_labels(ax, grp_dup, 'IBS0', 'KINSHIP', 'IID1', 'IID2',
                         'Dup/MZ', REL_COLORS['Duplicate / MZ twin'])

        out2 = os.path.join(output_dir, 'Relatedness_KING_scatter.pdf')
        plt.savefig(out2, bbox_inches='tight')
        plt.close()
        log(f"KING scatter saved: {out2}")


##########################################################################
### MAIN
##########################################################################

def main():
    parser = argparse.ArgumentParser(
        description='CRISP Step 7: Relatedness plotting'
    )
    parser.add_argument('output_dir')
    parser.add_argument('relatedness_file')
    parser.add_argument('method')
    parser.add_argument('pihat_cutoff', type=float)
    parser.add_argument('king_cutoff',  type=float)
    args = parser.parse_args()

    method = args.method.upper()

    log(f"CRISP Step 7: Relatedness plotting")
    log(f"Method  : {method}")
    log(f"File    : {args.relatedness_file}")

    if not os.path.isfile(args.relatedness_file):
        log(f"ERROR: Relatedness file not found: {args.relatedness_file}")
        sys.exit(1)

    df = pd.read_csv(args.relatedness_file, sep=r'\s+')

    if len(df) == 0:
        log("No pairs above threshold. Skipping plots.")
        sys.exit(0)

    os.makedirs(args.output_dir, exist_ok=True)

    if method == 'IBD':
        plot_ibd(df, args.output_dir, args.pihat_cutoff)
    elif method == 'KING':
        plot_king(df, args.output_dir, args.king_cutoff)
    else:
        log(f"ERROR: Unknown method '{method}'. Use IBD or KING.")
        sys.exit(1)

    log("Relatedness plotting complete.")


if __name__ == '__main__':
    main()
