#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 6: Heterozygosity and Homozygosity, Plotting Script
# Version: 0.2.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
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
import matplotlib.ticker as ticker


# ─────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────
COLORS = {
    'Normal'              : '#1D9E75',
    'High Heterozygosity' : '#febc2e',
    'Excess Homozygosity' : '#ff5f57',
}

CAPTION = "CRISP | Comprehensive Robust Integrated SNP Processing"


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
def log(msg: str):
    print(f"[HOMO] {msg}", flush=True)

def abort(msg: str):
    print(f"[HOMO] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# LOAD FILES
# ─────────────────────────────────────────────
def load_file(filepath: str, label: str) -> pd.DataFrame:
    if not os.path.isfile(filepath):
        abort(f"{label} not found: {filepath}")
    df = pd.read_csv(filepath, sep=r'\s+')
    log(f"Loaded {label}: {len(df):,} rows")
    return df


# ─────────────────────────────────────────────
# PREPARE DATA
# ─────────────────────────────────────────────
def prepare_data(df_sex: pd.DataFrame,
                 df_het: pd.DataFrame,
                 df_hom: pd.DataFrame) -> pd.DataFrame:
    """
    Merge sexcheck, het and hom.indiv into a single DataFrame.
    Calculates Z-scores for autosomal and X chromosome F-statistics.
    """
    # Sex index from sexcheck
    sex_index = df_sex[['IID', 'PEDSEX']].copy()

    # Merge sex into het
    df_het = df_het.merge(sex_index, on='IID', how='left')

    # Calculate autosomal Z-score
    f_mean = df_het['F'].mean()
    f_sd   = df_het['F'].std()
    df_het['Z_Score'] = (df_het['F'] - f_mean) / f_sd

    log(f"Autosomal F-stat mean : {f_mean:.4f}")
    log(f"Autosomal F-stat SD   : {f_sd:.4f}")

    # Merge with ROH data
    df_merged = df_het.merge(
        df_hom[['IID', 'KB', 'NSEG']],
        on  = 'IID',
        how = 'left'
    )

    return df_merged, f_mean, f_sd


# ─────────────────────────────────────────────
# CLASSIFY OUTLIERS
# ─────────────────────────────────────────────
def classify(df: pd.DataFrame,
             z_high: float,
             z_low: float) -> pd.DataFrame:
    df = df.copy()
    df['Category'] = 'Normal'
    df.loc[df['Z_Score'] > z_high, 'Category'] = 'Excess Homozygosity'
    df.loc[df['Z_Score'] < z_low,  'Category'] = 'High Heterozygosity'
    return df


# ─────────────────────────────────────────────
# PLOT
# ─────────────────────────────────────────────
def generate_plot(df: pd.DataFrame,
                  z_high: float,
                  z_low: float,
                  output_dir: str) -> str:

    n_excess = (df['Category'] == 'Excess Homozygosity').sum()
    n_het    = (df['Category'] == 'High Heterozygosity').sum()
    n_normal = (df['Category'] == 'Normal').sum()

    fig, ax = plt.subplots(figsize=(11, 7))

    for cat, grp in df.groupby('Category'):
        ax.scatter(
            grp['KB'], grp['Z_Score'],
            c      = COLORS[cat],
            s      = 18,
            alpha  = 0.75,
            label  = cat,
            zorder = 2
        )

    # Threshold lines
    ax.axhline(0,      color='black',   linestyle='dashed',
               linewidth=0.6, alpha=0.5)
    ax.axhline(z_high, color='#ff5f57', linestyle='dashed',
               linewidth=1.0, alpha=0.8)
    ax.axhline(z_low,  color='#febc2e', linestyle='dashed',
               linewidth=1.0, alpha=0.8)

    # Threshold annotations
    xlim = ax.get_xlim()
    xpos = xlim[0] + (xlim[1] - xlim[0]) * 0.02

    ax.text(xpos, z_high + 0.12,
            f'Z = {z_high} (excess homozygosity)',
            color='#ff5f57', fontsize=8, fontweight='bold')
    ax.text(xpos, z_low - 0.18,
            f'Z = {z_low} (high heterozygosity)',
            color='#febc2e', fontsize=8, fontweight='bold')

    # Axes formatting
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, _: f'{x:,.0f}'))
    ax.set_xlabel('Total runs of homozygosity (KB)',
                  fontweight='bold', fontsize=10)
    ax.set_ylabel('Homozygosity Z-score',
                  fontweight='bold', fontsize=10)
    ax.set_title(
        'CRISP: Step 6: Total Homozygosity Runs vs Homozygosity Z-score',
        fontweight='bold', fontsize=12
    )
    fig.suptitle(
        f'Z > {z_high}: {n_excess} excess homozygosity  |  '
        f'Z < {z_low}: {n_het} high heterozygosity  |  '
        f'{n_normal} normal',
        fontsize=9, color='#555555', y=0.96
    )

    # Style
    ax.set_facecolor('#f8f9fa')
    fig.patch.set_facecolor('white')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=9)

    # Legend
    patches = [mpatches.Patch(color=COLORS[c], label=c)
               for c in COLORS]
    ax.legend(handles=patches, loc='upper right',
              fontsize=9, framealpha=0.9)

    fig.text(0.5, 0.01, CAPTION,
             ha='center', fontsize=8, color='#888888')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    out_file = os.path.join(
        output_dir, "Homozygosity_Runs_vs_Zscore.pdf"
    )
    fig.savefig(out_file, bbox_inches='tight', dpi=150)
    plt.close(fig)

    log(f"Plot saved: {out_file}")
    return out_file


# ─────────────────────────────────────────────
# WRITE REPORT
# ─────────────────────────────────────────────
def write_report(df: pd.DataFrame,
                 z_high: float,
                 z_low: float,
                 f_mean: float,
                 f_sd: float,
                 output_dir: str) -> str:

    outliers_excess = df[df['Category'] == 'Excess Homozygosity']
    outliers_het    = df[df['Category'] == 'High Heterozygosity']

    report_path = os.path.join(output_dir, "report_homozygosity.txt")

    with open(report_path, 'w') as out:
        out.write("=" * 64 + "\n")
        out.write("  CRISP: STEP 6 HETEROZYGOSITY AND HOMOZYGOSITY REPORT\n")
        out.write("  Comprehensive Robust Integrated SNP Processing\n")
        out.write("=" * 64 + "\n\n")
        out.write(f"  Z-score upper : {z_high} (excess homozygosity)\n")
        out.write(f"  Z-score lower : {z_low} (high heterozygosity)\n\n")

        out.write("-" * 64 + "\n")
        out.write("AUTOSOMAL HETEROZYGOSITY\n")
        out.write("-" * 64 + "\n")
        out.write(f"  Samples total          : {len(df):,}\n")
        out.write(f"  F-stat mean            : {f_mean:.4f}\n")
        out.write(f"  F-stat SD              : {f_sd:.4f}\n")
        out.write(f"  Excess homozygosity    : {len(outliers_excess):,} "
                  f"(Z > {z_high})\n")
        out.write(f"  High heterozygosity    : {len(outliers_het):,} "
                  f"(Z < {z_low})\n\n")

        out.write("-" * 64 + "\n")
        out.write("EXCESS HOMOZYGOSITY SAMPLES\n")
        out.write("-" * 64 + "\n")
        if len(outliers_excess) > 0:
            out.write(outliers_excess[
                ['FID', 'IID', 'F', 'Z_Score', 'KB']
            ].to_string(index=False))
            out.write("\n\n")
        else:
            out.write("  None detected.\n\n")

        out.write("-" * 64 + "\n")
        out.write("HIGH HETEROZYGOSITY SAMPLES\n")
        out.write("-" * 64 + "\n")
        if len(outliers_het) > 0:
            out.write(outliers_het[
                ['FID', 'IID', 'F', 'Z_Score', 'KB']
            ].to_string(index=False))
            out.write("\n\n")
        else:
            out.write("  None detected.\n\n")

        out.write("=" * 64 + "\n")
        out.write("  END OF REPORT\n")
        out.write("=" * 64 + "\n")

    log(f"Report written: {report_path}")
    return report_path


# ─────────────────────────────────────────────
# WRITE EXCLUSION LISTS
# ─────────────────────────────────────────────
def write_exclusion_lists(df: pd.DataFrame, output_dir: str):

    outliers_excess = df[df['Category'] == 'Excess Homozygosity']
    outliers_het    = df[df['Category'] == 'High Heterozygosity']

    # High heterozygosity
    het_file = os.path.join(output_dir, "outliers_high_het.txt")
    outliers_het[['FID', 'IID']].to_csv(
        het_file, sep='\t', index=False, header=False
    )
    log(f"High het exclusion list     : {het_file} "
        f"({len(outliers_het):,} samples)")

    # Excess homozygosity
    homo_file = os.path.join(output_dir, "outliers_excess_homo.txt")
    outliers_excess[['FID', 'IID']].to_csv(
        homo_file, sep='\t', index=False, header=False
    )
    log(f"Excess homo exclusion list  : {homo_file} "
        f"({len(outliers_excess):,} samples)")

    # Combined
    combined = pd.concat(
        [outliers_het[['FID', 'IID']],
         outliers_excess[['FID', 'IID']]]
    ).drop_duplicates()

    combined_file = os.path.join(output_dir, "outliers_combined.txt")
    combined.to_csv(
        combined_file, sep='\t', index=False, header=False
    )
    log(f"Combined exclusion list     : {combined_file} "
        f"({len(combined):,} samples)")

    return het_file, homo_file, combined_file


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRISP Step 6: Heterozygosity and Homozygosity (Python)"
    )
    parser.add_argument("output_dir")
    parser.add_argument("sexcheck_file")
    parser.add_argument("het_file")
    parser.add_argument("hom_indiv_file")
    parser.add_argument("z_high", type=float)
    parser.add_argument("z_low",  type=float)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    log("CRISP Step 6: Heterozygosity and Homozygosity (Python)")
    log(f"Z-score upper threshold : {args.z_high} (excess homozygosity)")
    log(f"Z-score lower threshold : {args.z_low} (high heterozygosity)")
    log("")

    # Load files
    df_sex = load_file(args.sexcheck_file, ".sexcheck")
    df_het = load_file(args.het_file,      ".het")
    df_hom = load_file(args.hom_indiv_file, ".hom.indiv")
    log("")

    # Prepare and merge
    log("Preparing data...")
    df_merged, f_mean, f_sd = prepare_data(df_sex, df_het, df_hom)
    log("")

    # Classify outliers
    log("Classifying outliers...")
    df_classified = classify(df_merged, args.z_high, args.z_low)

    n_excess = (df_classified['Category'] == 'Excess Homozygosity').sum()
    n_het    = (df_classified['Category'] == 'High Heterozygosity').sum()
    n_normal = (df_classified['Category'] == 'Normal').sum()

    log(f"  Excess homozygosity (Z > {args.z_high}) : {n_excess:,}")
    log(f"  High heterozygosity (Z < {args.z_low}) : {n_het:,}")
    log(f"  Normal                                  : {n_normal:,}")
    log("")

    # Generate plot
    log("Generating plot...")
    generate_plot(df_classified, args.z_high, args.z_low,
                  args.output_dir)
    log("")

    # Write report
    log("Writing report...")
    write_report(df_classified, args.z_high, args.z_low,
                 f_mean, f_sd, args.output_dir)
    log("")

    # Write exclusion lists
    log("Writing exclusion lists...")
    write_exclusion_lists(df_classified, args.output_dir)

    log("\nHomozygosity check complete.")


if __name__ == "__main__":
    main()
