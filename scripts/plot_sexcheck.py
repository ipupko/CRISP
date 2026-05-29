#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 5: Sex Check, Plotting Script
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
import matplotlib.lines as mlines


# ─────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────
COLORS = {
    1 : '#e05c4b',   # male
    2 : '#1D9E75',   # female
    0 : '#aaaaaa',   # unknown
}

ANOMALY_COLORS = {
    'mismatch_f'  : '#ff8c00',
    'mismatch_m'  : '#9400d3',
    'turner'      : '#0000ff',
    'klinefelter' : '#ff00ff',
    'xxx'         : '#00ced1',
}

CAPTION = "CRISP | Comprehensive Robust Integrated SNP Processing"


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
def log(msg: str):
    print(f"[SEXCHECK] {msg}", flush=True)

def abort(msg: str):
    print(f"[SEXCHECK] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────
def load_data(filepath: str) -> pd.DataFrame:
    if not os.path.isfile(filepath):
        abort(f"Input file not found: {filepath}")
    df = pd.read_csv(filepath, sep=r'\s+')
    for col in ['FID', 'IID', 'PEDSEX', 'F', 'YCOUNT']:
        if col not in df.columns:
            abort(f"Required column '{col}' not found in: {filepath}")
    return df


# ─────────────────────────────────────────────
# SHARED PLOT STYLE
# ─────────────────────────────────────────────
def apply_style(ax):
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.tick_params(labelsize=9)


def add_threshold_line(ax, x=None, y=None, color='grey',
                       label='', linestyle='--'):
    if x is not None:
        ax.axvline(x, color=color, linestyle=linestyle,
                   linewidth=1.2, alpha=0.8)
        if label:
            ylim = ax.get_ylim()
            ax.text(x + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.01,
                    ylim[1] * 0.95, label,
                    color=color, fontsize=8, fontweight='bold')
    if y is not None:
        ax.axhline(y, color=color, linestyle=':', linewidth=1.0, alpha=0.7)
        if label:
            ax.text(ax.get_xlim()[0], y + 2, label,
                    color=color, fontsize=8)


# ─────────────────────────────────────────────
# SCATTER PLOT HELPER
# ─────────────────────────────────────────────
def draw_scatter(ax, df_plot: pd.DataFrame,
                 anomaly_masks: dict,
                 title: str, y_threshold: float):
    """
    Draw F-stat vs Y-count scatter on ax.
    anomaly_masks: dict of label -> boolean mask Series
    """
    # Plot normal samples
    for sex in [1, 2, 0]:
        grp = df_plot[df_plot['PEDSEX'] == sex]
        if len(grp) == 0:
            continue
        label = {1: 'Male', 2: 'Female', 0: 'Unknown'}[sex]
        ax.scatter(grp['F'], grp['YCOUNT'],
                   c=COLORS[sex], s=12, alpha=0.65,
                   label=label, zorder=2)

    # Overlay anomalies
    for atype, mask in anomaly_masks.items():
        grp = df_plot[mask]
        if len(grp) == 0:
            continue
        ax.scatter(grp['F'], grp['YCOUNT'],
                   c=ANOMALY_COLORS.get(atype, '#ffcc00'),
                   s=150, marker='*', zorder=5,
                   edgecolors='#333333', linewidths=0.5,
                   label=atype.replace('_', ' ').title())

    # Y threshold line
    ax.axhline(y_threshold, color='#666666',
               linestyle=':', linewidth=1.0, alpha=0.7)
    ax.text(ax.get_xlim()[0],
            y_threshold + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02,
            f'Y mean = {y_threshold:.1f}',
            color='#666666', fontsize=8)

    ax.set_xlabel('F Statistic', fontweight='bold', fontsize=10)
    ax.set_ylabel('Y Count',     fontweight='bold', fontsize=10)
    ax.set_title(title, fontweight='bold', fontsize=12, pad=8)
    apply_style(ax)


# ─────────────────────────────────────────────
# DETECTION
# ─────────────────────────────────────────────
def detect_anomalies(df: pd.DataFrame,
                     f_female_max: float,
                     f_male_min: float,
                     f_turner: float,
                     f_klinefelter: float,
                     f_xxx: float,
                     y_threshold: float) -> dict:
    """
    Returns dict of anomaly label -> DataFrame of flagged samples.
    """
    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]

    results = {
        'sex_mismatch_females' : df_femal[
            (df_femal['F'] > f_male_min) &
            (df_femal['YCOUNT'] > y_threshold)
        ],
        'sex_mismatch_males'   : df_mal[
            (df_mal['F'] < f_female_max) &
            (df_mal['YCOUNT'] < y_threshold)
        ],
        'turner_x0'            : df_femal[
            (df_femal['F'] > f_turner) &
            (df_femal['YCOUNT'] < y_threshold)
        ],
        'klinefelter_xxy'      : df_mal[
            (df_mal['F'] < f_klinefelter) &
            (df_mal['YCOUNT'] > y_threshold)
        ],
        'triple_x'             : df_femal[
            (df_femal['F'] < f_xxx) &
            (df_femal['YCOUNT'] < y_threshold)
        ],
    }
    return results


# ─────────────────────────────────────────────
# WRITE REPORT
# ─────────────────────────────────────────────
def write_report(df: pd.DataFrame,
                 anomalies: dict,
                 y_threshold: float,
                 y_use_mean: str,
                 params: dict,
                 output_dir: str):

    report_path = os.path.join(
        output_dir, "report_sex.mismatch_aneuploidies.txt"
    )

    with open(report_path, 'w') as out:
        out.write("=" * 64 + "\n")
        out.write("  CRISP: STEP 5 SEX CHECK AND ANEUPLOIDY REPORT\n")
        out.write("  Comprehensive Robust Integrated SNP Processing\n")
        out.write("=" * 64 + "\n\n")
        out.write(f"  F female max   : {params['f_female_max']}\n")
        out.write(f"  F male min     : {params['f_male_min']}\n")
        out.write(f"  Y threshold    : {y_threshold:.2f} "
                  f"({'mean' if y_use_mean == 'YES' else 'manual'})\n\n")

        df_mal   = df[df['PEDSEX'] == 1]
        df_femal = df[df['PEDSEX'] == 2]

        out.write("-" * 64 + "\n")
        out.write("SAMPLE COUNTS\n")
        out.write("-" * 64 + "\n")
        out.write(f"  Reported males   : {len(df_mal)}\n")
        out.write(f"  Reported females : {len(df_femal)}\n\n")

        out.write("-" * 64 + "\n")
        out.write("SEX MISMATCHES\n")
        out.write("-" * 64 + "\n")
        for key in ['sex_mismatch_females', 'sex_mismatch_males']:
            grp = anomalies[key]
            label = key.replace('_', ' ').title()
            out.write(f"  {label} : {len(grp)}\n")
            if len(grp) > 0:
                out.write("\n  Details:\n")
                out.write(grp.to_string(index=False))
                out.write("\n\n")

        out.write("-" * 64 + "\n")
        out.write("CHROMOSOMAL ANEUPLOIDIES\n")
        out.write("-" * 64 + "\n")
        labels = {
            'turner_x0'       : 'Turner syndrome (X0)',
            'klinefelter_xxy' : 'Klinefelter syndrome (XXY)',
            'triple_x'        : 'Triple-X syndrome (XXX)',
        }
        for key, label in labels.items():
            grp = anomalies[key]
            out.write(f"  {label} : {len(grp)}\n")
            if len(grp) > 0:
                out.write("\n  Details:\n")
                out.write(grp.to_string(index=False))
                out.write("\n\n")

        # Combined aneuploidies
        df_aneu = pd.concat([
            anomalies['turner_x0'],
            anomalies['klinefelter_xxy'],
            anomalies['triple_x']
        ]).drop_duplicates()

        out.write("-" * 64 + "\n")
        out.write(f"  Total aneuploidies flagged : {len(df_aneu)}\n")
        if len(df_aneu) > 0:
            out.write("\n  All aneuploidy samples:\n")
            out.write(df_aneu.to_string(index=False))
            out.write("\n")

        out.write("\n" + "=" * 64 + "\n")
        out.write("  END OF REPORT\n")
        out.write("=" * 64 + "\n")

    log(f"Report written: {report_path}")
    return report_path


# ─────────────────────────────────────────────
# WRITE EXCLUSION LISTS
# ─────────────────────────────────────────────
def write_exclusion_lists(anomalies: dict, output_dir: str):

    # Sex mismatch list
    df_mismatch = pd.concat([
        anomalies['sex_mismatch_females'],
        anomalies['sex_mismatch_males']
    ]).drop_duplicates()

    mismatch_file = os.path.join(output_dir, "id_list_sex_mismatch.txt")
    df_mismatch[['FID', 'IID']].to_csv(
        mismatch_file, sep='\t', index=False, header=False
    )
    log(f"Sex mismatch exclusion list: {mismatch_file} "
        f"({len(df_mismatch)} samples)")

    # Aneuploidy list
    df_aneu = pd.concat([
        anomalies['turner_x0'],
        anomalies['klinefelter_xxy'],
        anomalies['triple_x']
    ]).drop_duplicates()

    aneuploidy_file = os.path.join(output_dir, "id_list_aneuploidies.txt")
    df_aneu[['FID', 'IID']].to_csv(
        aneuploidy_file, sep='\t', index=False, header=False
    )
    log(f"Aneuploidy exclusion list  : {aneuploidy_file} "
        f"({len(df_aneu)} samples)")

    return mismatch_file, aneuploidy_file


# ─────────────────────────────────────────────
# GENERATE PLOTS
# ─────────────────────────────────────────────
def generate_plots(df: pd.DataFrame,
                   anomalies: dict,
                   params: dict,
                   y_threshold: float,
                   output_dir: str):

    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]

    pdf_file = os.path.join(output_dir, "Sex_check.pdf")

    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(pdf_file) as pdf:

        # ── OVERALL ──────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 7))

        anomaly_masks_overall = {
            'mismatch_f'  : df.index.isin(anomalies['sex_mismatch_females'].index),
            'mismatch_m'  : df.index.isin(anomalies['sex_mismatch_males'].index),
            'turner'      : df.index.isin(anomalies['turner_x0'].index),
            'klinefelter' : df.index.isin(anomalies['klinefelter_xxy'].index),
            'xxx'         : df.index.isin(anomalies['triple_x'].index),
        }

        draw_scatter(ax, df, anomaly_masks_overall,
                     "Overall Sex Check", y_threshold)
        add_threshold_line(ax, x=params['f_female_max'],
                           color='#1D9E75',
                           label=f"F={params['f_female_max']}")
        add_threshold_line(ax, x=params['f_male_min'],
                           color='#e05c4b',
                           label=f"F={params['f_male_min']}")

        n_aneu = sum(len(v) for k,v in anomalies.items()
                     if k not in ('sex_mismatch_females','sex_mismatch_males'))
        n_mis  = (len(anomalies['sex_mismatch_females']) +
                  len(anomalies['sex_mismatch_males']))

        fig.suptitle(
            f"500 samples  |  {n_mis} sex mismatches  |  "
            f"{n_aneu} aneuploidies flagged",
            fontsize=9, color='#555555', y=0.97
        )
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── MALE ─────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 6))
        anomaly_masks_male = {
            'mismatch_m'  : df_mal.index.isin(anomalies['sex_mismatch_males'].index),
            'klinefelter' : df_mal.index.isin(anomalies['klinefelter_xxy'].index),
        }
        draw_scatter(ax, df_mal, anomaly_masks_male,
                     "Male Sex Check", y_threshold)
        add_threshold_line(ax, x=params['f_klinefelter'],
                           color='darkred',
                           label=f"Klinefelter F={params['f_klinefelter']}")
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 1.0])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── FEMALE ───────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 6))
        anomaly_masks_female = {
            'mismatch_f' : df_femal.index.isin(
                anomalies['sex_mismatch_females'].index),
            'turner'     : df_femal.index.isin(anomalies['turner_x0'].index),
            'xxx'        : df_femal.index.isin(anomalies['triple_x'].index),
        }
        draw_scatter(ax, df_femal, anomaly_masks_female,
                     "Female Sex Check", y_threshold)
        add_threshold_line(ax, x=params['f_turner'],
                           color='purple',
                           label=f"Turner F={params['f_turner']}")
        add_threshold_line(ax, x=params['f_xxx'],
                           color='blue',
                           label=f"XXX F={params['f_xxx']}")
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 1.0])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    log(f"Plots saved: {pdf_file} (3 pages: overall, male, female)")
    return pdf_file


# ─────────────────────────────────────────────
# SEX AMENDMENT
# ─────────────────────────────────────────────
def apply_sex_amendment(df: pd.DataFrame, output_dir: str):
    """
    Updates PEDSEX to match SNPSEX for mismatched samples.
    Writes an amended FAM file and a detailed amendment log.
    Original FAM file is never overwritten.
    """
    if 'SNPSEX' not in df.columns:
        log("WARNING: SNPSEX column not found. Cannot apply sex amendment.")
        return

    mismatches = df[
        (df['PEDSEX'] != df['SNPSEX']) &
        (df['SNPSEX'].isin([1, 2]))
    ][['FID', 'IID', 'PEDSEX', 'SNPSEX']]

    if len(mismatches) == 0:
        log("No sex mismatches with valid SNPSEX found. No amendments made.")
        return

    amendment_log = os.path.join(output_dir, "sex_amendments.txt")

    sex_labels = {1: 'Male', 2: 'Female', 0: 'Unknown'}

    with open(amendment_log, 'w') as out:
        out.write("=" * 64 + "\n")
        out.write("  CRISP: SEX AMENDMENT LOG\n")
        out.write("  Comprehensive Robust Integrated SNP Processing\n")
        out.write("=" * 64 + "\n\n")
        out.write("  WARNING: These samples had their recorded sex updated\n")
        out.write("  to match the genotype-inferred sex (SNPSEX).\n")
        out.write("  Consult the data provider before proceeding.\n\n")
        out.write(f"  {'FID':<12} {'IID':<12} {'Original':<12} {'Amended':<12}\n")
        out.write("-" * 64 + "\n")

        for _, row in mismatches.iterrows():
            orig = sex_labels.get(int(row['PEDSEX']), 'Unknown')
            new  = sex_labels.get(int(row['SNPSEX']), 'Unknown')
            out.write(
                f"  {row['FID']:<12} {row['IID']:<12} "
                f"{orig:<12} {new:<12}\n"
            )
            log(f"WARNING: Amending {row['FID']} {row['IID']} "
                f"SEX {int(row['PEDSEX'])} -> {int(row['SNPSEX'])}")

        out.write("\n" + "=" * 64 + "\n")

    log(f"Amendment log written : {amendment_log}")
    log(f"Samples amended       : {len(mismatches)}")
    log("")
    log("!" * 60)
    log("Please review the amendment log before continuing.")
    log("!" * 60)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRISP Step 5: Sex Check and Aneuploidy Detection (Python)"
    )
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    parser.add_argument("f_stat_y_count_file")
    parser.add_argument("f_female_max",          type=float)
    parser.add_argument("f_male_min",             type=float)
    parser.add_argument("f_turner",               type=float)
    parser.add_argument("f_klinefelter",          type=float)
    parser.add_argument("f_xxx",                  type=float)
    parser.add_argument("y_use_mean")
    parser.add_argument("y_manual",               type=float)
    parser.add_argument("--exclude_mismatch",     default="YES")
    parser.add_argument("--exclude_aneuploidy",   default="YES")
    parser.add_argument("--amend_sex",            default="NO")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    log("CRISP Step 5: Sex Check and Aneuploidy Detection (Python)")
    log(f"Parameters:")
    log(f"  F female max              : {args.f_female_max}")
    log(f"  F male min                : {args.f_male_min}")
    log(f"  F Turner                  : {args.f_turner}")
    log(f"  F Klinefelter             : {args.f_klinefelter}")
    log(f"  F XXX threshold           : {args.f_xxx}")
    log(f"  Y use mean                : {args.y_use_mean}")
    log(f"  Exclude mismatches        : {args.exclude_mismatch}")
    log(f"  Exclude aneuploidies      : {args.exclude_aneuploidy}")
    log(f"  Amend sex in FAM          : {args.amend_sex}")
    if args.y_use_mean.upper() != "YES":
        log(f"  Y manual cutoff           : {args.y_manual}")
    log("")

    # Mutual exclusion check
    if (args.amend_sex.upper() == "YES" and
            args.exclude_mismatch.upper() == "YES"):
        abort("SEXCHECK_AMEND_SEX = YES and SEXCHECK_EXCLUDE_MISMATCH = YES "
              "cannot both be active. Choose one approach.")

    # Load data
    log(f"Loading: {args.f_stat_y_count_file}")
    df = load_data(args.f_stat_y_count_file)
    log(f"Samples loaded: {len(df)}")

    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]
    log(f"Reported males   : {len(df_mal)}")
    log(f"Reported females : {len(df_femal)}")
    log("")

    # Y threshold
    if args.y_use_mean.upper() == "YES":
        y_threshold = float(df['YCOUNT'].mean())
        log(f"Y count threshold: mean = {y_threshold:.2f}")
    else:
        y_threshold = args.y_manual
        log(f"Y count threshold: manual = {y_threshold:.2f}")
    log("")

    params = {
        'f_female_max'  : args.f_female_max,
        'f_male_min'    : args.f_male_min,
        'f_turner'      : args.f_turner,
        'f_klinefelter' : args.f_klinefelter,
        'f_xxx'         : args.f_xxx,
    }

    # Detect anomalies
    log("Running detection...")
    anomalies = detect_anomalies(
        df, args.f_female_max, args.f_male_min,
        args.f_turner, args.f_klinefelter, args.f_xxx,
        y_threshold
    )

    log(f"  Sex mismatch (females)    : {len(anomalies['sex_mismatch_females'])}")
    log(f"  Sex mismatch (males)      : {len(anomalies['sex_mismatch_males'])}")
    log(f"  Turner syndrome (X0)      : {len(anomalies['turner_x0'])}")
    log(f"  Klinefelter syndrome (XXY): {len(anomalies['klinefelter_xxy'])}")
    log(f"  Triple-X syndrome (XXX)   : {len(anomalies['triple_x'])}")
    log("")

    # Generate plots
    log("Generating plots...")
    generate_plots(df, anomalies, params, y_threshold, args.output_dir)

    # Write detection report
    log("Writing detection report...")
    write_report(df, anomalies, y_threshold,
                 args.y_use_mean, params, args.output_dir)

    # Write exclusion lists
    log("Writing exclusion lists...")
    write_exclusion_lists(anomalies, args.output_dir)

    # Sex amendment
    if args.amend_sex.upper() == "YES":
        log("")
        log("!" * 60)
        log("SEX AMENDMENT REQUESTED.")
        log("Updating PEDSEX in FAM file to match genotype-inferred sex.")
        log("It is strongly recommended to consult the data provider")
        log("for each amended sample before proceeding.")
        log("!" * 60)
        log("")
        apply_sex_amendment(df, args.output_dir)

    log("\nSex check complete.")


if __name__ == "__main__":
    main()ANOMALY_COLORS = {
    'mismatch_f'  : '#ff8c00',
    'mismatch_m'  : '#9400d3',
    'turner'      : '#0000ff',
    'klinefelter' : '#ff00ff',
    'xxx'         : '#00ced1',
}

CAPTION = "CRISP | Comprehensive Robust Integrated SNP Processing"


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
def log(msg: str):
    print(f"[SEXCHECK] {msg}", flush=True)

def abort(msg: str):
    print(f"[SEXCHECK] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────
def load_data(filepath: str) -> pd.DataFrame:
    if not os.path.isfile(filepath):
        abort(f"Input file not found: {filepath}")
    df = pd.read_csv(filepath, sep=r'\s+')
    for col in ['FID', 'IID', 'PEDSEX', 'F', 'YCOUNT']:
        if col not in df.columns:
            abort(f"Required column '{col}' not found in: {filepath}")
    return df


# ─────────────────────────────────────────────
# SHARED PLOT STYLE
# ─────────────────────────────────────────────
def apply_style(ax):
    ax.set_facecolor('#f8f9fa')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.tick_params(labelsize=9)


def add_threshold_line(ax, x=None, y=None, color='grey',
                       label='', linestyle='--'):
    if x is not None:
        ax.axvline(x, color=color, linestyle=linestyle,
                   linewidth=1.2, alpha=0.8)
        if label:
            ylim = ax.get_ylim()
            ax.text(x + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.01,
                    ylim[1] * 0.95, label,
                    color=color, fontsize=8, fontweight='bold')
    if y is not None:
        ax.axhline(y, color=color, linestyle=':', linewidth=1.0, alpha=0.7)
        if label:
            ax.text(ax.get_xlim()[0], y + 2, label,
                    color=color, fontsize=8)


# ─────────────────────────────────────────────
# SCATTER PLOT HELPER
# ─────────────────────────────────────────────
def draw_scatter(ax, df_plot: pd.DataFrame,
                 anomaly_masks: dict,
                 title: str, y_threshold: float):
    """
    Draw F-stat vs Y-count scatter on ax.
    anomaly_masks: dict of label -> boolean mask Series
    """
    # Plot normal samples
    for sex in [1, 2, 0]:
        grp = df_plot[df_plot['PEDSEX'] == sex]
        if len(grp) == 0:
            continue
        label = {1: 'Male', 2: 'Female', 0: 'Unknown'}[sex]
        ax.scatter(grp['F'], grp['YCOUNT'],
                   c=COLORS[sex], s=12, alpha=0.65,
                   label=label, zorder=2)

    # Overlay anomalies
    for atype, mask in anomaly_masks.items():
        grp = df_plot[mask]
        if len(grp) == 0:
            continue
        ax.scatter(grp['F'], grp['YCOUNT'],
                   c=ANOMALY_COLORS.get(atype, '#ffcc00'),
                   s=150, marker='*', zorder=5,
                   edgecolors='#333333', linewidths=0.5,
                   label=atype.replace('_', ' ').title())

    # Y threshold line
    ax.axhline(y_threshold, color='#666666',
               linestyle=':', linewidth=1.0, alpha=0.7)
    ax.text(ax.get_xlim()[0],
            y_threshold + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.02,
            f'Y mean = {y_threshold:.1f}',
            color='#666666', fontsize=8)

    ax.set_xlabel('F Statistic', fontweight='bold', fontsize=10)
    ax.set_ylabel('Y Count',     fontweight='bold', fontsize=10)
    ax.set_title(title, fontweight='bold', fontsize=12, pad=8)
    apply_style(ax)


# ─────────────────────────────────────────────
# DETECTION
# ─────────────────────────────────────────────
def detect_anomalies(df: pd.DataFrame,
                     f_female_max: float,
                     f_male_min: float,
                     f_turner: float,
                     f_klinefelter: float,
                     f_xxx: float,
                     y_threshold: float) -> dict:
    """
    Returns dict of anomaly label -> DataFrame of flagged samples.
    """
    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]

    results = {
        'sex_mismatch_females' : df_femal[
            (df_femal['F'] > f_male_min) &
            (df_femal['YCOUNT'] > y_threshold)
        ],
        'sex_mismatch_males'   : df_mal[
            (df_mal['F'] < f_female_max) &
            (df_mal['YCOUNT'] < y_threshold)
        ],
        'turner_x0'            : df_femal[
            (df_femal['F'] > f_turner) &
            (df_femal['YCOUNT'] < y_threshold)
        ],
        'klinefelter_xxy'      : df_mal[
            (df_mal['F'] < f_klinefelter) &
            (df_mal['YCOUNT'] > y_threshold)
        ],
        'triple_x'             : df_femal[
            (df_femal['F'] < f_xxx) &
            (df_femal['YCOUNT'] < y_threshold)
        ],
    }
    return results


# ─────────────────────────────────────────────
# WRITE REPORT
# ─────────────────────────────────────────────
def write_report(df: pd.DataFrame,
                 anomalies: dict,
                 y_threshold: float,
                 y_use_mean: str,
                 params: dict,
                 output_dir: str):

    report_path = os.path.join(
        output_dir, "report_sex.mismatch_aneuploidies.txt"
    )

    with open(report_path, 'w') as out:
        out.write("=" * 64 + "\n")
        out.write("  CRISP -- STEP 5 SEX CHECK AND ANEUPLOIDY REPORT\n")
        out.write("  Comprehensive Robust Integrated SNP Processing\n")
        out.write("=" * 64 + "\n\n")
        out.write(f"  F female max   : {params['f_female_max']}\n")
        out.write(f"  F male min     : {params['f_male_min']}\n")
        out.write(f"  Y threshold    : {y_threshold:.2f} "
                  f"({'mean' if y_use_mean == 'YES' else 'manual'})\n\n")

        df_mal   = df[df['PEDSEX'] == 1]
        df_femal = df[df['PEDSEX'] == 2]

        out.write("-" * 64 + "\n")
        out.write("SAMPLE COUNTS\n")
        out.write("-" * 64 + "\n")
        out.write(f"  Reported males   : {len(df_mal)}\n")
        out.write(f"  Reported females : {len(df_femal)}\n\n")

        out.write("-" * 64 + "\n")
        out.write("SEX MISMATCHES\n")
        out.write("-" * 64 + "\n")
        for key in ['sex_mismatch_females', 'sex_mismatch_males']:
            grp = anomalies[key]
            label = key.replace('_', ' ').title()
            out.write(f"  {label} : {len(grp)}\n")
            if len(grp) > 0:
                out.write("\n  Details:\n")
                out.write(grp.to_string(index=False))
                out.write("\n\n")

        out.write("-" * 64 + "\n")
        out.write("CHROMOSOMAL ANEUPLOIDIES\n")
        out.write("-" * 64 + "\n")
        labels = {
            'turner_x0'       : 'Turner syndrome (X0)',
            'klinefelter_xxy' : 'Klinefelter syndrome (XXY)',
            'triple_x'        : 'Triple-X syndrome (XXX)',
        }
        for key, label in labels.items():
            grp = anomalies[key]
            out.write(f"  {label} : {len(grp)}\n")
            if len(grp) > 0:
                out.write("\n  Details:\n")
                out.write(grp.to_string(index=False))
                out.write("\n\n")

        # Combined aneuploidies
        df_aneu = pd.concat([
            anomalies['turner_x0'],
            anomalies['klinefelter_xxy'],
            anomalies['triple_x']
        ]).drop_duplicates()

        out.write("-" * 64 + "\n")
        out.write(f"  Total aneuploidies flagged : {len(df_aneu)}\n")
        if len(df_aneu) > 0:
            out.write("\n  All aneuploidy samples:\n")
            out.write(df_aneu.to_string(index=False))
            out.write("\n")

        out.write("\n" + "=" * 64 + "\n")
        out.write("  END OF REPORT\n")
        out.write("=" * 64 + "\n")

    log(f"Report written: {report_path}")
    return report_path


# ─────────────────────────────────────────────
# WRITE EXCLUSION LISTS
# ─────────────────────────────────────────────
def write_exclusion_lists(anomalies: dict, output_dir: str):

    # Sex mismatch list
    df_mismatch = pd.concat([
        anomalies['sex_mismatch_females'],
        anomalies['sex_mismatch_males']
    ]).drop_duplicates()

    mismatch_file = os.path.join(output_dir, "id_list_sex_mismatch.txt")
    df_mismatch[['FID', 'IID']].to_csv(
        mismatch_file, sep='\t', index=False, header=False
    )
    log(f"Sex mismatch exclusion list: {mismatch_file} "
        f"({len(df_mismatch)} samples)")

    # Aneuploidy list
    df_aneu = pd.concat([
        anomalies['turner_x0'],
        anomalies['klinefelter_xxy'],
        anomalies['triple_x']
    ]).drop_duplicates()

    aneuploidy_file = os.path.join(output_dir, "id_list_aneuploidies.txt")
    df_aneu[['FID', 'IID']].to_csv(
        aneuploidy_file, sep='\t', index=False, header=False
    )
    log(f"Aneuploidy exclusion list  : {aneuploidy_file} "
        f"({len(df_aneu)} samples)")

    return mismatch_file, aneuploidy_file


# ─────────────────────────────────────────────
# GENERATE PLOTS
# ─────────────────────────────────────────────
def generate_plots(df: pd.DataFrame,
                   anomalies: dict,
                   params: dict,
                   y_threshold: float,
                   output_dir: str):

    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]

    pdf_file = os.path.join(output_dir, "Sex_check.pdf")

    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(pdf_file) as pdf:

        # ── OVERALL ──────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 7))

        anomaly_masks_overall = {
            'mismatch_f'  : df.index.isin(anomalies['sex_mismatch_females'].index),
            'mismatch_m'  : df.index.isin(anomalies['sex_mismatch_males'].index),
            'turner'      : df.index.isin(anomalies['turner_x0'].index),
            'klinefelter' : df.index.isin(anomalies['klinefelter_xxy'].index),
            'xxx'         : df.index.isin(anomalies['triple_x'].index),
        }

        draw_scatter(ax, df, anomaly_masks_overall,
                     "Overall Sex Check", y_threshold)
        add_threshold_line(ax, x=params['f_female_max'],
                           color='#1D9E75',
                           label=f"F={params['f_female_max']}")
        add_threshold_line(ax, x=params['f_male_min'],
                           color='#e05c4b',
                           label=f"F={params['f_male_min']}")

        n_aneu = sum(len(v) for k,v in anomalies.items()
                     if k not in ('sex_mismatch_females','sex_mismatch_males'))
        n_mis  = (len(anomalies['sex_mismatch_females']) +
                  len(anomalies['sex_mismatch_males']))

        fig.suptitle(
            f"500 samples  |  {n_mis} sex mismatches  |  "
            f"{n_aneu} aneuploidies flagged",
            fontsize=9, color='#555555', y=0.97
        )
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── MALE ─────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 6))
        anomaly_masks_male = {
            'mismatch_m'  : df_mal.index.isin(anomalies['sex_mismatch_males'].index),
            'klinefelter' : df_mal.index.isin(anomalies['klinefelter_xxy'].index),
        }
        draw_scatter(ax, df_mal, anomaly_masks_male,
                     "Male Sex Check", y_threshold)
        add_threshold_line(ax, x=params['f_klinefelter'],
                           color='darkred',
                           label=f"Klinefelter F={params['f_klinefelter']}")
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 1.0])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── FEMALE ───────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 6))
        anomaly_masks_female = {
            'mismatch_f' : df_femal.index.isin(
                anomalies['sex_mismatch_females'].index),
            'turner'     : df_femal.index.isin(anomalies['turner_x0'].index),
            'xxx'        : df_femal.index.isin(anomalies['triple_x'].index),
        }
        draw_scatter(ax, df_femal, anomaly_masks_female,
                     "Female Sex Check", y_threshold)
        add_threshold_line(ax, x=params['f_turner'],
                           color='purple',
                           label=f"Turner F={params['f_turner']}")
        add_threshold_line(ax, x=params['f_xxx'],
                           color='blue',
                           label=f"XXX F={params['f_xxx']}")
        ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color='#888888')
        plt.tight_layout(rect=[0, 0.03, 1, 1.0])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    log(f"Plots saved: {pdf_file} (3 pages: overall, male, female)")
    return pdf_file


# ─────────────────────────────────────────────
# SEX AMENDMENT
# ─────────────────────────────────────────────
def apply_sex_amendment(df: pd.DataFrame, output_dir: str):
    """
    Updates PEDSEX to match SNPSEX for mismatched samples.
    Writes an amended FAM file and a detailed amendment log.
    Original FAM file is never overwritten.
    """
    if 'SNPSEX' not in df.columns:
        log("WARNING: SNPSEX column not found. Cannot apply sex amendment.")
        return

    mismatches = df[
        (df['PEDSEX'] != df['SNPSEX']) &
        (df['SNPSEX'].isin([1, 2]))
    ][['FID', 'IID', 'PEDSEX', 'SNPSEX']]

    if len(mismatches) == 0:
        log("No sex mismatches with valid SNPSEX found. No amendments made.")
        return

    amendment_log = os.path.join(output_dir, "sex_amendments.txt")

    sex_labels = {1: 'Male', 2: 'Female', 0: 'Unknown'}

    with open(amendment_log, 'w') as out:
        out.write("=" * 64 + "\n")
        out.write("  CRISP -- SEX AMENDMENT LOG\n")
        out.write("  Comprehensive Robust Integrated SNP Processing\n")
        out.write("=" * 64 + "\n\n")
        out.write("  WARNING: These samples had their recorded sex updated\n")
        out.write("  to match the genotype-inferred sex (SNPSEX).\n")
        out.write("  Consult the data provider before proceeding.\n\n")
        out.write(f"  {'FID':<12} {'IID':<12} {'Original':<12} {'Amended':<12}\n")
        out.write("-" * 64 + "\n")

        for _, row in mismatches.iterrows():
            orig = sex_labels.get(int(row['PEDSEX']), 'Unknown')
            new  = sex_labels.get(int(row['SNPSEX']), 'Unknown')
            out.write(
                f"  {row['FID']:<12} {row['IID']:<12} "
                f"{orig:<12} {new:<12}\n"
            )
            log(f"WARNING: Amending {row['FID']} {row['IID']} "
                f"SEX {int(row['PEDSEX'])} -> {int(row['SNPSEX'])}")

        out.write("\n" + "=" * 64 + "\n")

    log(f"Amendment log written : {amendment_log}")
    log(f"Samples amended       : {len(mismatches)}")
    log("")
    log("!" * 60)
    log("Please review the amendment log before continuing.")
    log("!" * 60)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRISP Step 5: Sex Check and Aneuploidy Detection (Python)"
    )
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    parser.add_argument("f_stat_y_count_file")
    parser.add_argument("f_female_max",          type=float)
    parser.add_argument("f_male_min",             type=float)
    parser.add_argument("f_turner",               type=float)
    parser.add_argument("f_klinefelter",          type=float)
    parser.add_argument("f_xxx",                  type=float)
    parser.add_argument("y_use_mean")
    parser.add_argument("y_manual",               type=float)
    parser.add_argument("--exclude_mismatch",     default="YES")
    parser.add_argument("--exclude_aneuploidy",   default="YES")
    parser.add_argument("--amend_sex",            default="NO")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    log("CRISP Step 5: Sex Check and Aneuploidy Detection (Python)")
    log(f"Parameters:")
    log(f"  F female max              : {args.f_female_max}")
    log(f"  F male min                : {args.f_male_min}")
    log(f"  F Turner                  : {args.f_turner}")
    log(f"  F Klinefelter             : {args.f_klinefelter}")
    log(f"  F XXX threshold           : {args.f_xxx}")
    log(f"  Y use mean                : {args.y_use_mean}")
    log(f"  Exclude mismatches        : {args.exclude_mismatch}")
    log(f"  Exclude aneuploidies      : {args.exclude_aneuploidy}")
    log(f"  Amend sex in FAM          : {args.amend_sex}")
    if args.y_use_mean.upper() != "YES":
        log(f"  Y manual cutoff           : {args.y_manual}")
    log("")

    # Mutual exclusion check
    if (args.amend_sex.upper() == "YES" and
            args.exclude_mismatch.upper() == "YES"):
        abort("SEXCHECK_AMEND_SEX = YES and SEXCHECK_EXCLUDE_MISMATCH = YES "
              "cannot both be active. Choose one approach.")

    # Load data
    log(f"Loading: {args.f_stat_y_count_file}")
    df = load_data(args.f_stat_y_count_file)
    log(f"Samples loaded: {len(df)}")

    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]
    log(f"Reported males   : {len(df_mal)}")
    log(f"Reported females : {len(df_femal)}")
    log("")

    # Y threshold
    if args.y_use_mean.upper() == "YES":
        y_threshold = float(df['YCOUNT'].mean())
        log(f"Y count threshold: mean = {y_threshold:.2f}")
    else:
        y_threshold = args.y_manual
        log(f"Y count threshold: manual = {y_threshold:.2f}")
    log("")

    params = {
        'f_female_max'  : args.f_female_max,
        'f_male_min'    : args.f_male_min,
        'f_turner'      : args.f_turner,
        'f_klinefelter' : args.f_klinefelter,
        'f_xxx'         : args.f_xxx,
    }

    # Detect anomalies
    log("Running detection...")
    anomalies = detect_anomalies(
        df, args.f_female_max, args.f_male_min,
        args.f_turner, args.f_klinefelter, args.f_xxx,
        y_threshold
    )

    log(f"  Sex mismatch (females)    : {len(anomalies['sex_mismatch_females'])}")
    log(f"  Sex mismatch (males)      : {len(anomalies['sex_mismatch_males'])}")
    log(f"  Turner syndrome (X0)      : {len(anomalies['turner_x0'])}")
    log(f"  Klinefelter syndrome (XXY): {len(anomalies['klinefelter_xxy'])}")
    log(f"  Triple-X syndrome (XXX)   : {len(anomalies['triple_x'])}")
    log("")

    # Generate plots
    log("Generating plots...")
    generate_plots(df, anomalies, params, y_threshold, args.output_dir)

    # Write detection report
    log("Writing detection report...")
    write_report(df, anomalies, y_threshold,
                 args.y_use_mean, params, args.output_dir)

    # Write exclusion lists
    log("Writing exclusion lists...")
    write_exclusion_lists(anomalies, args.output_dir)

    # Sex amendment
    if args.amend_sex.upper() == "YES":
        log("")
        log("!" * 60)
        log("SEX AMENDMENT REQUESTED.")
        log("Updating PEDSEX in FAM file to match genotype-inferred sex.")
        log("It is strongly recommended to consult the data provider")
        log("for each amended sample before proceeding.")
        log("!" * 60)
        log("")
        apply_sex_amendment(df, args.output_dir)

    log("\nSex check complete.")


if __name__ == "__main__":
    main()
