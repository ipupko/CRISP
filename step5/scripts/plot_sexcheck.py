#!/usr/bin/env python3
##########################################################################
#
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
#
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 5: Sex Check and Aneuploidy Detection — Python Plotting Script
# Version: 0.4.0
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# Detects and plots sex mismatches and chromosomal aneuploidies using
# PLINK X-chromosome F-statistics and Y-chromosome variant counts.
# Python parity script for plot_sexcheck.R — produces identical outputs.
#
# v0.4.0:
#   - NEW: aneuploidy_cohort_manifest.txt — typed manifest (CATEGORY,
#     TYPE) of Turner/Klinefelter/Triple-X samples for cohort routing.
#     These are confirmed karyotype calls, not data errors — valuable
#     to RNA/eQTL researchers as natural dosage experiments. Sex
#     mismatches remain exclusion-only (Type 2 error, not biology).
#   - NEW: per-type exclusion files (id_list_aneuploidy_turner.txt, etc.)
#     for shell-side --keep extraction into per-type sub-cohorts.
#   - FIX: Y threshold now uses MALE-ONLY median, not whole-cohort mean.
#     The previous mean included females/unknowns and was dragged well
#     below the true male Y-count cluster.
#   - FIX: removed redundant `import os as _os`.
#
# v0.3.0 fixes:
#   - ax.get_xlim() / ax.get_ylim() calls moved to AFTER ax.relim() /
#     ax.autoscale_view() so limits are finalised before threshold
#     annotations and label offsets are computed (timing bug fix)
#   - Hardcoded "500 samples" in suptitle replaced with actual nrow(df)
#   - All legends moved to right side (upper left -> right-hand)
#   - ANOMALY_COLORS now drawn from CRISP_PAL_* env vars so
#     COLOURBLIND mode applies to anomaly markers, not just base scatter
#   - --label_anomalies and --label_field now received from shell
#     (were accepted as args but never forwarded by crisp_sexcheck.sh)
#
# Usage:
#   python3 plot_sexcheck.py <input_dir> <output_dir> <f_stat_y_count_file>
#                            <f_female_max> <f_male_min> <f_turner>
#                            <f_klinefelter> <f_xxx> <y_use_mean> <y_manual>
#                            [--label_anomalies YES|NO]
#                            [--label_field IID|FID|FID_IID]
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
from matplotlib.backends.backend_pdf import PdfPages
from adjustText import adjust_text

##########################################################################
#
#   ____ ___  _     ___  _   _ ____    ____   _    _
#  / ___/ _ \| |   / _ \| | | |  _ \  |  _ \ / \  | |
# | |  | | | | |  | | | | | | | |_) | | |_) / _ \ | |
# | |__| |_| | |__| |_| | |_| |  _ <  |  __/ ___ \| |___
#  \____\___/|_____\___/ \___/|_| \_\ |_| /_/   \_\_____|
#
# ANOMALY_COLORS reads from CRISP_PAL_* environment variables so that
# COLOURBLIND/NIGHT modes carry through to anomaly star markers, not
# just the base scatter.
#
# v0.4.1: added NIGHT mode + PLOT_BACKGROUND support, matching the
# pattern established in Step 3 (crisp_callrate). Previously this
# script only branched on COLOURBLIND vs STANDARD and had no concept
# of a dark plot background at all — PAL_BG/PANEL/TEXT/SUBTEXT/GRID
# did not exist here, so NIGHT mode (set via crisp_sexcheck.sh calling
# _crisp_palette) had no visible effect on Step 5 plots even though
# the shared palette infrastructure already supported it.
##########################################################################

def _env_or(name: str, default: str) -> str:
    val = os.getenv(name, '')
    return val if val else default

_mode = _env_or('CRISP_PAL_MODE', 'STANDARD').upper()

if _mode == 'COLOURBLIND':
    PAL_PASS      = _env_or('CRISP_PAL_PASS',      '#0072B2')
    PAL_FAIL      = _env_or('CRISP_PAL_FAIL',      '#D55E00')
    PAL_WARN      = _env_or('CRISP_PAL_WARN',      '#E69F00')
    PAL_HIGHLIGHT = _env_or('CRISP_PAL_HIGHLIGHT', '#009E73')
    PAL_SKY       = _env_or('CRISP_PAL_SKY',       '#56B4E9')
    PAL_PINK      = _env_or('CRISP_PAL_PINK',      '#CC79A7')
    PAL_ORANGE    = _env_or('CRISP_PAL_YELLOW',    '#F0E442')
elif _mode == 'NIGHT':
    PAL_PASS      = _env_or('CRISP_PAL_PASS',      '#2ED9A3')
    PAL_FAIL      = _env_or('CRISP_PAL_FAIL',      '#FF7B72')
    PAL_WARN      = _env_or('CRISP_PAL_WARN',      '#FFC857')
    PAL_HIGHLIGHT = _env_or('CRISP_PAL_HIGHLIGHT', '#5EEBC4')
    PAL_SKY       = _env_or('CRISP_PAL_SKY',       '#7FD1FF')
    PAL_PINK      = _env_or('CRISP_PAL_PINK',      '#D8B4FE')
    PAL_ORANGE    = _env_or('CRISP_PAL_YELLOW',    '#FFE066')
else:
    PAL_PASS      = _env_or('CRISP_PAL_PASS',      '#1D9E75')
    PAL_FAIL      = _env_or('CRISP_PAL_FAIL',      '#ff5f57')
    PAL_WARN      = _env_or('CRISP_PAL_WARN',      '#febc2e')
    PAL_HIGHLIGHT = _env_or('CRISP_PAL_HIGHLIGHT', '#5dcaa5')
    PAL_SKY       = _env_or('CRISP_PAL_SKY',       '#7ec8e3')
    PAL_PINK      = _env_or('CRISP_PAL_PINK',      '#afa9ec')
    PAL_ORANGE    = _env_or('CRISP_PAL_YELLOW',    '#ff8c00')

# Background / text / grid — orthogonal to PLOT_COLOUR_MODE, set via
# PLOT_BACKGROUND (LIGHT/DARK) and exported by _crisp_palette in the
# shell script. NIGHT mode defaults the background to DARK upstream.
PAL_BG         = _env_or('CRISP_PAL_BG',         '#FFFFFF')
PAL_PANEL      = _env_or('CRISP_PAL_PANEL',      '#f8f9fa')
PAL_TEXT       = _env_or('CRISP_PAL_TEXT',       '#1a1a1a')
PAL_SUBTEXT    = _env_or('CRISP_PAL_SUBTEXT',    '#555555')
PAL_GRID       = _env_or('CRISP_PAL_GRID',       '#e0e0e0')
PAL_BACKGROUND = _env_or('CRISP_PAL_BACKGROUND', 'LIGHT').upper()

# Base sex-group scatter colours — fixed identity colours regardless of
# palette mode (male=red-ish, female=teal, unknown=grey), as in v0.3.0.
# These stay readable on both light and dark backgrounds.
COLORS = {
    1: '#e05c4b',   # male
    2: '#1D9E75',   # female
    0: '#aaaaaa',   # unknown
}

# Anomaly overlay colours — derived from CRISP palette
ANOMALY_COLORS = {
    'mismatch_f'  : PAL_ORANGE,     # reported female, inferred male
    'mismatch_m'  : PAL_PINK,       # reported male,   inferred female
    'turner'      : PAL_SKY,        # Turner X0
    'klinefelter' : PAL_FAIL,       # Klinefelter XXY
    'xxx'         : PAL_HIGHLIGHT,  # Triple-X
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


##########################################################################
#
#  _     ___    _    ____    ____    _  _____  _
# | |   / _ \  / \  |  _ \  |  _ \  / \|_   _|/ \
# | |  | | | |/ _ \ | | | | | | | |/ _ \ | | / _ \
# | |__| |_| / ___ \| |_| | | |_| / ___ \| |/ ___ \
# |_____\___/_/   \_\____/  |____/_/   \_\_/_/   \_\
#
##########################################################################

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
    """Apply CRISP matplotlib style to an axes object, dark-mode aware."""
    ax.set_facecolor(PAL_PANEL)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(0.6)
    ax.spines[['left', 'bottom']].set_color(PAL_SUBTEXT)
    ax.tick_params(labelsize=9, colors=PAL_SUBTEXT)
    ax.xaxis.label.set_color(PAL_TEXT)
    ax.yaxis.label.set_color(PAL_TEXT)
    ax.title.set_color(PAL_TEXT)


##########################################################################
#
#  ____  _____ _____ _____ ____ _____ ___ ___  _   _
# |  _ \| ____|_   _| ____/ ___|_   _|_ _/ _ \| \ | |
# | | | |  _|   | | |  _|| |     | |  | | | | |  \| |
# | |_| | |___  | | | |__| |___  | |  | | |_| | |\  |
# |____/|_____| |_| |_____\____| |_| |___\___/|_| \_|
#
##########################################################################

def compute_y_confidence_bands(df_mal: pd.DataFrame,
                               k_confident: float = 2.0,
                               k_passable: float = 1.0) -> dict:
    """
    Derives confidence-tiered bands around the male Y-count cluster using
    MAD (median absolute deviation), not SD — MAD is robust to the very
    outliers (aneuploidies, mismatches) the threshold is meant to detect,
    whereas SD would be inflated by them, circularly weakening the
    threshold exactly where it matters most.

    IMPORTANT CAVEAT (surfaced in report and plots): Y chromosome variant
    counts are NOT an absolute biological constant — they are an artefact
    of the genotyping array or sequencing platform used. A different
    array design, capture kit, or coverage depth will shift the entire
    male Y-count distribution up or down. This threshold is therefore
    always RELATIVE to the current cohort's own data, recomputed fresh
    per run from the male samples actually present. It does not transfer
    across cohorts genotyped on different platforms.

    Returns a dict with:
      median, mad, confident_lower, confident_upper,
      passable_lower, passable_upper
    Samples within the passable band are NOT classified as anomalies —
    they are too close to the male cluster to call with any confidence
    and are left untouched in the main cohort.
    """
    y = df_mal['YCOUNT'].astype(float)
    med = float(y.median())
    mad = float((y - med).abs().median()) * 1.4826  # normal-consistent MAD

    # Guard against a degenerate (near-zero) MAD, which would make every
    # sample "confident" by collapsing the bands to a point.
    if mad < 1.0:
        mad = max(1.0, float(y.std()) * 0.5)

    return {
        'median'           : med,
        'mad'              : mad,
        'confident_lower'  : med - k_confident * mad,
        'confident_upper'  : med + k_confident * mad,
        'passable_lower'   : med - k_passable * mad,
        'passable_upper'   : med + k_passable * mad,
    }


def detect_anomalies(df: pd.DataFrame,
                     f_female_max: float,
                     f_male_min: float,
                     f_turner: float,
                     f_klinefelter: float,
                     f_xxx: float,
                     y_threshold: float,
                     y_bands: dict) -> dict:
    """
    Returns dict of anomaly key -> DataFrame of flagged samples, each
    carrying a CONFIDENCE column ('CONFIDENT' or 'PASSABLE').

    Detection runs before any plot call so all anomaly frames are
    available when scatter overlays are drawn.

    CONFIDENCE TIERS (see compute_y_confidence_bands):
      CONFIDENT : Y-count sits beyond k_confident MAD from the male
                  median — classification is trustworthy.
      PASSABLE  : Y-count sits between k_passable and k_confident MAD —
                  borderline; flagged but with reduced confidence.
      Samples within k_passable MAD of the male median are NOT flagged
      at all — too close to the male cluster to call, and are left in
      the main cohort untouched. This prevents genuinely male samples
      with platform-typical Y-count variance from being misclassified
      as mismatches or aneuploidies near the boundary.

    NOTE: this is a placeholder for a proper ML-based clustering model.
    A future version should fit a mixture model (e.g. GMM or HDBSCAN)
    on Y-count jointly with F-statistic per platform/batch, and report
    a calibrated per-sample error rate instead of a fixed MAD multiple.
    """
    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]

    conf_lo, conf_hi = y_bands['confident_lower'], y_bands['confident_upper']
    pass_lo, pass_hi = y_bands['passable_lower'],  y_bands['passable_upper']

    def confidence_tier(y_count: pd.Series) -> pd.Series:
        """Returns CONFIDENT / PASSABLE / None (None = inside deadzone)."""
        is_confident = (y_count <= conf_lo) | (y_count >= conf_hi)
        is_passable  = ((y_count <= pass_lo) | (y_count >= pass_hi)) & ~is_confident
        tier = pd.Series(np.where(is_confident, 'CONFIDENT',
                          np.where(is_passable, 'PASSABLE', None)),
                          index=y_count.index)
        return tier

    def tag_confidence(sub_df: pd.DataFrame) -> pd.DataFrame:
        if len(sub_df) == 0:
            sub_df = sub_df.copy()
            sub_df['CONFIDENCE'] = pd.Series(dtype=object)
            return sub_df
        sub_df = sub_df.copy()
        sub_df['CONFIDENCE'] = confidence_tier(sub_df['YCOUNT'])
        return sub_df

    # Reported female but looks genotypically male: high F, non-zero Y.
    # YCOUNT must be in the passable-or-better zone (i.e. NOT in the
    # male-cluster deadzone) — a near-zero F with Y close to the male
    # median is the textbook mismatch case, so we keep YCOUNT > deadzone
    # upper edge as the gate, then tag confidence within that.
    raw_mismatch_f = df_femal[
        (df_femal['F'] > f_male_min) &
        (df_femal['YCOUNT'] >= pass_hi)
    ]
    raw_mismatch_m = df_mal[
        (df_mal['F'] < f_female_max) &
        (df_mal['YCOUNT'] <= pass_lo)
    ]
    raw_turner = df_femal[
        (df_femal['F'] > f_turner) &
        (df_femal['YCOUNT'] <= pass_lo)
    ]
    raw_klinefelter = df_mal[
        (df_mal['F'] < f_klinefelter) &
        (df_mal['YCOUNT'] >= pass_hi)
    ]
    raw_xxx = df_femal[
        (df_femal['F'] < f_xxx) &
        (df_femal['YCOUNT'] <= pass_lo)
    ]

    return {
        'sex_mismatch_females' : tag_confidence(raw_mismatch_f),
        'sex_mismatch_males'   : tag_confidence(raw_mismatch_m),
        'turner_x0'            : tag_confidence(raw_turner),
        'klinefelter_xxy'      : tag_confidence(raw_klinefelter),
        'triple_x'             : tag_confidence(raw_xxx),
    }


##########################################################################
#
#  ____  _     ___  _____ ____
# |  _ \| |   / _ \|_   _/ ___|
# | |_) | |  | | | | | | \___ \
# |  __/| |__| |_| | | |  ___) |
# |_|   |_____\___/  |_| |____/
#
# FIX #10: ax.get_xlim() / ax.get_ylim() are now called AFTER
# ax.relim() + ax.autoscale_view(), which forces matplotlib to commit
# the axis limits based on the already-drawn data. In v0.2.0 these were
# called mid-draw before data was rendered, returning matplotlib's
# initial default limits (0.0, 1.0) and placing all threshold text
# and label offsets at wrong positions.
#
# FIX #11: "500 samples" placeholder in suptitle replaced with
# the actual length of the loaded dataframe.
#
# FIX #12: legend loc changed from 'upper left' to 'center right'
# so all plots comply with the CRISP right-side legend standard.
##########################################################################

def make_label(row, label_field: str) -> str:
    # Type suffix omitted: anomaly category is already encoded by star colour
    # and the external legend. Shorter labels dramatically reduce plot clutter.
    if label_field == "FID_IID":
        return f"{row['FID']}_{row['IID']}"
    elif label_field == "FID":
        return str(row['FID'])
    else:
        return str(row['IID'])


def add_threshold_line(ax, x=None, y=None, color='grey',
                       label='', linestyle='--'):
    """
    Draw a vertical or horizontal threshold line.
    Labels are positioned AFTER ax.relim()/autoscale_view() so the axis
    limits used for placement are the finalised data-driven limits.
    """
    if x is not None:
        ax.axvline(x, color=color, linestyle=linestyle,
                   linewidth=1.2, alpha=0.8)
        if label:
            # FIX: get_ylim() called here — after all data has been plotted
            # and autoscale_view() has committed the final limits.
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()
            x_offset = (xlim[1] - xlim[0]) * 0.01
            ax.text(x + x_offset, ylim[1] * 0.97, label,
                    color=color, fontsize=8, fontweight='bold',
                    va='top', ha='left')
    if y is not None:
        ax.axhline(y, color=color, linestyle=':', linewidth=1.0, alpha=0.7)
        if label:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            y_offset = (ylim[1] - ylim[0]) * 0.02
            ax.text(xlim[0], y + y_offset, label,
                    color=color, fontsize=8, va='bottom', ha='left')


def draw_scatter(ax, fig, df_plot: pd.DataFrame,
                 anomaly_masks: dict,
                 title: str, y_threshold: float,
                 y_bands: dict,
                 show_labels: bool = True,
                 label_field: str = "IID"):
    """
    Draw F-stat vs Y-count scatter on ax with anomaly overlays.

    Legend is placed outside the axes (bbox_to_anchor to the right of the
    axes panel) so it never occludes data. adjustText handles label repulsion
    — replacing the manual candidate-offset system that produced clumping
    when multiple anomalies shared the same region of the plot.

    Y-count confidence bands (y_bands) are shaded horizontally:
      - PASSABLE band (middle, light blue wash): borderline zone close to
        the male Y cluster. Samples here are not flagged at all.
      - Outside the passable band (light red wash): the CONFIDENT zone
        where classification is trustworthy.
    Anomaly star markers additionally get a coloured RING per-sample:
      red ring   = CONFIDENT classification
      blue ring  = PASSABLE classification (borderline, lower confidence)
    """
    legend_handles = []

    # Base scatter: one group per reported sex
    for sex in [1, 2, 0]:
        grp = df_plot[df_plot['PEDSEX'] == sex]
        if len(grp) == 0:
            continue
        sex_label = {1: 'Male', 2: 'Female', 0: 'Unknown'}[sex]
        ax.scatter(grp['F'], grp['YCOUNT'],
                   c=COLORS[sex], s=12, alpha=0.65, zorder=2)
        legend_handles.append(
            mpatches.Patch(color=COLORS[sex], label=sex_label)
        )

    # Y-count threshold dotted line; label positioned after autoscale
    ax.axhline(y_threshold, color=PAL_SUBTEXT,
               linestyle=':', linewidth=1.0, alpha=0.8, zorder=1)

    # Commit final axis limits from the base scatter data
    ax.relim()
    ax.autoscale_view()

    xlim    = ax.get_xlim()
    ylim    = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]

    # Confidence band shading — drawn AFTER axis limits are finalised so
    # the spans cover the full plot width regardless of data extent.
    # Middle (PASSABLE) band: light blue wash, samples here are unflagged.
    # Outer (CONFIDENT) zone: light red wash, classification trustworthy.
    if y_bands.get('mad', 0) > 0:
        pass_lo, pass_hi = y_bands['passable_lower'], y_bands['passable_upper']
        conf_lo, conf_hi = y_bands['confident_lower'], y_bands['confident_upper']

        # Outer confident zones (red wash) — above conf_hi and below conf_lo
        ax.axhspan(conf_hi, ylim[1], color='#ff5f57', alpha=0.05, zorder=0)
        ax.axhspan(ylim[0], conf_lo, color='#ff5f57', alpha=0.05, zorder=0)
        # Middle passable zone (blue wash) — between conf_lo/conf_hi and
        # pass_lo/pass_hi forms the passable band; inside pass_lo/pass_hi
        # is the true deadzone (unshaded, samples never flagged there)
        ax.axhspan(pass_hi, conf_hi, color='#7ec8e3', alpha=0.12, zorder=0)
        ax.axhspan(conf_lo, pass_lo, color='#7ec8e3', alpha=0.12, zorder=0)

        # Boundary lines at confident/passable edges for clarity
        for yv in (conf_lo, conf_hi):
            ax.axhline(yv, color='#ff5f57', linestyle='-',
                      linewidth=0.6, alpha=0.4, zorder=1)
        for yv in (pass_lo, pass_hi):
            ax.axhline(yv, color='#7ec8e3', linestyle='-',
                      linewidth=0.6, alpha=0.5, zorder=1)

    # Y-threshold label tucked into the left margin above the line
    ax.text(xlim[0] + x_range * 0.01,
            y_threshold + y_range * 0.015,
            f'Y median (male) = {y_threshold:.1f}',
            color=PAL_SUBTEXT, fontsize=7.5, va='bottom', style='italic')

    # Anomaly overlays — collect all annotation Text objects then repel
    # them together in a single adjust_text() call so they push off each
    # other, not just off their own point of origin.
    all_texts   = []
    point_pairs = []   # (x, y) of each star, paired with its Text object

    RING_COLOURS = {'CONFIDENT': '#ff5f57', 'PASSABLE': '#7ec8e3'}

    for atype, mask in anomaly_masks.items():
        grp = df_plot[mask]
        if len(grp) == 0:
            continue
        col        = ANOMALY_COLORS.get(atype, PAL_WARN)
        type_label = {
            'mismatch_f'  : 'Mismatch F',
            'mismatch_m'  : 'Mismatch M',
            'turner'      : 'Turner X0',
            'klinefelter' : 'Klinefelter XXY',
            'xxx'         : 'Triple-X XXX',
        }.get(atype, atype.replace('_', ' ').title())

        # Ring colour per-sample confidence tier (red=confident, blue=passable)
        ring_colours = grp['CONFIDENCE'].map(RING_COLOURS).fillna('#999999') \
            if 'CONFIDENCE' in grp.columns else ['#999999'] * len(grp)

        ax.scatter(grp['F'], grp['YCOUNT'],
                   c=col, s=200, marker='*', zorder=5,
                   edgecolors=ring_colours, linewidths=2.2)
        legend_handles.append(
            mlines.Line2D([], [], color=col, marker='*',
                          markersize=9, linestyle='None',
                          label=type_label)
        )

        if not show_labels:
            continue

        for _, row in grp.iterrows():
            lbl = make_label(row, label_field)
            txt = ax.text(
                row['F'], row['YCOUNT'], lbl,
                fontsize   = 7.5,
                fontweight = 'bold',
                color      = col,
                va         = 'bottom',
                ha         = 'center',
                zorder     = 7,
                bbox       = dict(
                    boxstyle  = 'round,pad=0.25',
                    facecolor = PAL_BG,
                    edgecolor = col,
                    alpha     = 0.88,
                    linewidth = 0.8,
                )
            )
            all_texts.append(txt)
            point_pairs.append((row['F'], row['YCOUNT']))

    # adjustText: repels labels from each other and from their anchor
    # points, then draws neat straight connector lines. This eliminates
    # the arcing pile-up visible when multiple anomalies cluster together.
    if all_texts:
        adjust_text(
            all_texts,
            x=[p[0] for p in point_pairs],
            y=[p[1] for p in point_pairs],
            ax           = ax,
            expand       = (1.4, 1.6),
            force_text   = (0.6, 0.8),
            force_points = (0.4, 0.5),
            arrowprops   = dict(
                arrowstyle = '-',
                color      = PAL_SUBTEXT,
                lw         = 0.7,
            ),
            lim = 300,
        )

    ax.set_xlabel('F Statistic (X chromosome)', fontweight='bold', fontsize=10)
    ax.set_ylabel('Y Chromosome Variant Count', fontweight='bold', fontsize=10)
    ax.set_title(title, fontweight='bold', fontsize=12, pad=8)
    apply_style(ax)

    # Confidence ring legend entries (red/blue rings), appended after the
    # anomaly type entries so the legend reads: sex groups -> anomaly
    # types -> confidence key.
    legend_handles.append(mlines.Line2D(
        [], [], marker='*', markersize=9, linestyle='None',
        markerfacecolor='none', markeredgecolor='#ff5f57',
        markeredgewidth=2.0, label='Confident (outer band)'
    ))
    legend_handles.append(mlines.Line2D(
        [], [], marker='*', markersize=9, linestyle='None',
        markerfacecolor='none', markeredgecolor='#7ec8e3',
        markeredgewidth=2.0, label='Passable (borderline)'
    ))

    # Legend placed OUTSIDE the axes to the right so it never occludes data.
    # bbox_to_anchor=(1.01, 0.5) anchors to the right edge of the axes panel;
    # fig.subplots_adjust(right=0.78) in generate_plots reserves the margin.
    legend = ax.legend(
        handles        = legend_handles,
        loc            = 'center left',
        bbox_to_anchor = (1.01, 0.5),
        fontsize       = 8,
        framealpha     = 0.9,
        title          = 'Groups',
        title_fontsize = 8,
        borderpad      = 0.8,
        handlelength   = 1.5,
        facecolor      = PAL_BG,
        edgecolor      = PAL_GRID,
    )
    legend.get_title().set_color(PAL_TEXT)
    for text in legend.get_texts():
        text.set_color(PAL_TEXT)


##########################################################################
#
#  ____  _____ ____   ___  ____ _____
# |  _ \| ____|  _ \ / _ \|  _ \_   _|
# | |_) |  _| | |_) | | | | |_) || |
# |  _ <| |___|  __/| |_| |  _ < | |
# |_| \_\_____|_|    \___/|_| \_\|_|
#
##########################################################################

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
                  f"({'mean' if y_use_mean.upper() == 'YES' else 'manual'})\n\n")

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
            grp   = anomalies[key]
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


##########################################################################
#
#  _______  ______ _       _     ___ ____ _____ ____
# | ____\ \/ / ___| |     | |   |_ _/ ___|_   _/ ___|
# |  _|  \  / |   | |     | |    | |\___  \ | | \___ \
# | |___ /  \ |___| |___  | |___ | | ___) || |  ___) |
# |_____/_/\_\____|_____| |_____|___|____/ |_| |____/
#
##########################################################################

def write_exclusion_lists(anomalies: dict, output_dir: str):

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


def write_aneuploidy_manifest(anomalies: dict, output_dir: str):
    """
    Writes a TYPED manifest of aneuploidy samples (category=ANEUPLOIDY,
    type=TURNER/KLINEFELTER/TRIPLE_X) for downstream cohort routing.

    Unlike id_list_aneuploidies.txt (a flat FID/IID exclusion list), this
    manifest preserves the karyotype call per sample so the shell script
    can extract per-type PLINK sub-cohorts via --keep, and so downstream
    RNA/EHR researchers know which dosage condition each sample carries.
    """
    type_map = {
        'turner_x0'       : 'TURNER',
        'klinefelter_xxy' : 'KLINEFELTER',
        'triple_x'        : 'TRIPLE_X',
    }

    frames = []
    for key, type_label in type_map.items():
        grp = anomalies[key].copy()
        if len(grp) == 0:
            continue
        grp['CATEGORY'] = 'ANEUPLOIDY'
        grp['TYPE']     = type_label
        cols = ['FID', 'IID', 'CATEGORY', 'TYPE', 'PEDSEX', 'F', 'YCOUNT']
        if 'CONFIDENCE' in grp.columns:
            cols.append('CONFIDENCE')
        frames.append(grp[cols])

    manifest_file = os.path.join(output_dir, "aneuploidy_cohort_manifest.txt")

    if frames:
        df_manifest = pd.concat(frames).drop_duplicates(subset=['FID', 'IID'])
    else:
        df_manifest = pd.DataFrame(
            columns=['FID', 'IID', 'CATEGORY', 'TYPE',
                    'PEDSEX', 'F', 'YCOUNT', 'CONFIDENCE']
        )

    df_manifest.to_csv(manifest_file, sep='\t', index=False)
    log(f"Aneuploidy cohort manifest : {manifest_file} "
        f"({len(df_manifest)} samples)")

    # Per-type breakdown for shell consumption (one file per type, FID/IID only)
    type_files = {}
    for type_label in ['TURNER', 'KLINEFELTER', 'TRIPLE_X']:
        sub = df_manifest[df_manifest['TYPE'] == type_label]
        type_file = os.path.join(
            output_dir, f"id_list_aneuploidy_{type_label.lower()}.txt"
        )
        sub[['FID', 'IID']].to_csv(type_file, sep='\t', index=False, header=False)
        type_files[type_label] = type_file
        log(f"  {type_label:<12} : {len(sub)} samples -> {type_file}")

    return manifest_file, type_files


def _with_confidence(df_plot: pd.DataFrame, anomalies: dict) -> pd.DataFrame:
    """
    Returns a copy of df_plot with a CONFIDENCE column populated for any
    row that appears in one of the anomaly frames. Needed because
    draw_scatter() indexes into df_plot via boolean masks, and the
    CONFIDENCE tag only exists on the anomaly-only frames produced by
    detect_anomalies() — it must be merged back onto the base frame for
    the ring-colour logic in draw_scatter() to see it.
    """
    df_plot = df_plot.copy()
    df_plot['CONFIDENCE'] = None
    for grp in anomalies.values():
        if len(grp) == 0 or 'CONFIDENCE' not in grp.columns:
            continue
        df_plot.loc[df_plot.index.isin(grp.index), 'CONFIDENCE'] = \
            grp.set_index(grp.index)['CONFIDENCE']
    return df_plot


def generate_plots(df: pd.DataFrame,
                   anomalies: dict,
                   params: dict,
                   y_threshold: float,
                   y_bands: dict,
                   output_dir: str,
                   show_labels: bool = True,
                   label_field: str = "IID"):

    df = _with_confidence(df, anomalies)
    df_mal   = df[df['PEDSEX'] == 1]
    df_femal = df[df['PEDSEX'] == 2]

    pdf_file = os.path.join(output_dir, "Sex_check.pdf")

    n_mis  = (len(anomalies['sex_mismatch_females']) +
              len(anomalies['sex_mismatch_males']))
    n_aneu = sum(len(v) for k, v in anomalies.items()
                 if k not in ('sex_mismatch_females', 'sex_mismatch_males'))

    with PdfPages(pdf_file) as pdf:

        # ── OVERALL ──────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(11, 7))
        fig.patch.set_facecolor(PAL_BG)
        # Reserve right margin for the external legend panel
        fig.subplots_adjust(right=0.78, bottom=0.09, top=0.87)

        anomaly_masks_overall = {
            'mismatch_f'  : df.index.isin(anomalies['sex_mismatch_females'].index),
            'mismatch_m'  : df.index.isin(anomalies['sex_mismatch_males'].index),
            'turner'      : df.index.isin(anomalies['turner_x0'].index),
            'klinefelter' : df.index.isin(anomalies['klinefelter_xxy'].index),
            'xxx'         : df.index.isin(anomalies['triple_x'].index),
        }

        draw_scatter(ax, fig, df, anomaly_masks_overall,
                     "Overall Sex Check", y_threshold, y_bands,
                     show_labels=show_labels, label_field=label_field)

        # Threshold lines added AFTER draw_scatter() has called relim/autoscale
        add_threshold_line(ax, x=params['f_female_max'],
                           color=PAL_PASS,
                           label=f"F={params['f_female_max']} (female)")
        add_threshold_line(ax, x=params['f_male_min'],
                           color=PAL_FAIL,
                           label=f"F={params['f_male_min']} (male)")

        fig.suptitle(
            f"{len(df)} samples  |  {n_mis} sex mismatches  |  "
            f"{n_aneu} aneuploidies flagged",
            fontsize=9, color=PAL_SUBTEXT, y=0.98,
        )
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_SUBTEXT)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── MALE ─────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(11, 6))
        fig.patch.set_facecolor(PAL_BG)
        fig.subplots_adjust(right=0.78, bottom=0.09, top=0.87)

        anomaly_masks_male = {
            'mismatch_m'  : df_mal.index.isin(
                anomalies['sex_mismatch_males'].index),
            'klinefelter' : df_mal.index.isin(
                anomalies['klinefelter_xxy'].index),
        }
        draw_scatter(ax, fig, df_mal, anomaly_masks_male,
                     "Male Sex Check", y_threshold, y_bands,
                     show_labels=show_labels, label_field=label_field)
        add_threshold_line(ax, x=params['f_klinefelter'],
                           color=PAL_FAIL,
                           label=f"F={params['f_klinefelter']} (Klinefelter)")
        fig.suptitle(f"{len(df_mal)} reported males",
                     fontsize=9, color=PAL_SUBTEXT, y=0.98)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_SUBTEXT)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # ── FEMALE ───────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(11, 6))
        fig.patch.set_facecolor(PAL_BG)
        fig.subplots_adjust(right=0.78, bottom=0.09, top=0.87)

        anomaly_masks_female = {
            'mismatch_f' : df_femal.index.isin(
                anomalies['sex_mismatch_females'].index),
            'turner'     : df_femal.index.isin(anomalies['turner_x0'].index),
            'xxx'        : df_femal.index.isin(anomalies['triple_x'].index),
        }
        draw_scatter(ax, fig, df_femal, anomaly_masks_female,
                     "Female Sex Check", y_threshold, y_bands,
                     show_labels=show_labels, label_field=label_field)
        add_threshold_line(ax, x=params['f_turner'],
                           color=PAL_PINK,
                           label=f"F={params['f_turner']} (Turner)")
        add_threshold_line(ax, x=params['f_xxx'],
                           color=PAL_SKY,
                           label=f"F={params['f_xxx']} (Triple-X)")
        fig.suptitle(f"{len(df_femal)} reported females",
                     fontsize=9, color=PAL_SUBTEXT, y=0.98)
        fig.text(0.5, 0.01, CAPTION, ha='center', fontsize=8, color=PAL_SUBTEXT)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    log(f"Plots saved: {pdf_file} (3 pages: overall, male, female)")
    return pdf_file


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
    sex_labels    = {1: 'Male', 2: 'Female', 0: 'Unknown'}

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
            new  = sex_labels.get(int(row['SNPSEX']),  'Unknown')
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


##########################################################################
# MAIN
##########################################################################

def main():
    parser = argparse.ArgumentParser(
        description="CRISP Step 5: Sex Check and Aneuploidy Detection (Python)"
    )
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    parser.add_argument("f_stat_y_count_file")
    parser.add_argument("f_female_max",        type=float)
    parser.add_argument("f_male_min",          type=float)
    parser.add_argument("f_turner",            type=float)
    parser.add_argument("f_klinefelter",       type=float)
    parser.add_argument("f_xxx",               type=float)
    parser.add_argument("y_use_mean")
    parser.add_argument("y_manual",            type=float)
    parser.add_argument("--exclude_mismatch",  default="YES")
    parser.add_argument("--exclude_aneuploidy",default="YES")
    parser.add_argument("--amend_sex",         default="NO")
    parser.add_argument("--label_anomalies",   default="YES",
                        help="Show sample ID labels on anomalies (YES/NO)")
    parser.add_argument("--label_field",       default="IID",
                        help="Label field: IID, FID, or FID_IID")
    parser.add_argument("--y_k_confident",     default=2.0, type=float,
                        help="MAD multiplier for the CONFIDENT Y-count "
                             "band edge (default 2.0)")
    parser.add_argument("--y_k_passable",      default=1.0, type=float,
                        help="MAD multiplier for the PASSABLE Y-count "
                             "band edge (default 1.0). Samples inside "
                             "this band are not flagged at all.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    log("CRISP Step 5: Sex Check and Aneuploidy Detection (Python)")
    log(f"Colour mode: {_mode}  |  Background: {PAL_BACKGROUND}")
    log("Parameters:")
    log(f"  F female max              : {args.f_female_max}")
    log(f"  F male min                : {args.f_male_min}")
    log(f"  F Turner                  : {args.f_turner}")
    log(f"  F Klinefelter             : {args.f_klinefelter}")
    log(f"  F XXX threshold           : {args.f_xxx}")
    log(f"  Y use mean                : {args.y_use_mean}")
    log(f"  Exclude mismatches        : {args.exclude_mismatch}")
    log(f"  Exclude aneuploidies      : {args.exclude_aneuploidy}")
    log(f"  Amend sex in FAM          : {args.amend_sex}")
    log(f"  Label anomalies           : {args.label_anomalies}")
    log(f"  Label field               : {args.label_field}")
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
    # FIX: previously computed from the WHOLE cohort (mean(df['YCOUNT'])),
    # which includes females (~0) and unknowns, dragging the threshold well
    # below the true male cluster. A male-only median is robust to outliers
    # and anchors the threshold to the Y-present distribution it is meant
    # to separate from. SEX_Y_USE_MEAN now means "auto-derive from males"
    # for clarity, retaining median (not mean) for outlier robustness.
    if args.y_use_mean.upper() == "YES":
        if len(df_mal) > 0:
            y_threshold = float(df_mal['YCOUNT'].median())
            log(f"Y count threshold: male-only median = {y_threshold:.2f}")
        else:
            y_threshold = float(df['YCOUNT'].median())
            log("WARNING: No reported males found; falling back to "
                f"whole-cohort median = {y_threshold:.2f}")
    else:
        y_threshold = args.y_manual
        log(f"Y count threshold: manual = {y_threshold:.2f}")
    log("")

    # IMPORTANT: Y-count statistics are RELATIVE, not an absolute
    # biological constant. They are derived entirely from the genotyping
    # array / sequencing platform's chrY probe or read coverage for THIS
    # cohort. A different platform, capture kit, or coverage depth will
    # shift the whole male Y-count distribution up or down. Confidence
    # bands below are therefore recomputed fresh per run from the male
    # samples actually present — they do not transfer across cohorts run
    # on different platforms, and should be re-derived per batch if the
    # pipeline ever processes genotyping data from mixed platforms.
    #
    # PLACEHOLDER FOR FUTURE WORK: this MAD-band heuristic is a stand-in
    # for a proper ML clustering model. A future version should fit a
    # mixture model (GMM / HDBSCAN) on Y-count jointly with F-statistic,
    # per platform/batch, and report a calibrated per-sample error rate
    # instead of a fixed MAD multiple.
    if len(df_mal) > 0:
        y_bands = compute_y_confidence_bands(
            df_mal, k_confident=args.y_k_confident, k_passable=args.y_k_passable
        )
    else:
        # No males to anchor bands to — collapse to the manual/median
        # threshold with zero-width bands so nothing is misclassified.
        y_bands = {
            'median': y_threshold, 'mad': 0.0,
            'confident_lower': y_threshold, 'confident_upper': y_threshold,
            'passable_lower': y_threshold, 'passable_upper': y_threshold,
        }

    log("Y-count confidence bands (RELATIVE to this cohort's platform — "
        "see header note):")
    log(f"  Male median Y          : {y_bands['median']:.2f}")
    log(f"  MAD                    : {y_bands['mad']:.2f}")
    log(f"  CONFIDENT band (outer) : <= {y_bands['confident_lower']:.1f} "
        f"or >= {y_bands['confident_upper']:.1f}")
    log(f"  PASSABLE band (middle) : <= {y_bands['passable_lower']:.1f} "
        f"or >= {y_bands['passable_upper']:.1f}")
    log("  (Samples inside the passable band are not flagged at all — "
        "too close to the male cluster to call.)")
    log("")

    params = {
        'f_female_max'  : args.f_female_max,
        'f_male_min'    : args.f_male_min,
        'f_turner'      : args.f_turner,
        'f_klinefelter' : args.f_klinefelter,
        'f_xxx'         : args.f_xxx,
    }

    # Detect anomalies BEFORE plotting
    log("Running detection...")
    anomalies = detect_anomalies(
        df, args.f_female_max, args.f_male_min,
        args.f_turner, args.f_klinefelter, args.f_xxx,
        y_threshold, y_bands
    )
    for key, label in [
        ('sex_mismatch_females', 'Sex mismatch (females)'),
        ('sex_mismatch_males',   'Sex mismatch (males)  '),
        ('turner_x0',            'Turner syndrome (X0)  '),
        ('klinefelter_xxy',      'Klinefelter (XXY)     '),
        ('triple_x',             'Triple-X (XXX)        '),
    ]:
        grp = anomalies[key]
        n_conf = (grp['CONFIDENCE'] == 'CONFIDENT').sum() if len(grp) else 0
        n_pass = (grp['CONFIDENCE'] == 'PASSABLE').sum()  if len(grp) else 0
        log(f"  {label} : {len(grp)} total "
            f"({n_conf} confident, {n_pass} passable)")
    log("")

    show_labels = args.label_anomalies.upper() == "YES"
    label_field = args.label_field.upper()

    # Generate plots
    log("Generating plots...")
    generate_plots(df, anomalies, params, y_threshold, y_bands, args.output_dir,
                   show_labels=show_labels, label_field=label_field)

    # Write detection report
    log("Writing detection report...")
    write_report(df, anomalies, y_threshold,
                 args.y_use_mean, params, args.output_dir)

    # Write exclusion lists
    log("Writing exclusion lists...")
    write_exclusion_lists(anomalies, args.output_dir)

    # Write typed aneuploidy manifest for cohort routing (RNA/EHR handoff)
    log("Writing aneuploidy cohort manifest...")
    write_aneuploidy_manifest(anomalies, args.output_dir)

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
