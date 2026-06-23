#!/bin/bash
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
# CHUNK: Sex Check and Aneuploidy Detection (Step 5)
# Version: 0.4.0
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# Performs sex verification and chromosomal aneuploidy detection
# using PLINK F-statistics and Y chromosome variant counts.
#
# Steps:
#   5a: PLINK --check-sex (generates F-statistics per sample)
#   5b: PLINK Y chromosome count (variant missingness on chrY)
#   5c: Merge F-stat and Y-count data into single input for R/Python
#   5d: Run plot_sexcheck.R or plot_sexcheck.py for detection and plotting
#   5e: Count flagged samples
#   5f: Apply exclusions
#   5g: Route aneuploidy samples to dedicated sub-cohorts (optional)
#   5h: Sex amendment in FAM file (optional)
#   5i: Write Step 5 report
#
# Outputs:
#   Sex_check.pdf                           three scatter plots
#   report_sex.mismatch_aneuploidies.txt    detailed text report
#   id_list_sex_mismatch.txt                sex mismatch exclusions
#   id_list_aneuploidies.txt                aneuploidy exclusions
#   aneuploidy_cohort_manifest.txt          typed aneuploidy manifest
#   aneuploidy_cohort/                      routed sub-cohorts (if enabled)
#
# Thresholds (all overridable in crisp_instructions.txt):
#   SEX_F_FEMALE_MAX   : F threshold for female classification
#   SEX_F_MALE_MIN     : F threshold for male classification
#   SEX_F_TURNER       : F threshold for Turner (X0) detection
#   SEX_F_KLINEFELTER  : F threshold for Klinefelter (XXY) detection
#   SEX_F_XXX          : F threshold for triple-X detection
#   SEX_Y_USE_MEAN     : use male-only median Y count as threshold
#   SEX_Y_MANUAL       : manual Y count threshold if not using the median
#
# v0.4.0:
#   - NEW: SEXCHECK_ROUTE_ANEUPLOIDIES — when YES, Turner/Klinefelter/
#     Triple-X samples are extracted into per-type PLINK sub-cohorts
#     (both Step 4 QC'd and pre-QC raw genotypes, compressed) instead of
#     being discarded. These are confirmed karyotype calls, valuable to
#     RNA/eQTL researchers as natural dosage experiments. Sex mismatches
#     are NEVER routed — Type 2 error, exclusion-only, unchanged.
#   - FIX: set -e was making the Y-count PLINK failure fallback
#     unreachable; now bracketed with set +e/-e so it actually fires.
#   - FIX: N_AMENDED now correctly reflects the amendment count (was
#     always 0 in the report; the Python heredoc count never made it
#     back to the parent shell).
#
# Usage:
#   bash crisp_sexcheck.sh
#   bash crisp_sexcheck.sh --config my_project.txt
##########################################################################

set -euo pipefail

##########################################################################
#
#  _     ___    _    ____    _____ _        ___     _____  _   _ ____
# | |   / _ \  / \  |  _ \  |  ___| |      / \ \   / / _ \| | | |  _ \
# | |  | | | |/ _ \ | | | | | |_  | |     / _ \ \ / / | | | | | | |_) |
# | |__| |_| / ___ \| |_| | |  _| | |___ / ___ \ V /| |_| | |_| |  _ <
# |_____\___/_/   \_\____/  |_|   |_____/_/   \_\_/  \___/ \___/|_| \_\
#
# Source _crisp_flavour.sh at the top of every CRISP script.
# This sets PREP_MSG / LAUNCH_MSG and exports CRISP_PAL_* colour
# variables for the current PLOT_COLOUR_MODE. Non-negotiable.
##########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime

##########################################################################
#
#  ____   _    ____  ____  _____   ____   _    ____      _    __  __ ____
# |  _ \ / \  |  _ \/ ___|| ____| |  _ \ / \  |  _ \    / \  |  \/  / ___|
# | |_) / _ \ | |_) \___ \|  _|   | |_) / _ \ | |_) |  / _ \ | |\/| \___ \
# |  __/ ___ \|  _ < ___) | |___  |  __/ ___ \|  _ <  / ___ \| |  | |___) |
# |_| /_/   \_\_| \_\____/|_____| |_| /_/   \_\_| \_\/_/   \_\_|  |_|____/
#
##########################################################################

INSTRUCTION_FILE="crisp_instructions.txt"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            INSTRUCTION_FILE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: bash crisp_sexcheck.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_sexcheck.sh [--config <file>]"
            exit 1
            ;;
    esac
done

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_sexcheck.sh --config <path>"
    exit 1
fi

parse_param() {
    local key="$1"
    local default="${2:-}"
    local value
    value=$(grep -E "^[[:space:]]*${key}[[:space:]]*=" "${INSTRUCTION_FILE}" \
        | head -1 \
        | sed 's/^[^=]*=//; s/#.*//' \
        | tr -d '[:space:]')
    if [ -z "${value}" ]; then
        echo "${default}"
    else
        echo "${value}"
    fi
}

OUTPUT_DIR=$(parse_param "OUTPUT_DIR" "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME" "project")
KEEP_INTERMEDIATE=$(parse_param "KEEP_INTERMEDIATE" "YES")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# Parse PLOT_COLOUR_MODE and PLOT_BACKGROUND, forward both to _crisp_palette.
# _crisp_palette is shared infrastructure (sourced via _crisp_flavour.sh) and
# already supports STANDARD / COLOURBLIND / NIGHT modes plus an orthogonal
# LIGHT / DARK background, established in Step 3. Step 5 previously only
# ever called _crisp_palette with the mode argument, silently dropping
# PLOT_BACKGROUND and never reaching NIGHT mode's dark-background default.
PLOT_COLOUR_MODE=$(parse_param "PLOT_COLOUR_MODE" "STANDARD")
PLOT_BACKGROUND=$(parse_param "PLOT_BACKGROUND" "")
_crisp_palette "${PLOT_COLOUR_MODE}" "${PLOT_BACKGROUND}"

# FIX #2: parse label parameters so they can be forwarded to plot scripts.
# These were parsed from crisp_instructions.txt but never read, so the
# instruction file settings had zero effect on plot output.
SEXCHECK_LABEL_ANOMALIES=$(parse_param "SEXCHECK_LABEL_ANOMALIES" "YES")
SEXCHECK_LABEL_FIELD=$(parse_param "SEXCHECK_LABEL_FIELD" "IID")

# sex check thresholds
SEX_F_FEMALE_MAX=$(parse_param "SEX_F_FEMALE_MAX"   "0.2")
SEX_F_MALE_MIN=$(parse_param "SEX_F_MALE_MIN"       "0.8")
SEX_F_TURNER=$(parse_param "SEX_F_TURNER"           "0.8")
SEX_F_KLINEFELTER=$(parse_param "SEX_F_KLINEFELTER" "0.8")
SEX_F_XXX=$(parse_param "SEX_F_XXX"                 "-0.15")
SEX_Y_USE_MEAN=$(parse_param "SEX_Y_USE_MEAN"       "YES")
SEX_Y_MANUAL=$(parse_param "SEX_Y_MANUAL"           "0")

SEXCHECK_EXCLUDE_MISMATCH=$(parse_param "SEXCHECK_EXCLUDE_MISMATCH"   "YES")
SEXCHECK_EXCLUDE_ANEUPLOIDY=$(parse_param "SEXCHECK_EXCLUDE_ANEUPLOIDY" "YES")
SEXCHECK_AMEND_SEX=$(parse_param "SEXCHECK_AMEND_SEX" "NO")

# Aneuploidy cohort routing (v0.4.0): Turner/Klinefelter/Triple-X samples
# are confirmed karyotype calls, not data errors — valuable to RNA/eQTL
# researchers as natural dosage experiments. When enabled, instead of
# being silently discarded they are extracted into a dedicated per-type
# sub-cohort (both Step 4 QC'd and pre-QC raw genotypes, compressed)
# alongside a typed manifest for EHR linkage handoff. Sex mismatches are
# NEVER routed — they are a Type 2 error category and stay exclusion-only.
SEXCHECK_ROUTE_ANEUPLOIDIES=$(parse_param "SEXCHECK_ROUTE_ANEUPLOIDIES" "NO")
SEXCHECK_ANEUPLOIDY_COHORT_DIR=$(parse_param "SEXCHECK_ANEUPLOIDY_COHORT_DIR" \
    "${OUTPUT_DIR}/step5_sexcheck/aneuploidy_cohort")

# resolve input: Step 4 output if available, else Step 2
STEP4_FINAL=$(parse_param "STEP4_FINAL_PREFIX" "")
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")

if [ -n "${STEP4_FINAL}" ] && [ -f "${STEP4_FINAL}.bed" ]; then
    INPUT_PREFIX="${STEP4_FINAL}"
elif [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    INPUT_PREFIX="${CONVERTED_PREFIX}"
else
    INPUT_PREFIX=$(ls -t "${OUTPUT_DIR}/step4_snprate/${PROJECT_NAME}"*maf*.bed \
        2>/dev/null | head -1 | sed 's/\.bed$//')
    if [ -z "${INPUT_PREFIX}" ]; then
        INPUT_PREFIX=$(ls -t "${OUTPUT_DIR}/step4_snprate/${PROJECT_NAME}"*.bed \
            2>/dev/null | head -1 | sed 's/\.bed$//')
    fi
    if [ -z "${INPUT_PREFIX}" ]; then
        INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
    fi
fi

# Pre-QC raw genotype prefix — used only for aneuploidy cohort routing,
# so RNA/EHR researchers can choose to re-derive QC on this special-case
# cohort rather than inherit standard call-rate/MAF/HWE filtering, which
# was never designed with known aneuploidies in mind.
if [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    PREQC_PREFIX="${CONVERTED_PREFIX}"
else
    PREQC_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
fi

##########################################################################
# HEADER
##########################################################################

echo ""
echo "[CRISP] ${PREP_MSG}"
echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""
echo "[CRISP] Parameters:"
echo "        SEXCHECK_EXCLUDE_MISMATCH   : ${SEXCHECK_EXCLUDE_MISMATCH}"
echo "        SEXCHECK_EXCLUDE_ANEUPLOIDY : ${SEXCHECK_EXCLUDE_ANEUPLOIDY}"
echo "        SEXCHECK_ROUTE_ANEUPLOIDIES : ${SEXCHECK_ROUTE_ANEUPLOIDIES}"
echo "        SEXCHECK_AMEND_SEX          : ${SEXCHECK_AMEND_SEX}"
echo "        SEXCHECK_LABEL_ANOMALIES    : ${SEXCHECK_LABEL_ANOMALIES}"
echo "        SEXCHECK_LABEL_FIELD        : ${SEXCHECK_LABEL_FIELD}"
echo "        SEX_F_FEMALE_MAX            : ${SEX_F_FEMALE_MAX}"
echo "        SEX_F_MALE_MIN              : ${SEX_F_MALE_MIN}"
echo "        SEX_F_TURNER                : ${SEX_F_TURNER}"
echo "        SEX_F_KLINEFELTER           : ${SEX_F_KLINEFELTER}"
echo "        SEX_F_XXX                   : ${SEX_F_XXX}"
echo "        SEX_Y_USE_MEAN              : ${SEX_Y_USE_MEAN}"
echo "        SEX_Y_MANUAL                : ${SEX_Y_MANUAL}"
echo "        INPUT_PREFIX                : ${INPUT_PREFIX}"
echo "        PLOT_ENGINE                 : ${PLOT_ENGINE}"
echo "        PLOT_COLOUR_MODE            : ${PLOT_COLOUR_MODE}"
echo "        PLOT_BACKGROUND              : ${CRISP_PAL_BACKGROUND}"
echo ""

##########################################################################
#
# __     ___    _     ___ ____    _  _____ _____
# \ \   / / \  | |   |_ _|  _ \  / \|_   _| ____|
#  \ \ / / _ \ | |    | || | | |/ _ \ | | |  _|
#   \ V / ___ \| |___ | || |_| / ___ \| | | |___
#    \_/_/   \_\_____|___|____/_/   \_\_| |_____|
#
# Check for logical conflicts before touching any data.
# SEXCHECK_AMEND_SEX and SEXCHECK_EXCLUDE_MISMATCH are mutually exclusive:
# you either correct the FAM file or you exclude the sample — not both.
##########################################################################

if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ] && \
   [ "${SEXCHECK_EXCLUDE_MISMATCH^^}" = "YES" ]; then
    echo "[CRISP] ERROR: SEXCHECK_AMEND_SEX = YES and SEXCHECK_EXCLUDE_MISMATCH = YES"
    echo "[CRISP]        cannot both be active at the same time."
    echo "[CRISP]        Choose one approach in ${INSTRUCTION_FILE}:"
    echo "[CRISP]          To correct sex in FAM file : SEXCHECK_EXCLUDE_MISMATCH = NO"
    echo "[CRISP]          To exclude mismatches      : SEXCHECK_AMEND_SEX = NO"
    exit 1
fi

for ext in .bed .bim .fam; do
    if [ ! -f "${INPUT_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${INPUT_PREFIX}${ext}"
        echo "[CRISP]        Ensure Step 4 (crisp_snprate.sh) has been run first."
        exit 1
    fi
done

echo "[OK]    Input files verified."
echo ""

SAMPLES=$(wc -l < "${INPUT_PREFIX}.fam" | tr -d '[:space:]')
VARIANTS=$(wc -l < "${INPUT_PREFIX}.bim" | tr -d '[:space:]')

echo "[CRISP] Samples  : ${SAMPLES}"
echo "[CRISP] Variants : ${VARIANTS}"
echo ""

##########################################################################
# SET UP DIRECTORIES
##########################################################################

STEP5_DIR="${OUTPUT_DIR}/step5_sexcheck"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step5_sexcheck_report.txt"

mkdir -p "${STEP5_DIR}"
mkdir -p "${LOG_DIR}"

##########################################################################
#
#   ____ _   _ _____ ____ _  ______  _______  __
#  / ___| | | | ____/ ___| |/ / ___|| ____\ \/ /
# | |   | |_| |  _|| |   | ' /\___ \|  _|  \  /
# | |___|  _  | |__| |___| . \ ___) | |___ /  \
#  \____|_| |_|_____\____|_|\_\____/|_____/_/\_\
#
# PLINK --check-sex computes an X-chromosome inbreeding coefficient
# (F-statistic) per sample. Values near 1 indicate male (single X),
# values near 0 indicate female (two X chromosomes).
##########################################################################

SEXCHECK_PREFIX="${STEP5_DIR}/${PROJECT_NAME}_sexcheck"
LOG_5A="${LOG_DIR}/step5a_checksex.log"

echo "[CRISP] Running PLINK --check-sex..."

${PLINK1} \
    --bfile "${INPUT_PREFIX}" \
    --check-sex \
    --out "${SEXCHECK_PREFIX}" \
    >> "${LOG_5A}" 2>&1

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: PLINK --check-sex failed. Check log: ${LOG_5A}"
    exit 1
fi

echo "[OK]    Sex check file generated: ${SEXCHECK_PREFIX}.sexcheck"
echo ""

##########################################################################
#
# __   __   ____ _   _ ____     ____ ___  _   _ _   _ _____
# \ \ / /  / ___| | | |  _ \   / ___/ _ \| | | | \ | |_   _|
#  \ V /  | |   | |_| | |_) | | |  | | | | | | |  \| | | |
#   | |   | |___|  _  |  _ <  | |__| |_| | |_| | |\  | | |
#   |_|    \____|_| |_|_| \_\  \____\___/ \___/|_| \_| |_|
#
# Count the number of non-missing variants on the Y chromosome per
# sample. Males carry a Y and should show non-zero counts; females
# should show near-zero. Used jointly with the F-statistic to detect
# sex chromosome aneuploidies (Turner X0, Klinefelter XXY, etc.).
##########################################################################

YCHR_PREFIX="${STEP5_DIR}/${PROJECT_NAME}_ychr"
LOG_5B="${LOG_DIR}/step5b_ychr.log"

echo "[CRISP] Counting Y chromosome variants per sample..."

# FIX v0.4.0: set -e was making the fallback below unreachable, since the
# script would already have exited on PLINK failure before this check ran.
# Temporarily disable errexit around this single call so the graceful
# "chrY may not be present" fallback actually executes.
set +e
${PLINK1} \
    --bfile "${INPUT_PREFIX}" \
    --chr Y \
    --missing \
    --out "${YCHR_PREFIX}" \
    >> "${LOG_5B}" 2>&1
YCHR_STATUS=$?
set -e

if [ ${YCHR_STATUS} -ne 0 ]; then
    echo "[CRISP] WARNING: Y chromosome count failed, chrY may not be present."
    echo "[CRISP]          Creating empty Y count file."
    touch "${YCHR_PREFIX}.imiss"
fi

echo "[OK]    Y chromosome count complete."
echo ""

##########################################################################
#
#  __  __ _____ ____   ____ _____
# |  \/  | ____|  _ \ / ___| ____|
# | |\/| |  _| | |_) | |  _|  _|
# | |  | | |___|  _ <| |_| | |___
# |_|  |_|_____|_| \_\\____|_____|
#
# Join the .sexcheck F-statistic table with the Y-chromosome variant
# counts into a single TSV. This combined file is the sole input to
# both plot scripts, keeping the detection logic self-contained.
##########################################################################

MERGED_FILE="${STEP5_DIR}/${PROJECT_NAME}_f_stat_y_count.txt"

echo "[CRISP] Merging F-stat and Y-count data..."

python3 - << PYEOF
import pandas as pd
import os
import sys

sexcheck_file = "${SEXCHECK_PREFIX}.sexcheck"
ychr_file     = "${YCHR_PREFIX}.imiss"
output_file   = "${MERGED_FILE}"

if not os.path.isfile(sexcheck_file):
    print(f"[CRISP] ERROR: .sexcheck file not found: {sexcheck_file}")
    sys.exit(1)

df_sex = pd.read_csv(sexcheck_file, sep=r'\s+')

if os.path.isfile(ychr_file) and os.path.getsize(ychr_file) > 0:
    df_y = pd.read_csv(ychr_file, sep=r'\s+')
    df_y['YCOUNT'] = df_y['N_GENO'] - df_y['N_MISS']
    df_y = df_y[['IID', 'YCOUNT']]
    df_merged = df_sex.merge(df_y, on='IID', how='left')
    df_merged['YCOUNT'] = df_merged['YCOUNT'].fillna(0).astype(int)
else:
    print("[CRISP] WARNING: No Y count data found. Setting YCOUNT to 0.")
    df_sex['YCOUNT'] = 0
    df_merged = df_sex

df_merged.to_csv(output_file, sep='\t', index=False)
print(f"[CRISP] Merged file written: {output_file}")
print(f"[CRISP] Samples in merged file: {len(df_merged)}")
PYEOF

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Merge step failed."
    exit 1
fi

echo ""

##########################################################################
#
#  ____  _____ _____ _____ ____ _____   ____  _     ___ _____
# |  _ \| ____|_   _| ____/ ___|_   _| |  _ \| |   / _ \_   _|
# | | | |  _|   | | |  _|| |     | |   | |_) | |  | | | || |
# | |_| | |___  | | | |__| |___  | |   |  __/| |__| |_| || |
# |____/|_____| |_| |_____\____| |_|   |_|   |_____\___/ |_|
#
# Dispatch to R or Python depending on PLOT_ENGINE. Both scripts
# receive identical arguments and produce identical outputs.
#
# FIX #3: SEXCHECK_LABEL_ANOMALIES and SEXCHECK_LABEL_FIELD are now
# forwarded as positional args (R: args[11], args[12]) and named
# flags (Python: --label_anomalies, --label_field). Previously these
# instruction file settings were silently ignored.
##########################################################################

LOG_5D="${LOG_DIR}/step5d_plots.log"

if [ "${PLOT_ENGINE^^}" = "R" ]; then
    echo "[CRISP] Running plot_sexcheck.R..."
    ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_sexcheck.R" \
        "${STEP5_DIR}" \
        "${STEP5_DIR}" \
        "${MERGED_FILE}" \
        "${SEX_F_FEMALE_MAX}" \
        "${SEX_F_MALE_MIN}" \
        "${SEX_F_TURNER}" \
        "${SEX_F_KLINEFELTER}" \
        "${SEX_F_XXX}" \
        "${SEX_Y_USE_MEAN}" \
        "${SEX_Y_MANUAL}" \
        "${SEXCHECK_LABEL_FIELD}" \
        "${SEXCHECK_LABEL_ANOMALIES}" \
        >> "${LOG_5D}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: plot_sexcheck.R failed. Check log: ${LOG_5D}"
        exit 1
    fi

elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
    echo "[CRISP] Running plot_sexcheck.py..."
    python3 "${SCRIPT_DIR}/scripts/plot_sexcheck.py" \
        "${STEP5_DIR}" \
        "${STEP5_DIR}" \
        "${MERGED_FILE}" \
        "${SEX_F_FEMALE_MAX}" \
        "${SEX_F_MALE_MIN}" \
        "${SEX_F_TURNER}" \
        "${SEX_F_KLINEFELTER}" \
        "${SEX_F_XXX}" \
        "${SEX_Y_USE_MEAN}" \
        "${SEX_Y_MANUAL}" \
        --label_anomalies "${SEXCHECK_LABEL_ANOMALIES}" \
        --label_field     "${SEXCHECK_LABEL_FIELD}" \
        >> "${LOG_5D}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: plot_sexcheck.py failed. Check log: ${LOG_5D}"
        exit 1
    fi
fi

echo "[OK]    Detection and plotting complete."
echo ""

##########################################################################
#
#  _______  ______ _    _   _ ____ ___ ___  _   _ ____
# | ____\ \/ / ___| |  | | | / ___|_ _/ _ \| \ | / ___|
# |  _|  \  / |   | |  | | | \___ \| | | | |  \| \___ \
# | |___ /  \ |___| |__| |_| |___) | | |_| | |\  |___) |
# |_____/_/\_\____|_____\___/|____/___\___/|_| \_|____/
#
# Read the exclusion lists written by the plot scripts.
##########################################################################

MISMATCH_FILE="${STEP5_DIR}/id_list_sex_mismatch.txt"
ANEUPLOIDY_FILE="${STEP5_DIR}/id_list_aneuploidies.txt"

N_MISMATCH=0
N_ANEUPLOIDY=0

[ -f "${MISMATCH_FILE}" ]   && N_MISMATCH=$(wc -l < "${MISMATCH_FILE}" | tr -d '[:space:]')
[ -f "${ANEUPLOIDY_FILE}" ] && N_ANEUPLOIDY=$(wc -l < "${ANEUPLOIDY_FILE}" | tr -d '[:space:]')

echo "[CRISP] Sex mismatches flagged : ${N_MISMATCH}"
echo "[CRISP] Aneuploidies flagged   : ${N_ANEUPLOIDY}"
echo ""

EXCL_COMBINED="${STEP5_DIR}/${PROJECT_NAME}_step5_all_exclusions.txt"
> "${EXCL_COMBINED}"

CLEAN_PREFIX="${INPUT_PREFIX}"
N_EXCLUDED_MISMATCH=0
N_EXCLUDED_ANEUPLOIDY=0

if [ "${SEXCHECK_EXCLUDE_MISMATCH^^}" = "YES" ] && \
   [ -f "${MISMATCH_FILE}" ] && [ -s "${MISMATCH_FILE}" ]; then
    cat "${MISMATCH_FILE}" >> "${EXCL_COMBINED}"
    N_EXCLUDED_MISMATCH="${N_MISMATCH}"
    echo "[CRISP] Sex mismatch samples added to exclusion list : ${N_EXCLUDED_MISMATCH}"
else
    echo "[CRISP] SEXCHECK_EXCLUDE_MISMATCH = NO, sex mismatch samples retained."
fi

if [ "${SEXCHECK_EXCLUDE_ANEUPLOIDY^^}" = "YES" ] && \
   [ -f "${ANEUPLOIDY_FILE}" ] && [ -s "${ANEUPLOIDY_FILE}" ]; then
    cat "${ANEUPLOIDY_FILE}" >> "${EXCL_COMBINED}"
    N_EXCLUDED_ANEUPLOIDY="${N_ANEUPLOIDY}"
    echo "[CRISP] Aneuploidy samples added to exclusion list   : ${N_EXCLUDED_ANEUPLOIDY}"
else
    echo "[CRISP] SEXCHECK_EXCLUDE_ANEUPLOIDY = NO, aneuploidy samples retained."
fi

if [ -s "${EXCL_COMBINED}" ]; then
    sort -u "${EXCL_COMBINED}" -o "${EXCL_COMBINED}"
    N_TOTAL_EXCL=$(wc -l < "${EXCL_COMBINED}" | tr -d '[:space:]')
    echo ""
    echo "[CRISP] Applying exclusions to dataset..."

    CLEAN_PREFIX="${STEP5_DIR}/${PROJECT_NAME}_step5_clean"
    LOG_EXCL="${LOG_DIR}/step5_exclusions.log"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --remove "${EXCL_COMBINED}" \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_EXCL}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Sample exclusion failed. Check log: ${LOG_EXCL}"
        exit 1
    fi

    SAMPLES_AFTER=$(wc -l < "${CLEAN_PREFIX}.fam" | tr -d '[:space:]')
    echo "[OK]    Exclusions applied."
    echo "[CRISP] Samples before exclusion : ${SAMPLES}"
    echo "[CRISP] Samples excluded         : ${N_TOTAL_EXCL}"
    echo "[CRISP] Samples after exclusion  : ${SAMPLES_AFTER}"
    echo ""
else
    echo ""
    echo "[CRISP] No samples to exclude. Dataset unchanged."
    echo ""
    SAMPLES_AFTER="${SAMPLES}"
fi

##########################################################################
#
#   ____ ___  _   _  ___  ____ _____   ____   ___  _   _ _____ ___ _   _  ____
#  / ___/ _ \| | | |/ _ \|  _ \_   _| |  _ \ / _ \| | | |_   _|_ _| \ | |/ ___|
# | |  | | | | |_| | | | | |_) || |   | |_) | | | | | | | |  | ||  \| | |  _
# | |__| |_| |  _  | |_| |  _ < | |   |  _ <| |_| | |_| | |  | || |\  | |_| |
#  \____\___/|_| |_|\___/|_| \_\|_|   |_| \_\\___/ \___/|_| |___|_| \_|\____|
#
# Aneuploidies (Turner X0, Klinefelter XXY, Triple-X XXX) are confirmed
# karyotype calls, not data errors. When SEXCHECK_ROUTE_ANEUPLOIDIES = YES,
# each type is extracted into its own sub-cohort PLINK dataset — from BOTH
# the Step 4 QC'd data and the pre-QC raw genotypes — and compressed, ready
# for handoff to RNA/eQTL researchers via EHR linkage. Sex mismatches are
# never routed; they remain exclusion-only regardless of this flag.
##########################################################################

N_ROUTED_TOTAL=0
ANEUPLOIDY_TYPES=("turner" "klinefelter" "triple_x")
ANEUPLOIDY_TYPE_LABELS=("TURNER" "KLINEFELTER" "TRIPLE_X")

if [ "${SEXCHECK_ROUTE_ANEUPLOIDIES^^}" = "YES" ]; then

    echo "[CRISP] SEXCHECK_ROUTE_ANEUPLOIDIES = YES"
    echo "[CRISP] Routing aneuploidy samples to dedicated sub-cohorts..."
    echo "[CRISP] (These are confirmed karyotype calls, not exclusions.)"
    echo ""

    mkdir -p "${SEXCHECK_ANEUPLOIDY_COHORT_DIR}"

    MANIFEST_FILE="${STEP5_DIR}/aneuploidy_cohort_manifest.txt"

    if [ ! -f "${MANIFEST_FILE}" ]; then
        echo "[CRISP] WARNING: Aneuploidy manifest not found: ${MANIFEST_FILE}"
        echo "[CRISP]          Skipping cohort routing."
    else
        for i in "${!ANEUPLOIDY_TYPES[@]}"; do
            TYPE_KEY="${ANEUPLOIDY_TYPES[$i]}"
            TYPE_LABEL="${ANEUPLOIDY_TYPE_LABELS[$i]}"
            TYPE_KEEP_FILE="${STEP5_DIR}/id_list_aneuploidy_${TYPE_KEY}.txt"

            if [ ! -f "${TYPE_KEEP_FILE}" ] || [ ! -s "${TYPE_KEEP_FILE}" ]; then
                echo "[CRISP]   ${TYPE_LABEL}: no samples, skipping."
                continue
            fi

            N_TYPE=$(wc -l < "${TYPE_KEEP_FILE}" | tr -d '[:space:]')
            echo "[CRISP]   ${TYPE_LABEL}: ${N_TYPE} samples"

            TYPE_DIR="${SEXCHECK_ANEUPLOIDY_COHORT_DIR}/${TYPE_KEY}"
            mkdir -p "${TYPE_DIR}"

            # ── QC'd extraction (from Step 4 / current INPUT_PREFIX) ──────
            QC_OUT="${TYPE_DIR}/${PROJECT_NAME}_${TYPE_KEY}_qc"
            LOG_ROUTE_QC="${LOG_DIR}/step5_route_${TYPE_KEY}_qc.log"

            set +e
            ${PLINK1} \
                --bfile "${INPUT_PREFIX}" \
                --keep "${TYPE_KEEP_FILE}" \
                --make-bed \
                --out "${QC_OUT}" \
                >> "${LOG_ROUTE_QC}" 2>&1
            QC_STATUS=$?
            set -e

            if [ ${QC_STATUS} -ne 0 ]; then
                echo "[CRISP]   WARNING: QC'd extraction failed for ${TYPE_LABEL}."
                echo "[CRISP]            Check log: ${LOG_ROUTE_QC}"
            else
                echo "[CRISP]     QC'd cohort     : ${QC_OUT}.{bed,bim,fam}"
            fi

            # ── Pre-QC raw extraction ──────────────────────────────────────
            if [ -f "${PREQC_PREFIX}.bed" ]; then
                RAW_OUT="${TYPE_DIR}/${PROJECT_NAME}_${TYPE_KEY}_raw"
                LOG_ROUTE_RAW="${LOG_DIR}/step5_route_${TYPE_KEY}_raw.log"

                set +e
                ${PLINK1} \
                    --bfile "${PREQC_PREFIX}" \
                    --keep "${TYPE_KEEP_FILE}" \
                    --make-bed \
                    --out "${RAW_OUT}" \
                    >> "${LOG_ROUTE_RAW}" 2>&1
                RAW_STATUS=$?
                set -e

                if [ ${RAW_STATUS} -ne 0 ]; then
                    echo "[CRISP]   WARNING: Pre-QC extraction failed for ${TYPE_LABEL}."
                    echo "[CRISP]            Check log: ${LOG_ROUTE_RAW}"
                else
                    echo "[CRISP]     Pre-QC cohort   : ${RAW_OUT}.{bed,bim,fam}"
                fi
            else
                echo "[CRISP]   WARNING: Pre-QC source not found: ${PREQC_PREFIX}.bed"
                echo "[CRISP]            Skipping pre-QC extraction for ${TYPE_LABEL}."
            fi

            N_ROUTED_TOTAL=$((N_ROUTED_TOTAL + N_TYPE))
        done

        # Copy the typed manifest alongside the routed cohorts for handoff
        cp "${MANIFEST_FILE}" "${SEXCHECK_ANEUPLOIDY_COHORT_DIR}/aneuploidy_cohort_manifest.txt"

        echo ""
        echo "[CRISP] Compressing routed cohort directory..."
        ROUTE_ARCHIVE="${SEXCHECK_ANEUPLOIDY_COHORT_DIR}.tar.gz"
        tar -czf "${ROUTE_ARCHIVE}" \
            -C "$(dirname "${SEXCHECK_ANEUPLOIDY_COHORT_DIR}")" \
            "$(basename "${SEXCHECK_ANEUPLOIDY_COHORT_DIR}")"

        echo "[OK]    Aneuploidy cohort routing complete."
        echo "[CRISP] Total samples routed : ${N_ROUTED_TOTAL}"
        echo "[CRISP] Cohort directory     : ${SEXCHECK_ANEUPLOIDY_COHORT_DIR}"
        echo "[CRISP] Compressed archive   : ${ROUTE_ARCHIVE}"
        echo ""
    fi
else
    echo "[CRISP] SEXCHECK_ROUTE_ANEUPLOIDIES = NO, aneuploidy samples not routed."
    echo ""
fi

##########################################################################
#
#  ____  _______  __     _    __  __ _____ _   _ ____
# / ___|| ____\ \/ /    / \  |  \/  | ____| \ | |  _ \
# \___ \|  _|  \  /    / _ \ | |\/| |  _| |  \| | | | |
#  ___) | |___ /  \   / ___ \| |  | | |___| |\  | |_| |
# |____/|_____/_/\_\ /_/   \_\_|  |_|_____|_| \_|____/
#
# Optionally update PEDSEX in the FAM file to match genotype-inferred
# sex (SNPSEX). This is only appropriate after confirming with the data
# provider. The original FAM file is never overwritten; an amended copy
# and a full amendment log are written instead.
# Mutually exclusive with SEXCHECK_EXCLUDE_MISMATCH = YES (validated
# above at startup).
##########################################################################

AMENDMENT_LOG="${OUTPUT_DIR}/${PROJECT_NAME}_sex_amendments.txt"
N_AMENDED=0

if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ]; then

    echo "[CRISP WARNING] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "[CRISP WARNING] SEX AMENDMENT HAS BEEN REQUESTED."
    echo "[CRISP WARNING] The recorded sex in the FAM file will be updated to match"
    echo "[CRISP WARNING] the genotype-inferred sex for flagged samples."
    echo "[CRISP WARNING] Before proceeding, consult the data provider for each"
    echo "[CRISP WARNING] amended sample to confirm whether the recorded or"
    echo "[CRISP WARNING] genotyped sex is correct. Automatic amendment may not"
    echo "[CRISP WARNING] always be the right decision."
    echo "[CRISP WARNING] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""

    python3 - << PYEOF
import pandas as pd
import os
import sys

merged_file    = "${MERGED_FILE}"
input_fam      = "${INPUT_PREFIX}.fam"
amended_fam    = "${STEP5_DIR}/${PROJECT_NAME}_step5_amended.fam"
amendment_log  = "${AMENDMENT_LOG}"
count_sidecar  = "${STEP5_DIR}/.n_amended"

df_sex  = pd.read_csv(merged_file, sep=r'\s+')
df_fam  = pd.read_csv(input_fam, sep=r'\s+',
                       header=None,
                       names=['FID','IID','PAT','MAT','SEX','PHENO'])

mismatches = df_sex[
    (df_sex['PEDSEX'] != df_sex['SNPSEX']) &
    (df_sex['SNPSEX'].isin([1,2]))
][['FID','IID','PEDSEX','SNPSEX']]

if len(mismatches) == 0:
    print("[CRISP] No sex mismatches with valid SNPSEX found. No amendments made.")
    open(amendment_log,'w').close()
    with open(count_sidecar, 'w') as f:
        f.write("0\n")
    sys.exit(0)

with open(amendment_log, 'w') as log:
    log.write("=" * 64 + "\n")
    log.write("  CRISP: SEX AMENDMENT LOG\n")
    log.write("  Comprehensive Robust Integrated SNP Processing\n")
    log.write("=" * 64 + "\n\n")
    log.write("  WARNING: These samples had their recorded sex updated\n")
    log.write("  to match the genotype-inferred sex (SNPSEX).\n")
    log.write("  Consult the data provider before proceeding.\n\n")
    log.write(f"  {'FID':<12} {'IID':<12} {'Original SEX':<15} {'Amended SEX':<12}\n")
    log.write("-" * 64 + "\n")
    for _, row in mismatches.iterrows():
        orig = {1:'Male',2:'Female',0:'Unknown'}.get(int(row['PEDSEX']),'Unknown')
        new  = {1:'Male',2:'Female'}.get(int(row['SNPSEX']),'Unknown')
        log.write(f"  {row['FID']:<12} {row['IID']:<12} {orig:<15} {new:<12}\n")
        print(f"[CRISP WARNING] Amending: {row['FID']} {row['IID']} "
              f"SEX {int(row['PEDSEX'])} -> {int(row['SNPSEX'])}")
    log.write("\n" + "=" * 64 + "\n")

amendment_map = {
    (str(row['FID']), str(row['IID'])): int(row['SNPSEX'])
    for _, row in mismatches.iterrows()
}

df_fam['SEX'] = df_fam.apply(
    lambda r: amendment_map.get((str(r['FID']), str(r['IID'])), r['SEX']),
    axis=1
)

df_fam.to_csv(amended_fam, sep=' ', header=False, index=False)

# FIX v0.4.0: write count to a sidecar file so the calling shell script
# can read it back — heredocs run in a subshell and cannot directly set
# variables in the parent shell. Previously N_AMENDED stayed at its
# initialised 0 regardless of how many amendments were actually made.
with open(count_sidecar, 'w') as f:
    f.write(f"{len(mismatches)}\n")

print(f"\n[CRISP] Amended FAM written : {amended_fam}")
print(f"[CRISP] Amendment log       : {amendment_log}")
print(f"[CRISP] Samples amended     : {len(mismatches)}")
PYEOF

    if [ -f "${STEP5_DIR}/.n_amended" ]; then
        N_AMENDED=$(cat "${STEP5_DIR}/.n_amended" | tr -d '[:space:]')
        rm -f "${STEP5_DIR}/.n_amended"
    fi

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Sex amendment step failed."
        exit 1
    fi

    echo ""
    echo "[CRISP WARNING] Amendment log written to: ${AMENDMENT_LOG}"
    echo "[CRISP WARNING] Please review before continuing the pipeline."
    echo ""
fi

##########################################################################
#
#  ____  _____ ____   ___  ____ _____
# |  _ \| ____|  _ \ / _ \|  _ \_   _|
# | |_) |  _| | |_) | | | | |_) || |
# |  _ <| |___|  __/| |_| |  _ < | |
# |_| \_\_____|_|    \___/|_| \_\|_|
#
##########################################################################

{
    echo "=================================================================="
    echo "  CRISP: STEP 5 SEX CHECK REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Input prefix : ${INPUT_PREFIX}"
    echo "------------------------------------------------------------------"
    echo "  PARAMETERS"
    echo "  SEX_F_FEMALE_MAX            : ${SEX_F_FEMALE_MAX}"
    echo "  SEX_F_MALE_MIN              : ${SEX_F_MALE_MIN}"
    echo "  SEX_F_TURNER                : ${SEX_F_TURNER}"
    echo "  SEX_F_KLINEFELTER           : ${SEX_F_KLINEFELTER}"
    echo "  SEX_F_XXX                   : ${SEX_F_XXX}"
    echo "  SEX_Y_USE_MEAN              : ${SEX_Y_USE_MEAN}"
    if [ "${SEX_Y_USE_MEAN^^}" = "NO" ]; then
    echo "  SEX_Y_MANUAL                : ${SEX_Y_MANUAL}"
    fi
    echo "  PLOT_COLOUR_MODE            : ${PLOT_COLOUR_MODE}"
    echo "  PLOT_BACKGROUND              : ${CRISP_PAL_BACKGROUND}"
    echo "------------------------------------------------------------------"
    echo "  DETECTION RESULTS"
    echo "  Total samples               : ${SAMPLES}"
    echo "  Sex mismatches flagged      : ${N_MISMATCH}"
    echo "  Aneuploidies flagged        : ${N_ANEUPLOIDY}"
    echo "------------------------------------------------------------------"
    echo "  EXCLUSIONS"
    echo "  SEXCHECK_EXCLUDE_MISMATCH   : ${SEXCHECK_EXCLUDE_MISMATCH}"
    echo "  SEXCHECK_EXCLUDE_ANEUPLOIDY : ${SEXCHECK_EXCLUDE_ANEUPLOIDY}"
    if [ -s "${EXCL_COMBINED}" ]; then
    echo "  Samples excluded            : ${N_TOTAL_EXCL}"
    echo "  Samples after exclusion     : ${SAMPLES_AFTER}"
    echo "  Clean dataset prefix        : ${CLEAN_PREFIX}"
    else
    echo "  Samples excluded            : 0"
    fi
    echo "------------------------------------------------------------------"
    echo "  ANEUPLOIDY COHORT ROUTING"
    echo "  SEXCHECK_ROUTE_ANEUPLOIDIES : ${SEXCHECK_ROUTE_ANEUPLOIDIES}"
    if [ "${SEXCHECK_ROUTE_ANEUPLOIDIES^^}" = "YES" ]; then
    echo "  Samples routed               : ${N_ROUTED_TOTAL}"
    echo "  Cohort directory             : ${SEXCHECK_ANEUPLOIDY_COHORT_DIR}"
    echo "  Compressed archive           : ${SEXCHECK_ANEUPLOIDY_COHORT_DIR}.tar.gz"
    echo "  Manifest                     : ${STEP5_DIR}/aneuploidy_cohort_manifest.txt"
    echo ""
    echo "  NOTE: Aneuploidies routed here are confirmed karyotype calls,"
    echo "  not data errors — they remain valuable for RNA/eQTL research"
    echo "  as natural dosage experiments. Both Step 4 QC'd and pre-QC raw"
    echo "  genotypes were extracted per type for downstream flexibility."
    fi
    echo "------------------------------------------------------------------"
    echo "  SEX AMENDMENT"
    echo "  SEXCHECK_AMEND_SEX          : ${SEXCHECK_AMEND_SEX}"
    if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ]; then
    echo ""
    echo "  WARNING: Sex amendment was applied. Please review:"
    echo "  ${AMENDMENT_LOG}"
    echo "  Samples amended : ${N_AMENDED}"
    echo "  before continuing the pipeline."
    fi
    echo "------------------------------------------------------------------"
    echo "  OUTPUT FILES"
    echo "  Plots       : ${STEP5_DIR}/Sex_check.pdf"
    echo "  Full report : ${STEP5_DIR}/report_sex.mismatch_aneuploidies.txt"
    echo "  Mismatches  : ${MISMATCH_FILE}"
    echo "  Aneuploidies: ${ANEUPLOIDY_FILE}"
    if [ -s "${EXCL_COMBINED}" ]; then
    echo "  Excl. list  : ${EXCL_COMBINED}"
    fi
    if [ "${SEXCHECK_ROUTE_ANEUPLOIDIES^^}" = "YES" ]; then
    echo "  Aneu. cohort: ${SEXCHECK_ANEUPLOIDY_COHORT_DIR}.tar.gz"
    fi
    if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ]; then
    echo "  Amendments  : ${AMENDMENT_LOG}"
    fi
    echo "  Log         : ${LOG_DIR}"
    echo "=================================================================="
    echo "  END OF REPORT"
    echo "=================================================================="
} | tee "${REPORT_FILE}"

echo ""
echo "[CRISP] Report written to: ${REPORT_FILE}"
echo ""
echo "[CRISP] Step 5 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
