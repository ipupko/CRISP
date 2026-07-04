#!/bin/bash
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
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Variant Call Rate (Step 4)
### Version: 0.3.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
#
# WHAT THIS SCRIPT DOES
# ----------------------
# Step 4 is the variant-level quality gate. After sample-level call rate
# filtering in Step 3, we now filter on the variant axis: a SNP that is
# missing across too many samples is unreliable regardless of per-sample
# quality. Step 4 also removes monomorphic SNPs and optionally applies a
# minor allele frequency threshold, producing a clean, analysis-ready
# variant set for downstream GWAS.
#
# Sub-steps:
#   4a : Variant call rate filtering (SIMPLE / CASCADE / CUSTOM)
#   4b : Diagnostic frequency counts + MAF spectrum (diagnostic checkpoint)
#   4c : Monomorphic SNP removal   (REMOVE_MONO = YES / NO)
#   4d : MAF filtering             (FILTER_MAF  = YES / NO)
#   4e : Final MAF spectrum        (only if 4c/4d changed the variant set)
#
# WHY THREE CALL-RATE MODES?
# ---------------------------
# SIMPLE applies a single --geno threshold in one pass. This is the right
# choice for well-curated arrays where missingness is uniformly distributed.
#
# CASCADE applies four progressive passes (0.25 → 0.20 → 0.10 → 0.05).
# Progressive filtering avoids the pathological case where a single lenient
# pass retains a cluster of highly-correlated, high-missingness SNPs that
# would all survive together. Each pass removes the worst offenders before
# the next threshold tightens, producing a cleaner result than a single
# aggressive threshold.
#
# CUSTOM lets the user define their own tier sequence when the fixed CASCADE
# thresholds do not match their data's missingness profile.
#
# WHY TWO MAF SPECTRUM CHECKPOINTS?
# -----------------------------------
# The diagnostic checkpoint (post-4a) shows the MAF distribution of the
# surviving variant set before any monomorphic or MAF-based removal. The
# final checkpoint (post-4c/4d, only if those steps ran) shows the cleaned
# distribution. Comparing the two reveals how much of the original spectrum
# was contributed by very-low-frequency variants that were subsequently
# removed, which is critical context for interpreting downstream GWAS power.
#
# Output:
#   Filtered BED/BIM/FAM per pass
#   .lmiss per-variant missingness statistics
#   PDF plots via plot_snprate.R or plot_snprate.py
#   MAF spectrum text table and PDF (maf_spectrum.R)
#   Structured Step 4 report
#   Separate exclusion lists per sub-step
#
# Usage:
#   bash crisp_snprate.sh
#   bash crisp_snprate.sh --config my_project.txt
##########################################################################

# Exit immediately on any error, on use of unset variables, and propagate
# pipe failures. This prevents silent partial runs.
set -euo pipefail

# ####################################
#  _____ _                  ___
# /  ___| |                /   |
# \ `--.| |_ ___ _ __     / /| |
#  `--. \ __/ _ \ '_ \   / /_| |
# /\__/ / ||  __/ |_) |  \___  |
# \____/ \__\___| .__/       |_/
#               | |
#               |_|
# ####################################
##########################################################################
### STEP 4 : LOAD CRISP FLAVOUR
# WHY: _crisp_flavour.sh is the shared runtime initialiser for the
# entire pipeline. It exports colour palette variables (CRISP_PAL_*),
# loading messages, and the _init_runtime / _crisp_palette helpers.
# Every CRISP script must source it first -- before any echo, any
# parse_param call, and before set -euo pipefail can interfere with
# sourcing.
##########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime

# #############################################
#  _____      _ _   _       _ _
# |_   _|    (_) | (_)     | (_)
#   | | _ __  _| |_ _  __ _| |_ ___  ___
#   | || '_ \| | __| |/ _` | | / __|/ _ \
#  _| || | | | | |_| | (_| | | \__ \  __/
#  \___/_| |_|_|\__|_|\__,_|_|_|___/\___|
# #############################################
##########################################################################
### STEP 4 : INITIALISE -- PARSE COMMAND LINE ARGUMENTS
# WHY: We accept an optional --config flag so that multiple projects can
# share the same script with different instruction files. The default
# crisp_instructions.txt serves single-project setups without any flags.
##########################################################################

INSTRUCTION_FILE="crisp_instructions.txt"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            INSTRUCTION_FILE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: bash crisp_snprate.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_snprate.sh [--config <file>]"
            exit 1
            ;;
    esac
done

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_snprate.sh --config <path>"
    exit 1
fi

# ###########################################################
# ______                              _
# | ___ \                            | |
# | |_/ /_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___
# |  __/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | | | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# \_|  \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
# ###########################################################
##########################################################################
### STEP 4 : PARAMETERS -- PARSE INSTRUCTION FILE
# WHY: parse_param() is CRISP's standard key=value reader. It is grep-
# based rather than sourcing the instruction file directly. Sourcing would
# execute any shell syntax in the file, which is a security risk. grep +
# sed extracts values safely, ignoring comments and whitespace.
#
# Defaults are set conservatively:
#   GENO       = 0.05  : 5% missingness is the standard array-QC threshold
#   CALLRATE_MODE = SIMPLE : safest default; CASCADE / CUSTOM are opt-in
#   REMOVE_MONO   = YES : monomorphic SNPs carry no information for GWAS
#   FILTER_MAF    = YES : variants below 1% MAF are unreliable in typical
#                         array cohorts (low power, high genotyping error)
#   MAF           = 0.01: the conventional GWAS lower bound for common-
#                         variant analysis; see Step 4 docs for rare-variant
#                         guidance (separate pipeline path, future release)
#   PLOT_COLOUR_MODE = STANDARD: activates the default CRISP palette;
#                         COLOURBLIND switches to Okabe-Ito (Nature Methods
#                         recommended); NIGHT enables a dark-background theme
##########################################################################

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
CALLRATE_MODE=$(parse_param "CALLRATE_MODE" "SIMPLE")
GENO=$(parse_param "GENO" "0.05")
VARIANT_CUSTOM_TIERS=$(parse_param "VARIANT_CUSTOM_TIERS" "")
REMOVE_MONO=$(parse_param "REMOVE_MONO" "YES")
FILTER_MAF=$(parse_param "FILTER_MAF" "YES")
MAF=$(parse_param "MAF" "0.01")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")
PLOT_COLOUR_MODE=$(parse_param "PLOT_COLOUR_MODE" "STANDARD")

# MAF spectrum parameters (v0.3.0)
# WHY COMMON/CUSTOM only: RARE mode requires MAC-based binning which
# depends on cohort size and is being researched separately. It will be
# added in a future session once the bin design is validated.
MAF_SPECTRUM_MODE=$(parse_param "MAF_SPECTRUM_MODE" "COMMON")
MAF_SPECTRUM_BINS=$(parse_param "MAF_SPECTRUM_BINS" "0,0.0001,0.01,0.05,0.10,1")

# Export CRISP_PAL_* environment variables for all plot scripts in this
# session. This must happen before any plot script is called.
# WHY: R and Python plot scripts cannot source _crisp_flavour.sh directly
# (it is bash), so they read CRISP_PAL_* from the environment instead.
_crisp_palette "${PLOT_COLOUR_MODE}"

# Resolve the input prefix from the previous step.
# Priority: explicit STEP3_FINAL_PREFIX > CONVERTED_PREFIX > naming convention.
# WHY: Users may have renamed outputs or run steps out of order. Explicit
# overrides prevent silent fallback to the wrong file.
STEP3_FINAL=$(parse_param "STEP3_FINAL_PREFIX" "")
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")

if [ -n "${STEP3_FINAL}" ] && [ -f "${STEP3_FINAL}.bed" ]; then
    INPUT_PREFIX="${STEP3_FINAL}"
    echo "[CRISP] Using Step 3 output as input."
elif [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    INPUT_PREFIX="${CONVERTED_PREFIX}"
    echo "[CRISP] Step 3 output not found, using Step 2 output."
else
    INPUT_PREFIX="${OUTPUT_DIR}/step3_callrate/${PROJECT_NAME}_mind${GENO}"
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
    fi
fi

echo ""
echo "[CRISP] ${PREP_MSG}"
echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""
echo "[CRISP] Parameters:"
echo "        CALLRATE_MODE        : ${CALLRATE_MODE}"
echo "        GENO                 : ${GENO}"
if [ "${CALLRATE_MODE^^}" = "CUSTOM" ]; then
echo "        VARIANT_CUSTOM_TIERS : ${VARIANT_CUSTOM_TIERS}"
fi
echo "        REMOVE_MONO          : ${REMOVE_MONO}"
echo "        FILTER_MAF           : ${FILTER_MAF}  (MAF = ${MAF})"
echo "        INPUT_PREFIX         : ${INPUT_PREFIX}"
echo "        PLOT_ENGINE          : ${PLOT_ENGINE}"
echo "        PLOT_COLOUR_MODE     : ${PLOT_COLOUR_MODE}"
echo "        MAF_SPECTRUM_MODE    : ${MAF_SPECTRUM_MODE}"
echo "        MAF_SPECTRUM_BINS    : ${MAF_SPECTRUM_BINS}"
echo ""

# ###########################################
#  _   _       _ _     _       _
# | | | |     | (_)   | |     | |
# | | | | __ _| |_  __| | __ _| |_ ___
# | | | |/ _` | | |/ _` |/ _` | __/ _ \
# \ \_/ / (_| | | | (_| | (_| | ||  __/
#  \___/ \__,_|_|_|\__,_|\__,_|\__\___|
# ###########################################
##########################################################################
### STEP 4 : VALIDATE INPUT FILES
# WHY: We verify .bed/.bim/.fam before doing anything else. PLINK will
# produce a confusing error message if any of the three is missing. A
# clear early failure with an actionable message is better than a cryptic
# PLINK segfault or empty output file mid-run.
##########################################################################

for ext in .bed .bim .fam; do
    if [ ! -f "${INPUT_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${INPUT_PREFIX}${ext}"
        echo "[CRISP]        Ensure Step 3 (crisp_callrate.sh) has been run first."
        exit 1
    fi
done

echo "[OK]    Input files verified."
echo ""

# Count variants and samples up front so we can report them throughout
# and compute exclusion counts without re-scanning later.
VARIANTS_BEFORE=$(wc -l < "${INPUT_PREFIX}.bim" | tr -d '[:space:]')
SAMPLES=$(wc -l < "${INPUT_PREFIX}.fam" | tr -d '[:space:]')

echo "[CRISP] Variants before filtering : ${VARIANTS_BEFORE}"
echo "[CRISP] Samples (unchanged)       : ${SAMPLES}"
echo ""

# Set up output directories. LOG_DIR keeps PLINK logs separate from data
# outputs so the results directory stays clean.
STEP4_DIR="${OUTPUT_DIR}/step4_snprate"
LOG_DIR="${STEP4_DIR}/logs"
REPORT_FILE="${STEP4_DIR}/${PROJECT_NAME}_step4_report.txt"

mkdir -p "${STEP4_DIR}" "${LOG_DIR}"

# Temporary working directory -- cleaned up on EXIT via trap regardless
# of whether the script succeeds or fails.
CRISP_TMP="${OUTPUT_DIR}/.crisp_tmp_$$"
mkdir -p "${CRISP_TMP}"
trap 'rm -rf "${CRISP_TMP}"' EXIT

##########################################################################
### INTERNAL FUNCTIONS
##########################################################################

# write_report: tee output to both stdout and the report file.
# WHY: Users need a persistent record of what was done to their data.
# Plain text is preferred over binary formats for reproducibility --
# anyone can read it without special software.
write_report() {
    echo "$@" | tee -a "${REPORT_FILE}"
}

# validate_tiers: confirm CUSTOM tier list is valid before any PLINK
# calls are made, so we fail fast rather than partway through.
# WHY: A non-descending tier sequence would silently undo prior passes;
# tiers outside (0,1) are PLINK errors. Both conditions are caught here.
validate_tiers() {
    local tiers_str="$1"
    IFS=',' read -ra tiers <<< "${tiers_str}"
    local prev=1
    for t in "${tiers[@]}"; do
        if ! [[ "${t}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
            echo "[CRISP] ERROR: Non-numeric tier value: ${t}"
            exit 1
        fi
        if (( $(echo "${t} >= ${prev}" | bc -l) )); then
            echo "[CRISP] ERROR: VARIANT_CUSTOM_TIERS must be strictly descending. Got ${t} after ${prev}."
            exit 1
        fi
        if (( $(echo "${t} <= 0 || ${t} >= 1" | bc -l) )); then
            echo "[CRISP] ERROR: Tier ${t} is outside (0, 1)."
            exit 1
        fi
        prev="${t}"
    done
}

# generate_lmiss: produce per-variant missingness statistics.
# WHY: .lmiss files are the input to the missingness plots. We generate
# them after each filtering pass so the plot shows the distribution at
# that exact threshold, not a residual from a prior pass.
generate_lmiss() {
    local bfile="$1"
    local out_prefix="$2"
    local log_file="$3"

    ${PLINK1} \
        --bfile "${bfile}" \
        --missing \
        --allow-extra-chr \
        --out "${out_prefix}" \
        >> "${log_file}" 2>&1
}

# generate_plots: dispatch to R or Python plot engine.
# WHY: Dual-engine parity (R via ggplot2, Python via matplotlib) ensures
# the pipeline works in environments where only one is available.
# PLOT_ENGINE is set in the instruction file and exported at startup.
generate_plots() {
    local mode="$1"
    local plot_dir="$2"
    local log_file="$3"
    shift 3
    local lmiss_files=("$@")

    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_snprate.R" \
            "${mode}" "${GENO}" "${plot_dir}" \
            "${lmiss_files[@]}" \
            >> "${log_file}" 2>&1
    elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
        python3 "${SCRIPT_DIR}/scripts/plot_snprate.py" \
            "${mode}" "${GENO}" "${plot_dir}" \
            "${lmiss_files[@]}" \
            >> "${log_file}" 2>&1
    fi
}

# generate_freq_table: run plink --freq on a bfile and return the path
# to the resulting .frq file via stdout.
# WHY: --allow-extra-chr is required across all PLINK calls in CRISP so
# that non-standard chromosome codes (e.g. MT, chrX_PAR) do not cause
# PLINK to abort. All other PLINK calls in this script include it too.
generate_freq_table() {
    local bfile="$1"
    local out_prefix="$2"
    local log_file="$3"

    ${PLINK1} \
        --bfile "${bfile}" \
        --freq \
        --allow-extra-chr \
        --out "${out_prefix}" \
        >> "${log_file}" 2>&1

    echo "${out_prefix}.frq"
}

# MAF spectrum checkpoint accumulator.
# queue_maf_spectrum adds a labelled checkpoint to MAF_SPECTRUM_ARGS[].
# run_maf_spectrum fires maf_spectrum.R with all queued checkpoints at once.
# WHY: Separating queuing from execution lets us accumulate both the
# diagnostic and final checkpoints before making a single R call.
# maf_spectrum.R then produces one grouped bar chart comparing both
# checkpoints, which is more useful than two separate single-bar charts.
MAF_SPECTRUM_ARGS=()

queue_maf_spectrum() {
    local checkpoint_label="$1"
    local freq_file="$2"

    if [ ! -f "${freq_file}" ]; then
        echo "[CRISP] WARNING: MAF spectrum freq file not found, skipping: ${freq_file}"
        return
    fi

    MAF_SPECTRUM_ARGS+=("${checkpoint_label}:${freq_file}")
    echo "[CRISP] MAF spectrum queued: ${checkpoint_label} -> $(basename ${freq_file})"
}

run_maf_spectrum() {
    if [ "${#MAF_SPECTRUM_ARGS[@]}" -eq 0 ]; then
        echo "[CRISP] WARNING: No MAF spectrum checkpoints queued, skipping."
        return
    fi

    local log_file="${LOG_DIR}/step4_maf_spectrum.log"

    echo ""
    echo "##########################################################################"
    echo "###              STEP 4e: MAF SPECTRUM                                ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Running MAF spectrum (mode: ${MAF_SPECTRUM_MODE})..."
    echo "[CRISP] Checkpoints: ${MAF_SPECTRUM_ARGS[*]}"
    echo ""

    # Pass MODE, BINS, SAMPLES, output dir, project name, then all
    # labelled checkpoint args (e.g. "diagnostic:/path/to.frq").
    # maf_spectrum.R produces both a text table and a PDF grouped bar
    # chart comparing all checkpoints.
    ${RSCRIPT} "${SCRIPT_DIR}/scripts/maf_spectrum.R" \
        "${MAF_SPECTRUM_MODE}" \
        "${MAF_SPECTRUM_BINS}" \
        "${SAMPLES}" \
        "${STEP4_DIR}" \
        "${PROJECT_NAME}" \
        "${MAF_SPECTRUM_ARGS[@]}" \
        >> "${log_file}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] WARNING: MAF spectrum script failed. Check: ${log_file}"
        echo "[CRISP]          Pipeline continues -- spectrum is non-critical."
    else
        echo "[OK]    MAF spectrum: ${STEP4_DIR}/${PROJECT_NAME}_step4_maf_spectrum.txt"
        echo "[OK]    MAF spectrum plot: ${STEP4_DIR}/${PROJECT_NAME}_step4_maf_spectrum.pdf"
    fi
    echo ""
}

# ##########################################
#  _____ _                  ___
# /  ___| |                /   |
# \ `--.| |_ ___ _ __     / /| | __ _
#  `--. \ __/ _ \ '_ \   / /_| |/ _` |
# /\__/ / ||  __/ |_) |  \___  | (_| |
# \____/ \__\___| .__/       |_/\__,_|
#               | |
#               |_|
# ##########################################
##########################################################################
### STEP 4a : VARIANT CALL RATE FILTERING
# WHY: --geno removes any variant with a missingness rate above the
# threshold. It operates on genotype-level data, not individual-level
# (that is --mind, done in Step 3). The two are complementary: Step 3
# removes samples that are globally unreliable; Step 4 removes variants
# that are locally unreliable regardless of which samples cause the
# missingness.
#
# --allow-extra-chr: required to handle non-standard chromosome codes.
# Without it, PLINK halts on MT, chrX_PAR, or cohort-specific codes.
#
# The exclusion list (comm -23 on sorted BIM variant IDs) is written
# so that downstream tools can audit exactly which variants were removed
# and at which threshold.
##########################################################################

# --- SIMPLE MODE ---
run_simple() {
    local log_file="${LOG_DIR}/step4a_simple.log"
    local out_prefix="${STEP4_DIR}/${PROJECT_NAME}_geno${GENO}"

    echo "[CRISP] CALLRATE_MODE = SIMPLE"
    echo "[CRISP] Applying --geno ${GENO} ..."

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --geno "${GENO}" \
        --make-bed \
        --allow-extra-chr \
        --out "${out_prefix}" \
        >> "${log_file}" 2>&1

    POST_CALLRATE_PREFIX="${out_prefix}"

    # Compute exclusion list by comparing BIM files before and after.
    # comm -23 outputs lines only in the first file (i.e. removed variants).
    local excl_list="${STEP4_DIR}/${PROJECT_NAME}_excl_step4a.txt"
    comm -23 \
        <(awk '{print $2}' "${INPUT_PREFIX}.bim"  | sort) \
        <(awk '{print $2}' "${out_prefix}.bim" | sort) \
        > "${excl_list}"

    VARIANTS_AFTER_4A=$(wc -l < "${out_prefix}.bim" | tr -d '[:space:]')
    local removed=$(( VARIANTS_BEFORE - VARIANTS_AFTER_4A ))

    echo "[OK]    SIMPLE pass complete."
    echo "[CRISP] Variants removed : ${removed}"
    echo "[CRISP] Variants retained: ${VARIANTS_AFTER_4A}"
    echo ""

    generate_lmiss "${out_prefix}" \
        "${STEP4_DIR}/${PROJECT_NAME}_simple" "${log_file}"
    generate_plots "SIMPLE" "${STEP4_DIR}" "${log_file}" \
        "${STEP4_DIR}/${PROJECT_NAME}_simple.lmiss"
}

# --- MULTIPASS (CASCADE / CUSTOM) ---
# WHY: Each pass uses the output of the previous pass as input. This
# progressive approach removes the worst-missingness variants first so
# that later, tighter thresholds see a cleaner dataset. Running a single
# tight threshold on the original data is less effective because high-
# missingness variants inflate the effective threshold for borderline ones.
run_multipass() {
    local mode="$1"
    shift
    local thresholds=("$@")

    local current_input="${INPUT_PREFIX}"
    local cumulative_excl="${STEP4_DIR}/${PROJECT_NAME}_excl_step4a.txt"
    > "${cumulative_excl}"

    local lmiss_files=()
    local pass_num=0

    for threshold in "${thresholds[@]}"; do
        (( pass_num++ ))
        local log_file="${LOG_DIR}/step4a_pass${pass_num}.log"
        local out_prefix="${STEP4_DIR}/${PROJECT_NAME}_pass${pass_num}_geno${threshold}"

        echo "[CRISP] Pass ${pass_num} / ${#thresholds[@]} : --geno ${threshold}"

        ${PLINK1} \
            --bfile "${current_input}" \
            --geno "${threshold}" \
            --make-bed \
            --allow-extra-chr \
            --out "${out_prefix}" \
            >> "${log_file}" 2>&1

        # Accumulate removed variants from this pass into the master list.
        comm -23 \
            <(awk '{print $2}' "${current_input}.bim" | sort) \
            <(awk '{print $2}' "${out_prefix}.bim"    | sort) \
            >> "${cumulative_excl}"

        local after=$(wc -l < "${out_prefix}.bim" | tr -d '[:space:]')
        local before_pass=$(wc -l < "${current_input}.bim" | tr -d '[:space:]')
        echo "[OK]    Pass ${pass_num}: ${before_pass} -> ${after} variants"

        generate_lmiss "${out_prefix}" \
            "${STEP4_DIR}/${PROJECT_NAME}_pass${pass_num}" "${log_file}"
        lmiss_files+=("${STEP4_DIR}/${PROJECT_NAME}_pass${pass_num}.lmiss")

        # Clean up intermediate BED/BIM/FAM from prior passes unless the
        # user has asked to keep them (KEEP_INTERMEDIATE = YES).
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && [ "${current_input}" != "${INPUT_PREFIX}" ]; then
            rm -f "${current_input}.bed" "${current_input}.bim" "${current_input}.fam"
        fi

        current_input="${out_prefix}"
    done

    POST_CALLRATE_PREFIX="${current_input}"
    VARIANTS_AFTER_4A=$(wc -l < "${current_input}.bim" | tr -d '[:space:]')

    echo ""
    echo "[OK]    ${mode} complete. ${VARIANTS_BEFORE} -> ${VARIANTS_AFTER_4A} variants."
    echo ""

    generate_plots "${mode}" "${STEP4_DIR}" \
        "${LOG_DIR}/step4a_plots.log" "${lmiss_files[@]}"
}

run_cascade() {
    echo "[CRISP] CALLRATE_MODE = CASCADE"
    run_multipass "CASCADE" "0.25" "0.20" "0.10" "0.05"
}

run_custom() {
    echo "[CRISP] CALLRATE_MODE = CUSTOM"
    validate_tiers "${VARIANT_CUSTOM_TIERS}"
    IFS=',' read -ra tiers <<< "${VARIANT_CUSTOM_TIERS}"
    run_multipass "CUSTOM" "${tiers[@]}"
}

# Dispatch based on CALLRATE_MODE.
case "${CALLRATE_MODE^^}" in
    SIMPLE)  run_simple  ;;
    CASCADE) run_cascade ;;
    CUSTOM)  run_custom  ;;
    *)
        echo "[CRISP] ERROR: Unknown CALLRATE_MODE '${CALLRATE_MODE}'. Use SIMPLE, CASCADE, or CUSTOM."
        exit 1
        ;;
esac

# ##########################################
#  _____ _                  ___ _
# /  ___| |                /   | |
# \ `--.| |_ ___ _ __     / /| | |__
#  `--. \ __/ _ \ '_ \   / /_| | '_ \
# /\__/ / ||  __/ |_) |  \___  | |_) |
# \____/ \__\___| .__/       |_/_.__/
#               | |
#               |_|
# ##########################################
##########################################################################
### STEP 4b : DIAGNOSTIC FREQUENCY COUNTS + MAF SPECTRUM (DIAGNOSTIC)
# WHY: After call-rate filtering but before monomorphic/MAF removal, we
# take a snapshot of the MAF distribution. This serves two purposes:
#
#   1. The diagnostic counts (monomorphic, below-MAF, passing) give the
#      user a preview of what Steps 4c and 4d will remove, so they can
#      judge whether their thresholds are appropriate before committing
#      to the removal.
#
#   2. The diagnostic MAF spectrum (via maf_spectrum.R) provides the
#      "before" side of the before/after grouped bar chart. Without this
#      checkpoint, there is no comparison -- Step 4e's final spectrum
#      would be a bar chart with nothing to compare against.
#
# The --freq output is also reused later if FINAL_PREFIX == POST_CALLRATE_PREFIX
# (i.e. 4c and 4d both skipped), in which case the diagnostic spectrum
# is the only spectrum and is presented as-is.
##########################################################################

DIAG_LOG="${LOG_DIR}/step4b_diag.log"

echo ""
echo "##########################################################################"
echo "###              STEP 4b: DIAGNOSTIC COUNTS                           ###"
echo "##########################################################################"
echo ""

${PLINK1} \
    --bfile "${POST_CALLRATE_PREFIX}" \
    --freq \
    --allow-extra-chr \
    --out "${STEP4_DIR}/${PROJECT_NAME}_diag_freq" \
    >> "${DIAG_LOG}" 2>&1

FREQ_FILE="${STEP4_DIR}/${PROJECT_NAME}_diag_freq.frq"

if [ -f "${FREQ_FILE}" ]; then
    DIAG_MONO=$(awk 'NR>1 && $5==0' "${FREQ_FILE}" | wc -l | tr -d '[:space:]')
    DIAG_BELOW_MAF=$(awk -v maf="${MAF}" 'NR>1 && $5<maf && $5>0' \
        "${FREQ_FILE}" | wc -l | tr -d '[:space:]')
    DIAG_TOTAL=$(awk 'NR>1' "${FREQ_FILE}" | wc -l | tr -d '[:space:]')
    DIAG_PASSING=$((DIAG_TOTAL - DIAG_MONO - DIAG_BELOW_MAF))
else
    DIAG_MONO="N/A"
    DIAG_BELOW_MAF="N/A"
    DIAG_TOTAL="N/A"
    DIAG_PASSING="N/A"
fi

echo "[CRISP] Diagnostic counts (post call-rate filter, pre-removal):"
echo "        Total variants             : ${DIAG_TOTAL}"
echo "        Monomorphic (MAF = 0)      : ${DIAG_MONO}"
echo "        Below MAF ${MAF}              : ${DIAG_BELOW_MAF}"
echo "        Passing both filters       : ${DIAG_PASSING}"
echo ""

# Queue the diagnostic MAF spectrum checkpoint.
# A separate --freq call is made here (distinct from the diag_freq.frq
# above) so that maf_spectrum.R always receives a file it can trust has
# the correct format, regardless of any future changes to Step 4b logic.
DIAG_SPEC_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_diag_spectrum"
DIAG_SPEC_FREQ=$(generate_freq_table "${POST_CALLRATE_PREFIX}" \
    "${DIAG_SPEC_PREFIX}" "${DIAG_LOG}")
queue_maf_spectrum "diagnostic" "${DIAG_SPEC_FREQ}"

# #########################################
#  _____ _                  ___
# /  ___| |                /   |
# \ `--.| |_ ___ _ __     / /| | ___
#  `--. \ __/ _ \ '_ \   / /_| |/ __|
# /\__/ / ||  __/ |_) |  \___  | (__
# \____/ \__\___| .__/       |_/\___|
#               | |
#               |_|
# #########################################
##########################################################################
### STEP 4c : MONOMORPHIC SNP REMOVAL
# WHY: A monomorphic SNP has MAF = 0, meaning every sample in the cohort
# carries the same allele. Variants with MAF very close to zero — while
# technically polymorphic — are similarly problematic: they carry near-zero
# variance, making association testing unreliable and genotyping error rates
# comparable to allele frequency. Such near-zero-MAF variants are handled
# by Step 4d's --maf threshold; Step 4c handles only the exact MAF==0 case.
#
# --mac 1 retains only variants where the minor allele appears at least
# once across all samples. This is equivalent to MAF > 0 but expressed
# as an exact integer count rather than a floating-point frequency,
# making it robust to sample size and free of floating-point imprecision.
#
# WHY NOT --maf 0.000001: Floating-point comparisons at near-zero MAF are
# imprecise in PLINK. --mac 1 is an exact integer comparison and always
# correct. Removing monomorphic SNPs explicitly is also considered good
# practice in genotype QC regardless of cohort size -- it makes the
# pipeline's behaviour deterministic and prevents silent tool failures.
##########################################################################

CURRENT_WORKING="${POST_CALLRATE_PREFIX}"
MONO_EXCL_LIST=""

if [ "${REMOVE_MONO^^}" = "YES" ]; then
    echo ""
    echo "##########################################################################"
    echo "###              STEP 4c: MONOMORPHIC REMOVAL                         ###"
    echo "##########################################################################"
    echo ""

    local_log="${LOG_DIR}/step4c_mono.log"
    mono_out="${STEP4_DIR}/${PROJECT_NAME}_nomono"
    MONO_EXCL_LIST="${STEP4_DIR}/${PROJECT_NAME}_excl_step4c_mono.txt"

    ${PLINK1} \
        --bfile "${CURRENT_WORKING}" \
        --mac 1 \
        --make-bed \
        --allow-extra-chr \
        --out "${mono_out}" \
        >> "${local_log}" 2>&1

    comm -23 \
        <(awk '{print $2}' "${CURRENT_WORKING}.bim" | sort) \
        <(awk '{print $2}' "${mono_out}.bim"        | sort) \
        > "${MONO_EXCL_LIST}"

    MONO_REMOVED=$(wc -l < "${MONO_EXCL_LIST}" | tr -d '[:space:]')
    AFTER_MONO=$(wc -l < "${mono_out}.bim" | tr -d '[:space:]')

    echo "[OK]    Monomorphic SNPs removed : ${MONO_REMOVED}"
    echo "[CRISP] Variants retained        : ${AFTER_MONO}"
    echo ""

    CURRENT_WORKING="${mono_out}"
else
    echo "[CRISP] REMOVE_MONO = NO, skipping."
    echo ""
fi

# ##########################################
#  _____ _                  ___     _
# /  ___| |                /   |   | |
# \ `--.| |_ ___ _ __     / /| | __| |
#  `--. \ __/ _ \ '_ \   / /_| |/ _` |
# /\__/ / ||  __/ |_) |  \___  | (_| |
# \____/ \__\___| .__/       |_/\__,_|
#               | |
#               |_|
# ##########################################
##########################################################################
### STEP 4d : MAF FILTERING
# WHY: Common-variant GWAS is conventionally limited to variants with
# MAF >= 1% (the default here). Below this threshold, genotyping error
# rates become comparable to allele frequency, making association signals
# unreliable. Additionally, statistical power for detecting associations
# at very low MAF requires sample sizes in the hundreds of thousands.
#
# The 1% threshold is well-established in the literature (see e.g. the
# UK Biobank QC paper, Bycroft et al. 2018). A 5% threshold is also
# commonly used -- particularly for smaller cohorts where variants in
# the 1-5% MAF range have insufficient power to detect associations
# anyway, or for studies where array genotyping accuracy in this range
# is uncertain. Users should set MAF = 0.05 in the instruction file
# when working with cohorts of N < 5,000 or when following a study-
# specific QC protocol that specifies 5%. The 1% default is appropriate
# for large biobank-scale datasets (N > 10,000).
#
# NOTE on rare variants: a separate rare-variant pipeline path is under
# development. The MAF threshold here is explicitly for common-variant
# arrays. Do not lower MAF below 0.001 in this step for rare-variant
# work -- use the dedicated path when available.
##########################################################################

MAF_EXCL_LIST=""
VARIANTS_AFTER_MAF=""

if [ "${FILTER_MAF^^}" = "YES" ]; then
    echo ""
    echo "##########################################################################"
    echo "###              STEP 4d: MAF FILTERING                               ###"
    echo "##########################################################################"
    echo ""

    local_log="${LOG_DIR}/step4d_maf.log"
    maf_out="${STEP4_DIR}/${PROJECT_NAME}_maf${MAF}"
    MAF_EXCL_LIST="${STEP4_DIR}/${PROJECT_NAME}_excl_step4d_maf.txt"

    ${PLINK1} \
        --bfile "${CURRENT_WORKING}" \
        --maf "${MAF}" \
        --make-bed \
        --allow-extra-chr \
        --out "${maf_out}" \
        >> "${local_log}" 2>&1

    comm -23 \
        <(awk '{print $2}' "${CURRENT_WORKING}.bim" | sort) \
        <(awk '{print $2}' "${maf_out}.bim"         | sort) \
        > "${MAF_EXCL_LIST}"

    MAF_REMOVED=$(wc -l < "${MAF_EXCL_LIST}" | tr -d '[:space:]')
    VARIANTS_AFTER_MAF=$(wc -l < "${maf_out}.bim" | tr -d '[:space:]')

    echo "[OK]    Variants removed by MAF < ${MAF} : ${MAF_REMOVED}"
    echo "[CRISP] Variants retained                : ${VARIANTS_AFTER_MAF}"
    echo ""

    FINAL_PREFIX="${maf_out}"
    CURRENT_WORKING="${maf_out}"
else
    echo "[CRISP] FILTER_MAF = NO, skipping."
    echo ""
    VARIANTS_AFTER_MAF="N/A (skipped)"
    FINAL_PREFIX="${CURRENT_WORKING}"
fi

# #########################################
#  _____ _                  ___
# /  ___| |                /   |
# \ `--.| |_ ___ _ __     / /| | ___
#  `--. \ __/ _ \ '_ \   / /_| |/ _ \
# /\__/ / ||  __/ |_) |  \___  |  __/
# \____/ \__\___| .__/       |_/\___|
#               | |
#               |_|
# #########################################
##########################################################################
### STEP 4e : FINAL MAF SPECTRUM
# WHY: We only run a second spectrum if the variant set actually changed
# during Steps 4c and 4d. If both steps were skipped (or neither removed
# any variants), the post-call-rate prefix equals the final prefix and
# the diagnostic spectrum already represents the final state -- generating
# a second identical chart would be redundant and misleading.
#
# When the final spectrum does run, run_maf_spectrum() fires maf_spectrum.R
# with both the diagnostic and final checkpoints, producing a single
# grouped bar chart that places them side by side for direct comparison.
# This is the primary visual output of Step 4 for MAF distribution review.
##########################################################################

if [ "${FINAL_PREFIX}" != "${POST_CALLRATE_PREFIX}" ]; then
    FINAL_SPEC_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_final_spectrum"
    FINAL_SPEC_FREQ=$(generate_freq_table "${FINAL_PREFIX}" \
        "${FINAL_SPEC_PREFIX}" \
        "${LOG_DIR}/step4e_final_spectrum.log")
    queue_maf_spectrum "final" "${FINAL_SPEC_FREQ}"
else
    echo "[CRISP] Post-call-rate prefix equals final prefix (4c/4d made no changes or were skipped)."
    echo "[CRISP] Diagnostic spectrum is the final spectrum -- no second checkpoint needed."
    echo ""
fi

run_maf_spectrum

# ######################################
# ______                      _
# | ___ \                    | |
# | |_/ /___ _ __   ___  _ __| |_
# |    // _ \ '_ \ / _ \| '__| __|
# | |\ \  __/ |_) | (_) | |  | |_
# \_| \_\___| .__/ \___/|_|   \__|
#           | |
#           |_|
# ######################################
##########################################################################
### STEP 4 : REPORT
# WHY: A structured plain-text report captures the key numbers from every
# sub-step in one place. This is the record a bioinformatician would
# consult when reviewing the QC trail for a publication or audit. It is
# written via write_report() which tees to stdout and the file
# simultaneously, so the user sees it in real time and also has a
# persistent copy.
##########################################################################

{
write_report ""
write_report "=================================================================="
write_report "  CRISP: STEP 4 VARIANT CALL RATE REPORT"
write_report "  Project       : ${PROJECT_NAME}"
write_report "  Instruction   : ${INSTRUCTION_FILE}"
write_report "  CALLRATE_MODE : ${CALLRATE_MODE}"
write_report "  GENO          : ${GENO}"
write_report "  REMOVE_MONO   : ${REMOVE_MONO}"
write_report "  FILTER_MAF    : ${FILTER_MAF} (MAF = ${MAF})"
write_report "=================================================================="
write_report ""
write_report "------------------------------------------------------------------"
write_report "  VARIANT COUNTS"
write_report "  Before Step 4      : ${VARIANTS_BEFORE}"
write_report "  After 4a (call rate): ${VARIANTS_AFTER_4A}"
if [ "${REMOVE_MONO^^}" = "YES" ]; then
write_report "  After 4c (mono rm) : ${AFTER_MONO}"
fi
if [ "${FILTER_MAF^^}" = "YES" ]; then
write_report "  After 4d (MAF filt): ${VARIANTS_AFTER_MAF}"
fi
write_report "  Samples (unchanged): ${SAMPLES}"
write_report "------------------------------------------------------------------"
write_report "  DIAGNOSTIC COUNTS (pre-removal)"
write_report "  Total variants     : ${DIAG_TOTAL}"
write_report "  Monomorphic        : ${DIAG_MONO}"
write_report "  Below MAF ${MAF}      : ${DIAG_BELOW_MAF}"
write_report "  Passing            : ${DIAG_PASSING}"
write_report "------------------------------------------------------------------"
write_report "  OUTPUT"
write_report "  Final BED prefix   : ${FINAL_PREFIX}"
write_report "  Call rate excl.    : ${STEP4_DIR}/${PROJECT_NAME}_excl_step4a.txt"
if [ "${REMOVE_MONO^^}" = "YES" ]; then
write_report "  Monomorphic excl.  : ${MONO_EXCL_LIST}"
fi
if [ "${FILTER_MAF^^}" = "YES" ]; then
write_report "  MAF excl.          : ${MAF_EXCL_LIST}"
fi
write_report "  MAF spectrum txt   : ${STEP4_DIR}/${PROJECT_NAME}_step4_maf_spectrum.txt"
write_report "  MAF spectrum pdf   : ${STEP4_DIR}/${PROJECT_NAME}_step4_maf_spectrum.pdf"
write_report "  Plots              : ${STEP4_DIR}"
write_report "  Log                : ${LOG_DIR}"
write_report "=================================================================="
write_report "  END OF REPORT"
write_report "=================================================================="
write_report ""
}

echo ""
echo "[OK]    Step 4 complete."
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
