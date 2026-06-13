#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Sample Call Rate (Step 3)
### Version: 0.3.5
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
### Date Updated : 13/06/2026
##########################################################################
# Filters samples based on genotype missingness rate using PLINK
# --mind flag.
#
# Three modes controlled by CALLRATE_MODE in crisp_instructions.txt:
#
#   SIMPLE  (default)
#     Single pass at MIND threshold (default 0.05)
#
#   CASCADE
#     Four fixed progressive passes:
#       Pass 1 : --mind 0.25
#       Pass 2 : --mind 0.20
#       Pass 3 : --mind 0.10
#       Pass 4 : --mind 0.05
#
#   CUSTOM
#     User-defined tiers via SAMPLE_CUSTOM_TIERS
#     Example: SAMPLE_CUSTOM_TIERS = 0.30,0.20,0.10,0.05
#     Any number of tiers, strictly descending, between 0 and 1
#
# Output:
#   .imiss              per-sample missingness statistics
#   PDF plots via plot_callrate.R or plot_callrate.py
#   Structured Step 3 report
#   Cumulative exclusion list for amendments step
#
# Usage:
#   bash crisp_callrate.sh
#   bash crisp_callrate.sh --config my_project.txt
##########################################################################

set -euo pipefail

##########################################################################
### LOAD CRISP FLAVOUR (NON-NEGOTIABLE)
##########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime

# Print a clear CRISP-style error if any command fails under set -e
trap 'echo "[CRISP] ERROR: Step 3 failed at line ${LINENO}. Check logs in ${LOG_DIR:-./results/logs}."' ERR

##########################################################################
### PARSE COMMAND LINE ARGUMENTS
##########################################################################

INSTRUCTION_FILE="crisp_instructions.txt"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            INSTRUCTION_FILE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: bash crisp_callrate.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_callrate.sh [--config <file>]"
            exit 1
            ;;
    esac
done

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_callrate.sh --config <path>"
    exit 1
fi

##########################################################################
### PARSE INSTRUCTION FILE
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
MIND=$(parse_param "MIND" "0.05")
SAMPLE_CUSTOM_TIERS=$(parse_param "SAMPLE_CUSTOM_TIERS" "")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PYTHON3=$(parse_param "PYTHON3_PATH" "python3")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")
PLOT_COLOUR_MODE=$(parse_param "PLOT_COLOUR_MODE" "STANDARD")
PLOT_BACKGROUND=$(parse_param "PLOT_BACKGROUND" "")

# resolve full CRISP_PAL_* palette (handles NIGHT -> dark background default)
_crisp_palette "${PLOT_COLOUR_MODE}" "${PLOT_BACKGROUND}"
PLOT_BACKGROUND="${CRISP_PAL_BACKGROUND}"

# use Step 2 output, fall back to naming convention if not set
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")
if [ -z "${CONVERTED_PREFIX}" ]; then
    CONVERTED_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
fi

##########################################################################
### HEADER
##########################################################################

echo ""
echo "[CRISP] ${PREP_MSG}"
echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""
echo "[CRISP] Parameters:"
echo "        CALLRATE_MODE     : ${CALLRATE_MODE}"
echo "        MIND              : ${MIND}"
if [ "${CALLRATE_MODE^^}" = "CUSTOM" ]; then
echo "        SAMPLE_CUSTOM_TIERS: ${SAMPLE_CUSTOM_TIERS}"
fi
echo "        INPUT_PREFIX      : ${CONVERTED_PREFIX}"
echo "        OUTPUT_DIR        : ${OUTPUT_DIR}"
echo "        PROJECT_NAME      : ${PROJECT_NAME}"
echo "        KEEP_INTERMEDIATE : ${KEEP_INTERMEDIATE}"
echo "        PLOT_ENGINE       : ${PLOT_ENGINE}"
echo "        PLOT_COLOUR_MODE  : ${CRISP_PAL_MODE}"
echo "        PLOT_BACKGROUND   : ${PLOT_BACKGROUND}"
echo ""

##########################################################################
### VALIDATE INPUT
##########################################################################

for ext in .bed .bim .fam; do
    if [ ! -f "${CONVERTED_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${CONVERTED_PREFIX}${ext}"
        echo "[CRISP]        Ensure Step 2 (crisp_convert.sh) has been run first."
        exit 1
    fi
done

echo "[OK]    Input files verified."
echo ""

##########################################################################
### SET UP DIRECTORIES
##########################################################################

STEP3_DIR="${OUTPUT_DIR}/step3_callrate"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step3_callrate_report.txt"

mkdir -p "${STEP3_DIR}"
mkdir -p "${LOG_DIR}"

SAMPLES_BEFORE=$(wc -l < "${CONVERTED_PREFIX}.fam" | tr -d '[:space:]')
VARIANTS=$(wc -l < "${CONVERTED_PREFIX}.bim" | tr -d '[:space:]')

echo "[CRISP] Samples before filtering : ${SAMPLES_BEFORE}"
echo "[CRISP] Variants                 : ${VARIANTS}"
echo ""

##########################################################################
### VALIDATE CUSTOM TIERS
##########################################################################

validate_tiers() {
    local tiers_str="$1"
    local prev=1.0

    IFS=',' read -ra tiers <<< "${tiers_str}"

    if [ "${#tiers[@]}" -lt 2 ]; then
        echo "[CRISP] ERROR: CUSTOM mode requires at least 2 tiers."
        echo "[CRISP]        Got: ${tiers_str}"
        exit 1
    fi

    for tier in "${tiers[@]}"; do
        if ! [[ "${tier}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
            echo "[CRISP] ERROR: Invalid tier '${tier}', must be numeric between 0 and 1."
            exit 1
        fi
        if (( $(echo "${tier} <= 0 || ${tier} >= 1" | bc -l) )); then
            echo "[CRISP] ERROR: Tier '${tier}' out of range, must be between 0 and 1."
            exit 1
        fi
        if (( $(echo "${tier} >= ${prev}" | bc -l) )); then
            echo "[CRISP] ERROR: Tiers must be strictly descending."
            echo "[CRISP]        Got: ${tiers_str}"
            echo "[CRISP]        Example: SAMPLE_CUSTOM_TIERS = 0.30,0.20,0.10,0.05"
            exit 1
        fi
        prev="${tier}"
    done

    echo "[OK]    Custom tiers validated: ${tiers_str}"
}

##########################################################################
### GENERATE IMISS FILE
##########################################################################

generate_imiss() {
    local bfile="$1"
    local out_prefix="$2"
    local log_file="$3"

    ${PLINK1} \
        --bfile "${bfile}" \
        --missing \
        --out "${out_prefix}" \
        >> "${log_file}" 2>&1
}

##########################################################################
### GENERATE PLOTS
##########################################################################

generate_plots() {
    local mode="$1"
    local plot_dir="$2"
    local log_file="$3"
    shift 3
    local imiss_files=("$@")

    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        CRISP_PAL_PASS="${CRISP_PAL_PASS}" \
        CRISP_PAL_FAIL="${CRISP_PAL_FAIL}" \
        CRISP_PAL_BG="${CRISP_PAL_BG}" \
        CRISP_PAL_PANEL="${CRISP_PAL_PANEL}" \
        CRISP_PAL_TEXT="${CRISP_PAL_TEXT}" \
        CRISP_PAL_SUBTEXT="${CRISP_PAL_SUBTEXT}" \
        CRISP_PAL_GRID="${CRISP_PAL_GRID}" \
        CRISP_PAL_MODE="${CRISP_PAL_MODE}" \
        CRISP_PAL_BACKGROUND="${CRISP_PAL_BACKGROUND}" \
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_callrate.R" \
            "${mode}" "${MIND}" "${plot_dir}" \
            "${imiss_files[@]}" \
            >> "${log_file}" 2>&1
    elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
        CRISP_PAL_PASS="${CRISP_PAL_PASS}" \
        CRISP_PAL_FAIL="${CRISP_PAL_FAIL}" \
        CRISP_PAL_BG="${CRISP_PAL_BG}" \
        CRISP_PAL_PANEL="${CRISP_PAL_PANEL}" \
        CRISP_PAL_TEXT="${CRISP_PAL_TEXT}" \
        CRISP_PAL_SUBTEXT="${CRISP_PAL_SUBTEXT}" \
        CRISP_PAL_GRID="${CRISP_PAL_GRID}" \
        CRISP_PAL_MODE="${CRISP_PAL_MODE}" \
        CRISP_PAL_BACKGROUND="${CRISP_PAL_BACKGROUND}" \
        ${PYTHON3} "${SCRIPT_DIR}/scripts/plot_callrate.py" \
            "${mode}" "${MIND}" "${plot_dir}" \
            "${imiss_files[@]}" \
            >> "${log_file}" 2>&1
    fi
}

##########################################################################
### SIMPLE MODE
##########################################################################

run_simple() {

    echo "##########################################################################"
    echo "###            STEP 3: SAMPLE CALL RATE (SIMPLE MODE)                 ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Applying --mind ${MIND}..."
    echo ""

    OUT_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_mind${MIND}"
    LOG_FILE="${LOG_DIR}/step3_simple.log"

    ${PLINK1} \
        --bfile "${CONVERTED_PREFIX}" \
        --mind "${MIND}" \
        --make-bed \
        --out "${OUT_PREFIX}" \
        >> "${LOG_FILE}" 2>&1

    SAMPLES_AFTER=$(wc -l < "${OUT_PREFIX}.fam" | tr -d '[:space:]')
    SAMPLES_REMOVED=$((SAMPLES_BEFORE - SAMPLES_AFTER))

    echo "[OK]    SIMPLE filtering complete."
    echo ""
    echo "[CRISP] Samples before : ${SAMPLES_BEFORE}"
    echo "[CRISP] Samples removed: ${SAMPLES_REMOVED}"
    echo "[CRISP] Samples after  : ${SAMPLES_AFTER}"
    echo ""

    # missingness stats for plotting
    IMISS_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_simple_imiss"
    generate_imiss "${CONVERTED_PREFIX}" "${IMISS_PREFIX}" "${LOG_FILE}"
    IMISS_FILE="${IMISS_PREFIX}.imiss"

    echo "[CRISP] Generating call rate plots..."
    generate_plots "SIMPLE" "${STEP3_DIR}" "${LOG_DIR}/step3_plots.log" "${IMISS_FILE}"
    echo "[OK]    Plots saved to: ${STEP3_DIR}"
    echo ""

    # exclusion list
    # Format matches PLINK .irem: whitespace-separated FID IID per line,
    # no header. An empty file means zero samples were excluded.
    EXCLUSION_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step3.txt"
    IREM_FILE="${OUT_PREFIX}.irem"
    if [ -f "${IREM_FILE}" ] && [ -s "${IREM_FILE}" ]; then
        cp "${IREM_FILE}" "${EXCLUSION_LIST}"
    else
        touch "${EXCLUSION_LIST}"
    fi
    echo "[CRISP] Exclusion list: ${EXCLUSION_LIST}"
    echo ""

    write_report "SIMPLE" "${MIND}" \
        "${SAMPLES_BEFORE}" "${SAMPLES_AFTER}" "${SAMPLES_REMOVED}" \
        "" "" "" "" "${OUT_PREFIX}"
}

##########################################################################
### MULTI-PASS RUNNER (shared by CASCADE and CUSTOM)
##########################################################################

run_multipass() {

    local mode="$1"
    shift
    local thresholds=("$@")
    local n="${#thresholds[@]}"

    CURRENT_INPUT="${CONVERTED_PREFIX}"
    IMISS_FILES=()
    TOTAL_REMOVED=0
    REMOVED_PER_PASS=()
    SAMPLES_PER_PASS=()

    CUMULATIVE_EXCL="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step3.txt"
    > "${CUMULATIVE_EXCL}"

    for i in "${!thresholds[@]}"; do

        PASS=$((i + 1))
        THRESHOLD="${thresholds[$i]}"
        OUT_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_mind${THRESHOLD}"
        LOG_FILE="${LOG_DIR}/step3_${mode,,}_pass${PASS}.log"

        echo "------------------------------------------------------------------"
        echo "[CRISP] Pass ${PASS} of ${n}, --mind ${THRESHOLD}"
        echo "------------------------------------------------------------------"
        echo ""

        SAMPLES_THIS_PASS=$(wc -l < "${CURRENT_INPUT}.fam" | tr -d '[:space:]')

        ${PLINK1} \
            --bfile "${CURRENT_INPUT}" \
            --mind "${THRESHOLD}" \
            --make-bed \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        SAMPLES_AFTER_PASS=$(wc -l < "${OUT_PREFIX}.fam" | tr -d '[:space:]')
        REMOVED_THIS_PASS=$((SAMPLES_THIS_PASS - SAMPLES_AFTER_PASS))
        TOTAL_REMOVED=$((TOTAL_REMOVED + REMOVED_THIS_PASS))

        REMOVED_PER_PASS+=("${REMOVED_THIS_PASS}")
        SAMPLES_PER_PASS+=("${SAMPLES_AFTER_PASS}")

        echo "[OK]    Pass ${PASS} complete."
        echo "        Samples this pass : ${SAMPLES_THIS_PASS}"
        echo "        Samples removed   : ${REMOVED_THIS_PASS}"
        echo "        Samples remaining : ${SAMPLES_AFTER_PASS}"
        echo ""

        # append to cumulative exclusion list
        IREM_FILE="${OUT_PREFIX}.irem"
        if [ -f "${IREM_FILE}" ] && [ -s "${IREM_FILE}" ]; then
            cat "${IREM_FILE}" >> "${CUMULATIVE_EXCL}"
        fi

        # imiss for plotting
        IMISS_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_mind${THRESHOLD}_imiss"
        generate_imiss "${CURRENT_INPUT}" "${IMISS_PREFIX}" "${LOG_FILE}"
        IMISS_FILES+=("${IMISS_PREFIX}.imiss")

        # clean up intermediate passes if not keeping
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && [ "${PASS}" -lt "${n}" ]; then
            for ext in .bed .bim .fam .log .nosex; do
                [ -f "${OUT_PREFIX}${ext}" ] && rm "${OUT_PREFIX}${ext}"
            done
        fi

        CURRENT_INPUT="${OUT_PREFIX}"
    done

    # deduplicate cumulative exclusion list (a sample could in principle
    # appear at multiple pass boundaries; FID/IID pairs are kept unique)
    if [ -s "${CUMULATIVE_EXCL}" ]; then
        sort -u "${CUMULATIVE_EXCL}" -o "${CUMULATIVE_EXCL}"
    fi

    SAMPLES_FINAL="${SAMPLES_PER_PASS[$((n-1))]}"

    echo "------------------------------------------------------------------"
    echo "[CRISP] ${mode} complete."
    echo ""
    echo "[CRISP] Samples before : ${SAMPLES_BEFORE}"
    echo "[CRISP] Total removed  : ${TOTAL_REMOVED}"
    echo "[CRISP] Samples after  : ${SAMPLES_FINAL}"
    echo "[CRISP] Exclusion list : ${CUMULATIVE_EXCL}"
    echo ""

    echo "[CRISP] Generating ${mode} call rate plots..."
    generate_plots "${mode}" "${STEP3_DIR}" \
        "${LOG_DIR}/step3_plots.log" "${IMISS_FILES[@]}"
    echo "[OK]    Plots saved to: ${STEP3_DIR}"
    echo ""

    # build pass summary string for report
    PASS_SUMMARY=""
    for i in "${!thresholds[@]}"; do
        PASS_SUMMARY="${PASS_SUMMARY}  Pass $((i+1)) (mind=${thresholds[$i]}) : ${REMOVED_PER_PASS[$i]} samples removed\n"
    done

    write_report "${mode}" "${MIND}" \
        "${SAMPLES_BEFORE}" "${SAMPLES_FINAL}" "${TOTAL_REMOVED}" \
        "${REMOVED_PER_PASS[0]:-}" "${REMOVED_PER_PASS[1]:-}" \
        "${REMOVED_PER_PASS[2]:-}" "${REMOVED_PER_PASS[3]:-}" \
        "${CURRENT_INPUT}"
}

##########################################################################
### CASCADE MODE
##########################################################################

run_cascade() {

    echo "##########################################################################"
    echo "###           STEP 3: SAMPLE CALL RATE (CASCADE MODE)                 ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Running cascade: 0.25 >> 0.20 >> 0.10 >> 0.05"
    echo ""

    run_multipass "CASCADE" "0.25" "0.20" "0.10" "0.05"
}

##########################################################################
### CUSTOM MODE
##########################################################################

run_custom() {

    echo "##########################################################################"
    echo "###           STEP 3: SAMPLE CALL RATE (CUSTOM MODE)                  ###"
    echo "##########################################################################"
    echo ""

    if [ -z "${SAMPLE_CUSTOM_TIERS}" ]; then
        echo "[CRISP] ERROR: CALLRATE_MODE = CUSTOM but SAMPLE_CUSTOM_TIERS is not set."
        echo "[CRISP]        Add SAMPLE_CUSTOM_TIERS to ${INSTRUCTION_FILE}"
        echo "[CRISP]        Example: SAMPLE_CUSTOM_TIERS = 0.30,0.20,0.10,0.05"
        exit 1
    fi

    validate_tiers "${SAMPLE_CUSTOM_TIERS}"

    IFS=',' read -ra CUSTOM_THRESHOLDS <<< "${SAMPLE_CUSTOM_TIERS}"

    echo "[CRISP] Running custom cascade: $(echo "${SAMPLE_CUSTOM_TIERS}" | tr ',' ' >> ')"
    echo ""

    run_multipass "CUSTOM" "${CUSTOM_THRESHOLDS[@]}"
}

##########################################################################
### WRITE REPORT
##########################################################################

write_report() {
    local mode="$1"
    local mind="$2"
    local before="$3"
    local after="$4"
    local total_removed="$5"
    local r1="$6" r2="$7" r3="$8" r4="$9"
    local final_prefix="${10}"

    {
        echo "=================================================================="
        echo "  CRISP: STEP 3 SAMPLE CALL RATE REPORT"
        echo "  Comprehensive Robust Integrated SNP Processing"
        echo "=================================================================="
        echo "  Project      : ${PROJECT_NAME}"
        echo "  Date         : $(date)"
        echo "  Mode         : ${mode}"
        echo "  MIND         : ${mind}"
        if [ "${mode}" = "CUSTOM" ]; then
        echo "  Custom tiers : ${SAMPLE_CUSTOM_TIERS}"
        fi
        echo "------------------------------------------------------------------"
        echo "  SAMPLE COUNTS"
        echo "  Samples before filtering : ${before}"
        echo "  Samples after filtering  : ${after}"
        echo "  Samples removed          : ${total_removed}"
        echo "  Variants (unchanged)     : ${VARIANTS}"
        if [ "${mode}" = "CASCADE" ]; then
        echo "------------------------------------------------------------------"
        echo "  CASCADE PASS SUMMARY"
        echo "  Pass 1 (mind=0.25) : ${r1} samples removed"
        echo "  Pass 2 (mind=0.20) : ${r2} samples removed"
        echo "  Pass 3 (mind=0.10) : ${r3} samples removed"
        echo "  Pass 4 (mind=0.05) : ${r4} samples removed"
        fi
        echo "------------------------------------------------------------------"
        echo "  OUTPUT"
        echo "  Final BED prefix  : ${final_prefix}"
        echo "  Exclusion list    : ${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step3.txt"
        echo "  Plots             : ${STEP3_DIR}"
        echo "  Log               : ${LOG_DIR}"
        echo "  Plot engine       : ${PLOT_ENGINE}"
        echo "  Plot colour mode  : ${CRISP_PAL_MODE}"
        echo "  Plot background   : ${CRISP_PAL_BACKGROUND}"
        echo "------------------------------------------------------------------"
        echo "  PLOTS GENERATED"
        if [ "${mode}" = "SIMPLE" ]; then
            echo "  step3_callrate_simple.pdf"
        else
            echo "  step3_callrate_${mode,,}_pass[N].pdf (one per pass)"
            echo "  step3_callrate_${mode,,}_faceted.pdf"
        fi
        echo "=================================================================="
        echo "  END OF REPORT"
        echo "=================================================================="
    } | tee "${REPORT_FILE}"

    echo ""
    echo "[CRISP] Report written to: ${REPORT_FILE}"
}

##########################################################################
### DISPATCH
##########################################################################

case "${CALLRATE_MODE^^}" in
    SIMPLE)
        run_simple
        ;;
    CASCADE)
        run_cascade
        ;;
    CUSTOM)
        run_custom
        ;;
    *)
        echo "[CRISP] ERROR: Unknown CALLRATE_MODE '${CALLRATE_MODE}'"
        echo "[CRISP]        Valid options: SIMPLE, CASCADE, CUSTOM"
        exit 1
        ;;
esac

##########################################################################
### DONE
##########################################################################

echo "[CRISP] Step 3 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
