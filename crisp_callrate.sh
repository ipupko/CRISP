#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Sample Call Rate (Step 3)
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
##########################################################################
### DESCRIPTION
### Filters samples based on genotype missingness rate using PLINK
### --mind flag.
###
### Two modes controlled by CALLRATE_MODE in crisp_instructions.txt:
###
###   SIMPLE  (default)
###     Single pass at MIND threshold (default 0.05)
###
###   CASCADE
###     Four progressive passes:
###       Pass 1 : --mind 0.25
###       Pass 2 : --mind 0.20
###       Pass 3 : --mind 0.10
###       Pass 4 : --mind 0.05
###
### Output:
###   .irem  -- excluded sample IDs per pass
###   .imiss -- per-sample missingness statistics
###   PDF plots via plot_callrate.R
###   Structured Step 3 report
###   Cumulative exclusion list for amendments step
##########################################################################

set -euo pipefail

##########################################################################
### LOAD CRISP FLAVOUR -- NON-NEGOTIABLE
##########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

INSTRUCTION_FILE="crisp_instructions.txt"

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
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

INPUT_FORMAT=$(parse_param "INPUT_FORMAT")
OUTPUT_DIR=$(parse_param "OUTPUT_DIR" "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME" "project")
KEEP_INTERMEDIATE=$(parse_param "KEEP_INTERMEDIATE" "YES")
CALLRATE_MODE=$(parse_param "CALLRATE_MODE" "SIMPLE")
MIND=$(parse_param "MIND" "0.05")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# Resolve input -- use converted prefix from Step 2 if available
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")
if [ -z "${CONVERTED_PREFIX}" ]; then
    # Fall back to Step 2 output naming convention
    CONVERTED_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
fi

##########################################################################
### HEADER
##########################################################################

echo ""
echo "##########################################################################"
echo "###                     CRISP - Version 0.1                           ###"
echo "###        Comprehensive Robust Integrated SNP Processing            ###"
echo "###                https://github.com/ipupko/CRISP                  ###"
echo "##########################################################################"
echo ""
echo "[CRISP] ${PREP_MSG}"
echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""
echo "[CRISP] Parameters read from instruction file:"
echo "        CALLRATE_MODE     : ${CALLRATE_MODE}"
echo "        MIND              : ${MIND}"
echo "        INPUT_PREFIX      : ${CONVERTED_PREFIX}"
echo "        OUTPUT_DIR        : ${OUTPUT_DIR}"
echo "        PROJECT_NAME      : ${PROJECT_NAME}"
echo "        KEEP_INTERMEDIATE : ${KEEP_INTERMEDIATE}"
echo "        PLOT_ENGINE       : ${PLOT_ENGINE}"
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

# Count samples before filtering
SAMPLES_BEFORE=$(wc -l < "${CONVERTED_PREFIX}.fam" | tr -d '[:space:]')
VARIANTS=$(wc -l < "${CONVERTED_PREFIX}.bim" | tr -d '[:space:]')

echo "[CRISP] Samples before filtering : ${SAMPLES_BEFORE}"
echo "[CRISP] Variants                 : ${VARIANTS}"
echo ""

##########################################################################
### SIMPLE MODE
##########################################################################

run_simple() {

    echo "##########################################################################"
    echo "###            STEP 3: SAMPLE CALL RATE (SIMPLE MODE)                 ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Applying single --mind ${MIND} filter..."
    echo ""

    OUT_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_mind${MIND}"
    LOG_FILE="${LOG_DIR}/step3_simple.log"

    ${PLINK1} \
        --bfile "${CONVERTED_PREFIX}" \
        --mind "${MIND}" \
        --make-bed \
        --out "${OUT_PREFIX}" \
        >> "${LOG_FILE}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: PLINK --mind failed. Check log: ${LOG_FILE}"
        exit 1
    fi

    # Count results
    SAMPLES_AFTER=$(wc -l < "${OUT_PREFIX}.fam" | tr -d '[:space:]')
    SAMPLES_REMOVED=$((SAMPLES_BEFORE - SAMPLES_AFTER))

    # Check if .irem exists (PLINK only creates it if samples were removed)
    IREM_FILE="${OUT_PREFIX}.irem"
    if [ ! -f "${IREM_FILE}" ]; then
        touch "${IREM_FILE}"
        echo "[CRISP] No samples removed at --mind ${MIND}."
    fi

    echo "[OK]    SIMPLE filtering complete."
    echo ""
    echo "[CRISP] Samples before : ${SAMPLES_BEFORE}"
    echo "[CRISP] Samples removed: ${SAMPLES_REMOVED}"
    echo "[CRISP] Samples after  : ${SAMPLES_AFTER}"
    echo ""

    # Generate missingness stats for plotting
    IMISS_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_simple_imiss"
    ${PLINK1} \
        --bfile "${CONVERTED_PREFIX}" \
        --missing \
        --out "${IMISS_PREFIX}" \
        >> "${LOG_FILE}" 2>&1

    IMISS_FILE="${IMISS_PREFIX}.imiss"

    # Generate plots
    echo "[CRISP] Generating call rate plots..."
    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_callrate.R" \
            "SIMPLE" "${MIND}" "${STEP3_DIR}" "${IMISS_FILE}" \
            >> "${LOG_DIR}/step3_plots.log" 2>&1
        echo "[OK]    Plots saved to: ${STEP3_DIR}"
    fi
    echo ""

    # Write cumulative exclusion list
    EXCLUSION_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step3.txt"
    if [ -f "${IREM_FILE}" ] && [ -s "${IREM_FILE}" ]; then
        cp "${IREM_FILE}" "${EXCLUSION_LIST}"
    else
        touch "${EXCLUSION_LIST}"
    fi

    echo "[CRISP] Exclusion list written to: ${EXCLUSION_LIST}"
    echo ""

    # Write report
    write_report "SIMPLE" "${MIND}" \
        "${SAMPLES_BEFORE}" "${SAMPLES_AFTER}" "${SAMPLES_REMOVED}" \
        "" "" "" "" \
        "${OUT_PREFIX}"
}

##########################################################################
### CASCADE MODE
##########################################################################

run_cascade() {

    echo "##########################################################################"
    echo "###           STEP 3: SAMPLE CALL RATE (CASCADE MODE)                 ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Running cascade: --mind 0.25 >> 0.20 >> 0.10 >> 0.05"
    echo ""

    CASCADE_THRESHOLDS=("0.25" "0.20" "0.10" "0.05")
    CURRENT_INPUT="${CONVERTED_PREFIX}"
    IMISS_FILES=()

    TOTAL_REMOVED=0
    REMOVED_PER_PASS=()
    SAMPLES_PER_PASS=()

    CUMULATIVE_IREM="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step3.txt"
    > "${CUMULATIVE_IREM}"

    for i in "${!CASCADE_THRESHOLDS[@]}"; do

        PASS=$((i + 1))
        THRESHOLD="${CASCADE_THRESHOLDS[$i]}"
        OUT_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_cascade_pass${PASS}_mind${THRESHOLD}"
        LOG_FILE="${LOG_DIR}/step3_cascade_pass${PASS}.log"

        echo "------------------------------------------------------------------"
        echo "[CRISP] Pass ${PASS} of 4 -- --mind ${THRESHOLD}"
        echo "------------------------------------------------------------------"
        echo ""

        SAMPLES_THIS_PASS=$(wc -l < "${CURRENT_INPUT}.fam" | tr -d '[:space:]')

        ${PLINK1} \
            --bfile "${CURRENT_INPUT}" \
            --mind "${THRESHOLD}" \
            --make-bed \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: PLINK --mind ${THRESHOLD} failed. Check log: ${LOG_FILE}"
            exit 1
        fi

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

        # Append to cumulative exclusion list
        IREM_FILE="${OUT_PREFIX}.irem"
        if [ -f "${IREM_FILE}" ] && [ -s "${IREM_FILE}" ]; then
            cat "${IREM_FILE}" >> "${CUMULATIVE_IREM}"
        fi

        # Generate .imiss for this pass input
        IMISS_PREFIX="${STEP3_DIR}/${PROJECT_NAME}_cascade_pass${PASS}_imiss"
        ${PLINK1} \
            --bfile "${CURRENT_INPUT}" \
            --missing \
            --out "${IMISS_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        IMISS_FILES+=("${IMISS_PREFIX}.imiss")

        # Clean up intermediate if requested (keep final pass always)
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && [ "${PASS}" -lt 4 ]; then
            for ext in .bed .bim .fam .log .nosex; do
                [ -f "${OUT_PREFIX}${ext}" ] && rm "${OUT_PREFIX}${ext}"
            done
            echo "[CRISP] Intermediate Pass ${PASS} files removed."
        fi

        CURRENT_INPUT="${OUT_PREFIX}"
    done

    SAMPLES_FINAL="${SAMPLES_PER_PASS[3]}"

    echo "------------------------------------------------------------------"
    echo "[CRISP] CASCADE complete."
    echo ""
    echo "[CRISP] Samples before cascade : ${SAMPLES_BEFORE}"
    echo "[CRISP] Total samples removed  : ${TOTAL_REMOVED}"
    echo "[CRISP] Samples after cascade  : ${SAMPLES_FINAL}"
    echo "[CRISP] Exclusion list         : ${CUMULATIVE_IREM}"
    echo ""

    # Generate plots
    echo "[CRISP] Generating cascade call rate plots..."
    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_callrate.R" \
            "CASCADE" "${MIND}" "${STEP3_DIR}" \
            "${IMISS_FILES[0]}" "${IMISS_FILES[1]}" \
            "${IMISS_FILES[2]}" "${IMISS_FILES[3]}" \
            >> "${LOG_DIR}/step3_plots.log" 2>&1
        echo "[OK]    Individual pass plots saved to  : ${STEP3_DIR}"
        echo "[OK]    Faceted comparison plot saved to: ${STEP3_DIR}"
    fi
    echo ""

    write_report "CASCADE" "${MIND}" \
        "${SAMPLES_BEFORE}" "${SAMPLES_FINAL}" "${TOTAL_REMOVED}" \
        "${REMOVED_PER_PASS[0]}" "${REMOVED_PER_PASS[1]}" \
        "${REMOVED_PER_PASS[2]}" "${REMOVED_PER_PASS[3]}" \
        "${CURRENT_INPUT}"
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
        echo "  CRISP -- STEP 3 SAMPLE CALL RATE REPORT"
        echo "  Comprehensive Robust Integrated SNP Processing"
        echo "=================================================================="
        echo "  Project      : ${PROJECT_NAME}"
        echo "  Date         : $(date)"
        echo "  Mode         : ${mode}"
        echo "  MIND         : ${mind}"
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
        echo "------------------------------------------------------------------"
        echo "  PLOTS GENERATED"
        if [ "${mode}" = "SIMPLE" ]; then
            echo "  step3_callrate_simple.pdf"
        else
            echo "  step3_callrate_cascade_pass1.pdf"
            echo "  step3_callrate_cascade_pass2.pdf"
            echo "  step3_callrate_cascade_pass3.pdf"
            echo "  step3_callrate_cascade_pass4.pdf"
            echo "  step3_callrate_cascade_faceted.pdf"
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
    *)
        echo "[CRISP] ERROR: Unknown CALLRATE_MODE '${CALLRATE_MODE}'"
        echo "[CRISP]        Valid options: SIMPLE, CASCADE"
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
