#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Variant Call Rate (Step 4)
### Version: 0.2.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Filters variants based on genotype missingness rate using PLINK
# --geno flag. Also handles monomorphic SNP removal and MAF filtering.
#
# Three modes controlled by CALLRATE_MODE in crisp_instructions.txt:
#
#   SIMPLE  (default)
#     Single pass at GENO threshold (default 0.05)
#
#   CASCADE
#     Four fixed progressive passes:
#       Pass 1 : --geno 0.25
#       Pass 2 : --geno 0.20
#       Pass 3 : --geno 0.10
#       Pass 4 : --geno 0.05
#
#   CUSTOM
#     User-defined thresholds via VARIANT_CUSTOM_TIERS
#     Example: VARIANT_CUSTOM_TIERS = 0.30,0.20,0.15,0.05
#     Any number of tiers, strictly descending, between 0 and 1
#
# Sub-steps:
#   4a: Variant call rate filtering (SIMPLE, CASCADE or CUSTOM)
#   4b: Diagnostic counts (pre-filtering snapshot)
#   4c: Monomorphic SNP removal (REMOVE_MONO = YES/NO)
#   4d: MAF filtering (FILTER_MAF = YES/NO)
#
# Output:
#   Filtered BED/BIM/FAM per pass
#   .lmiss per-variant missingness statistics
#   PDF plots via plot_snprate.R or plot_snprate.py
#   Structured Step 4 report
#   Separate exclusion lists per sub-step
#
# Usage:
#   bash crisp_snprate.sh
#   bash crisp_snprate.sh --config my_project.txt
##########################################################################

set -euo pipefail

##########################################################################
### LOAD CRISP FLAVOUR (NON-NEGOTIABLE)
##########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime

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

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_snprate.sh --config <path>"
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
GENO=$(parse_param "GENO" "0.05")
VARIANT_CUSTOM_TIERS=$(parse_param "VARIANT_CUSTOM_TIERS" "")
REMOVE_MONO=$(parse_param "REMOVE_MONO" "YES")
FILTER_MAF=$(parse_param "FILTER_MAF" "YES")
MAF=$(parse_param "MAF" "0.01")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# resolve input, Step 3 output if available, else Step 2
STEP3_FINAL=$(parse_param "STEP3_FINAL_PREFIX" "")
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")

if [ -n "${STEP3_FINAL}" ] && [ -f "${STEP3_FINAL}.bed" ]; then
    INPUT_PREFIX="${STEP3_FINAL}"
    echo "[CRISP] Using Step 3 output as input."
elif [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    INPUT_PREFIX="${CONVERTED_PREFIX}"
    echo "[CRISP] Step 3 output not found, using Step 2 output."
else
    # fall back to naming convention
    INPUT_PREFIX="${OUTPUT_DIR}/step3_callrate/${PROJECT_NAME}_mind${GENO}"
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
    fi
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
echo "        CALLRATE_MODE        : ${CALLRATE_MODE}"
echo "        GENO                 : ${GENO}"
if [ "${CALLRATE_MODE^^}" = "CUSTOM" ]; then
echo "        VARIANT_CUSTOM_TIERS : ${VARIANT_CUSTOM_TIERS}"
fi
echo "        REMOVE_MONO          : ${REMOVE_MONO}"
echo "        FILTER_MAF           : ${FILTER_MAF}  (MAF = ${MAF})"
echo "        INPUT_PREFIX         : ${INPUT_PREFIX}"
echo "        PLOT_ENGINE          : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE INPUT
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

VARIANTS_BEFORE=$(wc -l < "${INPUT_PREFIX}.bim" | tr -d '[:space:]')
SAMPLES=$(wc -l < "${INPUT_PREFIX}.fam" | tr -d '[:space:]')

echo "[CRISP] Variants before filtering : ${VARIANTS_BEFORE}"
echo "[CRISP] Samples (unchanged)       : ${SAMPLES}"
echo ""

##########################################################################
### SET UP DIRECTORIES
##########################################################################

STEP4_DIR="${OUTPUT_DIR}/step4_snprate"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step4_snprate_report.txt"

mkdir -p "${STEP4_DIR}"
mkdir -p "${LOG_DIR}"

# per-run temp folder, avoids /tmp collisions on HPC shared nodes
CRISP_TMP="${OUTPUT_DIR}/.crisp_tmp_$$"
mkdir -p "${CRISP_TMP}"
trap "rm -rf ${CRISP_TMP}" EXIT

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
            echo "[CRISP]        Example: VARIANT_CUSTOM_TIERS = 0.30,0.20,0.10,0.05"
            exit 1
        fi
        prev="${tier}"
    done

    echo "[OK]    Custom tiers validated: ${tiers_str}"
}

##########################################################################
### GENERATE LMISS FILE
##########################################################################

generate_lmiss() {
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

##########################################################################
### STEP 4a: SIMPLE MODE
##########################################################################

run_simple() {

    echo "##########################################################################"
    echo "###          STEP 4a: VARIANT CALL RATE (SIMPLE MODE)                 ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Applying --geno ${GENO}..."
    echo ""

    OUT_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_geno${GENO}"
    LOG_FILE="${LOG_DIR}/step4_simple.log"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --geno "${GENO}" \
        --make-bed \
        --out "${OUT_PREFIX}" \
        >> "${LOG_FILE}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: PLINK --geno failed. Check log: ${LOG_FILE}"
        exit 1
    fi

    VARIANTS_AFTER=$(wc -l < "${OUT_PREFIX}.bim" | tr -d '[:space:]')
    VARIANTS_REMOVED=$((VARIANTS_BEFORE - VARIANTS_AFTER))

    echo "[OK]    SIMPLE filtering complete."
    echo ""
    echo "[CRISP] Variants before : ${VARIANTS_BEFORE}"
    echo "[CRISP] Variants removed: ${VARIANTS_REMOVED}"
    echo "[CRISP] Variants after  : ${VARIANTS_AFTER}"
    echo ""

    # lmiss for plotting
    LMISS_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_simple_lmiss"
    generate_lmiss "${INPUT_PREFIX}" "${LMISS_PREFIX}" "${LOG_FILE}"
    LMISS_FILE="${LMISS_PREFIX}.lmiss"

    echo "[CRISP] Generating variant call rate plots..."
    generate_plots "SIMPLE" "${STEP4_DIR}" \
        "${LOG_DIR}/step4_plots.log" "${LMISS_FILE}"
    echo "[OK]    Plots saved to: ${STEP4_DIR}"
    echo ""

    # exclusion list
    EXCLUSION_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
    awk '{print $2}' "${INPUT_PREFIX}.bim" | sort > ${CRISP_TMP}/crisp_vars_before.txt
    awk '{print $2}' "${OUT_PREFIX}.bim"   | sort > ${CRISP_TMP}/crisp_vars_after.txt
    comm -23 ${CRISP_TMP}/crisp_vars_before.txt \
             ${CRISP_TMP}/crisp_vars_after.txt > "${EXCLUSION_LIST}"

    echo "[CRISP] Exclusion list: ${EXCLUSION_LIST}"
    echo ""

    # set DIAG_PREFIX for downstream sub-steps
    DIAG_PREFIX="${OUT_PREFIX}"

    write_report "SIMPLE" "${GENO}" \
        "${VARIANTS_BEFORE}" "${VARIANTS_AFTER}" "${VARIANTS_REMOVED}" \
        "" "${OUT_PREFIX}"
}

##########################################################################
### STEP 4a: MULTI-PASS RUNNER (CASCADE and CUSTOM)
##########################################################################

run_multipass() {

    local mode="$1"
    shift
    local thresholds=("$@")
    local n="${#thresholds[@]}"

    CURRENT_INPUT="${INPUT_PREFIX}"
    LMISS_FILES=()
    TOTAL_REMOVED=0
    REMOVED_PER_PASS=()
    VARIANTS_PER_PASS=()
    CUMULATIVE_EXCL="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
    > "${CUMULATIVE_EXCL}"

    # track original variants for cumulative exclusion list
    awk '{print $2}' "${INPUT_PREFIX}.bim" | sort > ${CRISP_TMP}/crisp_vars_original.txt

    for i in "${!thresholds[@]}"; do

        PASS=$((i + 1))
        THRESHOLD="${thresholds[$i]}"
        OUT_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_geno${THRESHOLD}"
        LOG_FILE="${LOG_DIR}/step4_${mode,,}_pass${PASS}.log"

        echo "------------------------------------------------------------------"
        echo "[CRISP] Pass ${PASS} of ${n}, --geno ${THRESHOLD}"
        echo "------------------------------------------------------------------"
        echo ""

        VARIANTS_THIS_PASS=$(wc -l < "${CURRENT_INPUT}.bim" | tr -d '[:space:]')

        ${PLINK1} \
            --bfile "${CURRENT_INPUT}" \
            --geno "${THRESHOLD}" \
            --make-bed \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: PLINK --geno ${THRESHOLD} failed. Check log: ${LOG_FILE}"
            exit 1
        fi

        VARIANTS_AFTER_PASS=$(wc -l < "${OUT_PREFIX}.bim" | tr -d '[:space:]')
        REMOVED_THIS_PASS=$((VARIANTS_THIS_PASS - VARIANTS_AFTER_PASS))
        TOTAL_REMOVED=$((TOTAL_REMOVED + REMOVED_THIS_PASS))

        REMOVED_PER_PASS+=("${REMOVED_THIS_PASS}")
        VARIANTS_PER_PASS+=("${VARIANTS_AFTER_PASS}")

        echo "[OK]    Pass ${PASS} complete."
        echo "        Variants this pass : ${VARIANTS_THIS_PASS}"
        echo "        Variants removed   : ${REMOVED_THIS_PASS}"
        echo "        Variants remaining : ${VARIANTS_AFTER_PASS}"
        echo ""

        # lmiss for this pass
        LMISS_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_geno${THRESHOLD}_lmiss"
        generate_lmiss "${CURRENT_INPUT}" "${LMISS_PREFIX}" "${LOG_FILE}"
        LMISS_FILES+=("${LMISS_PREFIX}.lmiss")

        # update cumulative exclusion list
        awk '{print $2}' "${OUT_PREFIX}.bim" | sort \
            > ${CRISP_TMP}/crisp_vars_pass${PASS}.txt
        comm -23 ${CRISP_TMP}/crisp_vars_original.txt \
                 ${CRISP_TMP}/crisp_vars_pass${PASS}.txt > "${CUMULATIVE_EXCL}"

        # clean intermediate passes if not keeping
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && [ "${PASS}" -lt "${n}" ]; then
            for ext in .bed .bim .fam .log .nosex; do
                [ -f "${OUT_PREFIX}${ext}" ] && rm "${OUT_PREFIX}${ext}"
            done
        fi

        CURRENT_INPUT="${OUT_PREFIX}"
    done

    VARIANTS_FINAL="${VARIANTS_PER_PASS[$((n-1))]}"

    echo "------------------------------------------------------------------"
    echo "[CRISP] ${mode} complete."
    echo ""
    echo "[CRISP] Variants before : ${VARIANTS_BEFORE}"
    echo "[CRISP] Total removed   : ${TOTAL_REMOVED}"
    echo "[CRISP] Variants after  : ${VARIANTS_FINAL}"
    echo "[CRISP] Exclusion list  : ${CUMULATIVE_EXCL}"
    echo ""

    echo "[CRISP] Generating variant call rate plots..."
    generate_plots "${mode}" "${STEP4_DIR}" \
        "${LOG_DIR}/step4_plots.log" "${LMISS_FILES[@]}"
    echo "[OK]    Plots saved to: ${STEP4_DIR}"
    echo ""

    # set DIAG_PREFIX for downstream sub-steps
    DIAG_PREFIX="${CURRENT_INPUT}"

    write_report "${mode}" "${GENO}" \
        "${VARIANTS_BEFORE}" "${VARIANTS_FINAL}" "${TOTAL_REMOVED}" \
        "${REMOVED_PER_PASS[*]:-}" "${CURRENT_INPUT}"
}

##########################################################################
### STEP 4a: CASCADE MODE
##########################################################################

run_cascade() {

    echo "##########################################################################"
    echo "###          STEP 4a: VARIANT CALL RATE (CASCADE MODE)                ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Running cascade: 0.25 >> 0.20 >> 0.10 >> 0.05"
    echo ""

    run_multipass "CASCADE" "0.25" "0.20" "0.10" "0.05"
}

##########################################################################
### STEP 4a: CUSTOM MODE
##########################################################################

run_custom() {

    echo "##########################################################################"
    echo "###          STEP 4a: VARIANT CALL RATE (CUSTOM MODE)                 ###"
    echo "##########################################################################"
    echo ""

    if [ -z "${VARIANT_CUSTOM_TIERS}" ]; then
        echo "[CRISP] ERROR: CALLRATE_MODE = CUSTOM but VARIANT_CUSTOM_TIERS is not set."
        echo "[CRISP]        Add VARIANT_CUSTOM_TIERS to ${INSTRUCTION_FILE}"
        echo "[CRISP]        Example: VARIANT_CUSTOM_TIERS = 0.30,0.20,0.10,0.05"
        exit 1
    fi

    validate_tiers "${VARIANT_CUSTOM_TIERS}"

    IFS=',' read -ra CUSTOM_THRESHOLDS <<< "${VARIANT_CUSTOM_TIERS}"

    echo "[CRISP] Running custom cascade: $(echo "${VARIANT_CUSTOM_TIERS}" | tr ',' ' >> ')"
    echo ""

    run_multipass "CUSTOM" "${CUSTOM_THRESHOLDS[@]}"
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
### STEP 4b: DIAGNOSTIC COUNTS
### Pre-filtering snapshot, monomorphic and MAF counts before removal.
##########################################################################

echo "##########################################################################"
echo "###                   STEP 4b: DIAGNOSTIC COUNTS                      ###"
echo "##########################################################################"
echo ""

echo "[CRISP] Running frequency check for diagnostics..."

DIAG_LOG="${LOG_DIR}/step4_diagnostics.log"

${PLINK1} \
    --bfile "${DIAG_PREFIX}" \
    --freq \
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

echo ""
echo "[CRISP] Diagnostic counts (pre-removal):"
echo "        Total variants             : ${DIAG_TOTAL}"
echo "        Monomorphic (MAF = 0)      : ${DIAG_MONO}"
echo "        Below MAF ${MAF}              : ${DIAG_BELOW_MAF}"
echo "        Passing both filters       : ${DIAG_PASSING}"
echo ""

##########################################################################
### STEP 4c: MONOMORPHIC SNP REMOVAL
##########################################################################

MONO_REMOVED=0
MONO_EXCL_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4_monomorphic.txt"
VARIANTS_AFTER_MONO="${DIAG_TOTAL}"

if [ "${REMOVE_MONO^^}" = "YES" ]; then

    echo "##########################################################################"
    echo "###              STEP 4c: MONOMORPHIC SNP REMOVAL                     ###"
    echo "##########################################################################"
    echo ""

    MONO_IN="${DIAG_PREFIX}"
    MONO_OUT="${STEP4_DIR}/${PROJECT_NAME}_nomono"
    MONO_LOG="${LOG_DIR}/step4_mono.log"

    VARIANTS_PRE_MONO=$(wc -l < "${MONO_IN}.bim" | tr -d '[:space:]')

    echo "[CRISP] Removing monomorphic SNPs from: $(basename ${MONO_IN})"
    echo ""

    ${PLINK1} \
        --bfile "${MONO_IN}" \
        --mac 1 \
        --make-bed \
        --out "${MONO_OUT}" \
        >> "${MONO_LOG}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Monomorphic removal failed. Check log: ${MONO_LOG}"
        exit 1
    fi

    VARIANTS_AFTER_MONO=$(wc -l < "${MONO_OUT}.bim" | tr -d '[:space:]')
    MONO_REMOVED=$((VARIANTS_PRE_MONO - VARIANTS_AFTER_MONO))

    echo "[OK]    Monomorphic removal complete."
    echo ""
    echo "[CRISP] Before : ${VARIANTS_PRE_MONO}  Removed : ${MONO_REMOVED}  After : ${VARIANTS_AFTER_MONO}"
    echo ""

    # build exclusion list from diff
    awk '{print $2}' "${MONO_IN}.bim"  | sort > ${CRISP_TMP}/crisp_pre_mono.txt
    awk '{print $2}' "${MONO_OUT}.bim" | sort > ${CRISP_TMP}/crisp_post_mono.txt
    comm -23 ${CRISP_TMP}/crisp_pre_mono.txt \
             ${CRISP_TMP}/crisp_post_mono.txt > "${MONO_EXCL_LIST}"

    echo "[CRISP] Monomorphic exclusion list: ${MONO_EXCL_LIST}"
    echo ""

    if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && \
       [ "${MONO_IN}" != "${INPUT_PREFIX}" ]; then
        for ext in .bed .bim .fam .log .nosex; do
            [ -f "${MONO_IN}${ext}" ] && rm "${MONO_IN}${ext}"
        done
    fi

    CURRENT_WORKING="${MONO_OUT}"

else
    echo "[CRISP] REMOVE_MONO = NO, skipping."
    echo ""
    VARIANTS_AFTER_MONO="N/A (skipped)"
    CURRENT_WORKING="${DIAG_PREFIX}"
fi

##########################################################################
### STEP 4d: MAF FILTERING
##########################################################################

MAF_REMOVED=0
MAF_EXCL_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4_maf.txt"
VARIANTS_AFTER_MAF="N/A"

if [ "${FILTER_MAF^^}" = "YES" ]; then

    echo "##########################################################################"
    echo "###                    STEP 4d: MAF FILTERING                         ###"
    echo "##########################################################################"
    echo ""

    MAF_IN="${CURRENT_WORKING}"
    MAF_OUT="${STEP4_DIR}/${PROJECT_NAME}_maf${MAF}"
    MAF_LOG="${LOG_DIR}/step4_maf.log"

    VARIANTS_PRE_MAF=$(wc -l < "${MAF_IN}.bim" | tr -d '[:space:]')

    echo "[CRISP] Applying --maf ${MAF}..."
    echo ""

    ${PLINK1} \
        --bfile "${MAF_IN}" \
        --maf "${MAF}" \
        --make-bed \
        --out "${MAF_OUT}" \
        >> "${MAF_LOG}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: MAF filtering failed. Check log: ${MAF_LOG}"
        exit 1
    fi

    VARIANTS_AFTER_MAF=$(wc -l < "${MAF_OUT}.bim" | tr -d '[:space:]')
    MAF_REMOVED=$((VARIANTS_PRE_MAF - VARIANTS_AFTER_MAF))

    echo "[OK]    MAF filtering complete."
    echo ""
    echo "[CRISP] Before : ${VARIANTS_PRE_MAF}  Removed : ${MAF_REMOVED}  After : ${VARIANTS_AFTER_MAF}"
    echo ""

    # build exclusion list
    awk '{print $2}' "${MAF_IN}.bim"  | sort > ${CRISP_TMP}/crisp_pre_maf.txt
    awk '{print $2}' "${MAF_OUT}.bim" | sort > ${CRISP_TMP}/crisp_post_maf.txt
    comm -23 ${CRISP_TMP}/crisp_pre_maf.txt \
             ${CRISP_TMP}/crisp_post_maf.txt > "${MAF_EXCL_LIST}"

    echo "[CRISP] MAF exclusion list: ${MAF_EXCL_LIST}"
    echo ""

    if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && \
       [ "${MAF_IN}" != "${INPUT_PREFIX}" ]; then
        for ext in .bed .bim .fam .log .nosex; do
            [ -f "${MAF_IN}${ext}" ] && rm "${MAF_IN}${ext}"
        done
    fi

    FINAL_PREFIX="${MAF_OUT}"

else
    echo "[CRISP] FILTER_MAF = NO, skipping."
    echo ""
    VARIANTS_AFTER_MAF="N/A (skipped)"
    FINAL_PREFIX="${CURRENT_WORKING}"
fi

##########################################################################
### WRITE REPORT
##########################################################################

write_report() {
    local mode="$1"
    local geno="$2"
    local before="$3"
    local after="$4"
    local total_removed="$5"
    local pass_summary="$6"
    local final_prefix="$7"

    {
        echo "=================================================================="
        echo "  CRISP: STEP 4 VARIANT CALL RATE REPORT"
        echo "  Comprehensive Robust Integrated SNP Processing"
        echo "=================================================================="
        echo "  Project      : ${PROJECT_NAME}"
        echo "  Date         : $(date)"
        echo "  Mode         : ${mode}"
        echo "  GENO         : ${geno}"
        if [ "${mode}" = "CUSTOM" ]; then
        echo "  Custom tiers : ${VARIANT_CUSTOM_TIERS}"
        fi
        echo "------------------------------------------------------------------"
        echo "  DIAGNOSTIC COUNTS (pre-removal)"
        echo "  Total variants             : ${DIAG_TOTAL}"
        echo "  Monomorphic (MAF = 0)      : ${DIAG_MONO}"
        echo "  Below MAF ${MAF}              : ${DIAG_BELOW_MAF}"
        echo "  Passing both filters       : ${DIAG_PASSING}"
        echo "------------------------------------------------------------------"
        echo "  VARIANT CALL RATE FILTERING"
        echo "  Variants before  : ${before}"
        echo "  Variants after   : ${after}"
        echo "  Variants removed : ${total_removed}"
        echo "  Samples (unchanged) : ${SAMPLES}"
        if [ -n "${pass_summary}" ]; then
        echo "------------------------------------------------------------------"
        echo "  PASS SUMMARY"
        echo -e "${pass_summary}"
        fi
        echo "------------------------------------------------------------------"
        echo "  MONOMORPHIC SNP REMOVAL"
        echo "  REMOVE_MONO              : ${REMOVE_MONO}"
        if [ "${REMOVE_MONO^^}" = "YES" ]; then
        echo "  Monomorphic SNPs removed : ${MONO_REMOVED}"
        echo "  Variants after removal   : ${VARIANTS_AFTER_MONO}"
        echo "  Exclusion list           : ${MONO_EXCL_LIST}"
        fi
        echo "------------------------------------------------------------------"
        echo "  MAF FILTERING"
        echo "  FILTER_MAF               : ${FILTER_MAF}"
        echo "  MAF threshold            : ${MAF}"
        if [ "${FILTER_MAF^^}" = "YES" ]; then
        echo "  Variants removed (MAF)   : ${MAF_REMOVED}"
        echo "  Variants after MAF       : ${VARIANTS_AFTER_MAF}"
        echo "  Exclusion list           : ${MAF_EXCL_LIST}"
        fi
        echo "------------------------------------------------------------------"
        echo "  OUTPUT"
        echo "  Final BED prefix  : ${final_prefix}"
        echo "  Call rate excl.   : ${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
        if [ "${REMOVE_MONO^^}" = "YES" ]; then
        echo "  Monomorphic excl. : ${MONO_EXCL_LIST}"
        fi
        if [ "${FILTER_MAF^^}" = "YES" ]; then
        echo "  MAF excl.         : ${MAF_EXCL_LIST}"
        fi
        echo "  Plots             : ${STEP4_DIR}"
        echo "  Log               : ${LOG_DIR}"
        echo "=================================================================="
        echo "  END OF REPORT"
        echo "=================================================================="
    } | tee "${REPORT_FILE}"

    echo ""
    echo "[CRISP] Report written to: ${REPORT_FILE}"
}

##########################################################################
### DONE
##########################################################################

echo "[CRISP] Step 4 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
