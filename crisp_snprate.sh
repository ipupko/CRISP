#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Variant Call Rate (Step 4)
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
##########################################################################
### DESCRIPTION
### Filters variants based on genotype missingness rate using PLINK
### --geno flag.
###
### Three modes controlled by CALLRATE_MODE in crisp_instructions.txt:
###
###   SIMPLE  (default)
###     Single pass at GENO threshold (default 0.05)
###
###   CASCADE
###     Four fixed progressive passes:
###       Pass 1 : --geno 0.25
###       Pass 2 : --geno 0.20
###       Pass 3 : --geno 0.10
###       Pass 4 : --geno 0.05
###
###   CUSTOM
###     User-defined thresholds via VARIANT_CUSTOM_TIERS
###     Example: VARIANT_CUSTOM_TIERS = 0.30,0.20,0.15,0.05
###     Any number of tiers, must be descending, between 0 and 1
###
### Output:
###   Filtered BED/BIM/FAM per pass
###   .lmiss -- per-variant missingness statistics
###   PDF plots via plot_snprate.R
###   Structured Step 4 report
###   Variant exclusion list for amendments step
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

OUTPUT_DIR=$(parse_param "OUTPUT_DIR" "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME" "project")
KEEP_INTERMEDIATE=$(parse_param "KEEP_INTERMEDIATE" "YES")
CALLRATE_MODE=$(parse_param "CALLRATE_MODE" "SIMPLE")
GENO=$(parse_param "GENO" "0.05")
VARIANT_CUSTOM_TIERS=$(parse_param "VARIANT_CUSTOM_TIERS" "")
REMOVE_MONO=$(parse_param "REMOVE_MONO" "YES")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# Resolve input -- use Step 3 output if available, else Step 2
STEP3_FINAL=$(parse_param "STEP3_FINAL_PREFIX" "")
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")

if [ -n "${STEP3_FINAL}" ] && [ -f "${STEP3_FINAL}.bed" ]; then
    INPUT_PREFIX="${STEP3_FINAL}"
    echo "[CRISP] Using Step 3 output as input."
elif [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    INPUT_PREFIX="${CONVERTED_PREFIX}"
    echo "[CRISP] Step 3 output not found. Using Step 2 converted output."
else
    # Fall back to convention
    INPUT_PREFIX="${OUTPUT_DIR}/step3_callrate/${PROJECT_NAME}_mind${GENO}"
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
    fi
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
echo "        CALLRATE_MODE        : ${CALLRATE_MODE}"
echo "        GENO                 : ${GENO}"
echo "        VARIANT_CUSTOM_TIERS : ${VARIANT_CUSTOM_TIERS:-N/A}"
echo "        REMOVE_MONO          : ${REMOVE_MONO}"
echo "        INPUT_PREFIX         : ${INPUT_PREFIX}"
echo "        OUTPUT_DIR           : ${OUTPUT_DIR}"
echo "        PROJECT_NAME         : ${PROJECT_NAME}"
echo "        KEEP_INTERMEDIATE    : ${KEEP_INTERMEDIATE}"
echo "        PLOT_ENGINE          : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE INPUT FILES
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

##########################################################################
### VALIDATE CUSTOM TIERS
##########################################################################

validate_tiers() {
    local tiers_str="$1"
    local prev=1.0
    local valid=1

    IFS=',' read -ra tiers <<< "${tiers_str}"

    if [ "${#tiers[@]}" -lt 2 ]; then
        echo "[CRISP] ERROR: CUSTOM mode requires at least 2 tiers."
        echo "[CRISP]        Got: ${tiers_str}"
        exit 1
    fi

    for tier in "${tiers[@]}"; do
        # Check numeric
        if ! [[ "${tier}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
            echo "[CRISP] ERROR: Invalid tier value '${tier}' -- must be numeric between 0 and 1."
            exit 1
        fi
        # Check between 0 and 1
        if (( $(echo "${tier} <= 0 || ${tier} >= 1" | bc -l) )); then
            echo "[CRISP] ERROR: Tier '${tier}' out of range -- must be between 0 and 1 (exclusive)."
            exit 1
        fi
        # Check descending
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
### SIMPLE MODE
##########################################################################

run_simple() {

    echo "##########################################################################"
    echo "###            STEP 4: VARIANT CALL RATE (SIMPLE MODE)                ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Applying single --geno ${GENO} filter..."
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

    # Generate lmiss for plotting
    LMISS_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_simple_lmiss"
    generate_lmiss "${INPUT_PREFIX}" "${LMISS_PREFIX}" "${LOG_FILE}"
    LMISS_FILE="${LMISS_PREFIX}.lmiss"

    # Generate plots
    echo "[CRISP] Generating variant call rate plots..."
    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_snprate.R" \
            "SIMPLE" "${GENO}" "${STEP4_DIR}" "${LMISS_FILE}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Plot saved to: ${STEP4_DIR}/step4_snprate_simple.pdf"
    elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
        python3 "${SCRIPT_DIR}/scripts/plot_snprate.py" \
            "SIMPLE" "${GENO}" "${STEP4_DIR}" "${LMISS_FILE}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Plot saved to: ${STEP4_DIR}/step4_snprate_simple.pdf"
    fi
    echo ""

    # Write exclusion list -- variants in input not in output
    EXCLUSION_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
    awk '{print $2}' "${INPUT_PREFIX}.bim" > /tmp/crisp_vars_before.txt
    awk '{print $2}' "${OUT_PREFIX}.bim"   > /tmp/crisp_vars_after.txt
    comm -23 <(sort /tmp/crisp_vars_before.txt) \
             <(sort /tmp/crisp_vars_after.txt) > "${EXCLUSION_LIST}"

    echo "[CRISP] Exclusion list written to: ${EXCLUSION_LIST}"
    echo ""

    write_report "SIMPLE" "${GENO}" \
        "${VARIANTS_BEFORE}" "${VARIANTS_AFTER}" "${VARIANTS_REMOVED}" \
        "" "" "" "" \
        "${OUT_PREFIX}"
}

##########################################################################
### CASCADE MODE
##########################################################################

run_cascade() {

    echo "##########################################################################"
    echo "###           STEP 4: VARIANT CALL RATE (CASCADE MODE)                ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Running cascade: --geno 0.25 >> 0.20 >> 0.10 >> 0.05"
    echo ""

    CASCADE_THRESHOLDS=("0.25" "0.20" "0.10" "0.05")
    run_multipass "CASCADE" "${CASCADE_THRESHOLDS[@]}"
}

##########################################################################
### CUSTOM MODE
##########################################################################

run_custom() {

    echo "##########################################################################"
    echo "###           STEP 4: VARIANT CALL RATE (CUSTOM MODE)                 ###"
    echo "##########################################################################"
    echo ""

    if [ -z "${VARIANT_CUSTOM_TIERS}" ]; then
        echo "[CRISP] ERROR: CALLRATE_MODE = CUSTOM but VARIANT_CUSTOM_TIERS is not set."
        echo "[CRISP]        Please add VARIANT_CUSTOM_TIERS to ${INSTRUCTION_FILE}"
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
### MULTI-PASS RUNNER (shared by CASCADE and CUSTOM)
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

    # Track all variants from input for exclusion list
    awk '{print $2}' "${INPUT_PREFIX}.bim" | sort > /tmp/crisp_vars_original.txt

    for i in "${!thresholds[@]}"; do

        PASS=$((i + 1))
        THRESHOLD="${thresholds[$i]}"
        OUT_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_geno${THRESHOLD}"
        LOG_FILE="${LOG_DIR}/step4_${mode,,}_pass${PASS}.log"

        echo "------------------------------------------------------------------"
        echo "[CRISP] Pass ${PASS} of ${n} -- --geno ${THRESHOLD}"
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

        # Generate lmiss for this pass input (pre-filter distribution)
        LMISS_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_geno${THRESHOLD}_lmiss"
        generate_lmiss "${CURRENT_INPUT}" "${LMISS_PREFIX}" "${LOG_FILE}"
        LMISS_FILES+=("${LMISS_PREFIX}.lmiss")

        # Update cumulative exclusion list
        awk '{print $2}' "${OUT_PREFIX}.bim" | sort > /tmp/crisp_vars_pass${PASS}.txt
        comm -23 /tmp/crisp_vars_original.txt \
                 /tmp/crisp_vars_pass${PASS}.txt > "${CUMULATIVE_EXCL}"

        # Clean intermediate if requested (keep final pass)
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && [ "${PASS}" -lt "${n}" ]; then
            for ext in .bed .bim .fam .log .nosex; do
                [ -f "${OUT_PREFIX}${ext}" ] && rm "${OUT_PREFIX}${ext}"
            done
            echo "[CRISP] Intermediate Pass ${PASS} files removed."
        fi

        CURRENT_INPUT="${OUT_PREFIX}"
    done

    VARIANTS_FINAL="${VARIANTS_PER_PASS[$((n-1))]}"

    echo "------------------------------------------------------------------"
    echo "[CRISP] ${mode} complete."
    echo ""
    echo "[CRISP] Variants before  : ${VARIANTS_BEFORE}"
    echo "[CRISP] Total removed    : ${TOTAL_REMOVED}"
    echo "[CRISP] Variants after   : ${VARIANTS_FINAL}"
    echo "[CRISP] Exclusion list   : ${CUMULATIVE_EXCL}"
    echo ""

    # Generate plots
    echo "[CRISP] Generating variant call rate plots..."
    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_snprate.R" \
            "${mode}" "${GENO}" "${STEP4_DIR}" \
            "${LMISS_FILES[@]}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Individual pass plots saved to  : ${STEP4_DIR}"
        echo "[OK]    Faceted comparison plot saved to: ${STEP4_DIR}"
    elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
        python3 "${SCRIPT_DIR}/scripts/plot_snprate.py" \
            "${mode}" "${GENO}" "${STEP4_DIR}" \
            "${LMISS_FILES[@]}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Individual pass plots saved to  : ${STEP4_DIR}"
        echo "[OK]    Faceted comparison plot saved to: ${STEP4_DIR}"
    fi
    echo ""

    # Build pass summary string for report
    PASS_SUMMARY=""
    for i in "${!thresholds[@]}"; do
        PASS_SUMMARY="${PASS_SUMMARY}  Pass $((i+1)) (geno=${thresholds[$i]}) : ${REMOVED_PER_PASS[$i]} variants removed\n"
    done

    write_report "${mode}" "${GENO}" \
        "${VARIANTS_BEFORE}" "${VARIANTS_FINAL}" "${TOTAL_REMOVED}" \
        "${PASS_SUMMARY}" \
        "${CURRENT_INPUT}"
}

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
        echo "  CRISP -- STEP 4 VARIANT CALL RATE REPORT"
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
        echo "  VARIANT COUNTS (CALL RATE FILTERING)"
        echo "  Variants before filtering : ${before}"
        echo "  Variants after filtering  : ${after}"
        echo "  Variants removed          : ${total_removed}"
        echo "  Samples (unchanged)       : ${SAMPLES}"
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
        echo "  OUTPUT"
        echo "  Final BED prefix  : ${final_prefix}"
        echo "  Call rate excl.   : ${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
        if [ "${REMOVE_MONO^^}" = "YES" ]; then
        echo "  Monomorphic excl. : ${MONO_EXCL_LIST}"
        fi
        echo "  Plots             : ${STEP4_DIR}"
        echo "  Log               : ${LOG_DIR}"
        echo "------------------------------------------------------------------"
        echo "  PLOTS GENERATED"
        if [ "${mode}" = "SIMPLE" ]; then
            echo "  step4_snprate_simple.pdf"
        else
            echo "  step4_snprate_${mode,,}_pass[N].pdf (one per pass)"
            echo "  step4_snprate_${mode,,}_faceted.pdf"
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
### MONOMORPHIC SNP REMOVAL (OPTIONAL)
### Applied after call rate filtering, before final report.
### Uses --mac 1 to remove variants with minor allele count < 1.
### Controlled by REMOVE_MONO in crisp_instructions.txt.
##########################################################################

MONO_REMOVED=0
MONO_EXCL_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4_monomorphic.txt"

if [ "${REMOVE_MONO^^}" = "YES" ]; then

    echo "##########################################################################"
    echo "###              STEP 4b: MONOMORPHIC SNP REMOVAL                     ###"
    echo "##########################################################################"
    echo ""

    # Resolve the final prefix from call rate step
    FINAL_CALLRATE_PREFIX=$(ls -t "${STEP4_DIR}"/${PROJECT_NAME}*.bed 2>/dev/null \
        | head -1 | sed 's/\.bed$//')

    if [ -z "${FINAL_CALLRATE_PREFIX}" ]; then
        echo "[CRISP] WARNING: Could not resolve call rate output. Skipping mono removal."
    else
        MONO_IN="${FINAL_CALLRATE_PREFIX}"
        MONO_OUT="${STEP4_DIR}/${PROJECT_NAME}_nomono"
        MONO_LOG="${LOG_DIR}/step4_mono.log"

        VARIANTS_PRE_MONO=$(wc -l < "${MONO_IN}.bim" | tr -d '[:space:]')

        echo "[CRISP] Removing monomorphic SNPs from: $(basename ${MONO_IN})"
        echo "[CRISP] Running: ${PLINK1} --bfile ${MONO_IN} --mac 1 --make-bed --out ${MONO_OUT}"
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
        echo "[CRISP] Variants before mono removal : ${VARIANTS_PRE_MONO}"
        echo "[CRISP] Monomorphic SNPs removed     : ${MONO_REMOVED}"
        echo "[CRISP] Variants after mono removal  : ${VARIANTS_AFTER_MONO}"
        echo ""

        # Separate monomorphic exclusion list
        awk '{print $2}' "${MONO_IN}.bim"  | sort > /tmp/crisp_pre_mono.txt
        awk '{print $2}' "${MONO_OUT}.bim" | sort > /tmp/crisp_post_mono.txt
        comm -23 /tmp/crisp_pre_mono.txt \
                 /tmp/crisp_post_mono.txt > "${MONO_EXCL_LIST}"

        echo "[CRISP] Monomorphic exclusion list   : ${MONO_EXCL_LIST}"
        echo ""

        # Clean up if requested
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && \
           [ "${MONO_IN}" != "${INPUT_PREFIX}" ]; then
            for ext in .bed .bim .fam .log .nosex; do
                [ -f "${MONO_IN}${ext}" ] && rm "${MONO_IN}${ext}"
            done
            echo "[CRISP] Intermediate pre-mono files removed."
            echo ""
        fi
    fi

else
    echo "[CRISP] REMOVE_MONO = NO -- monomorphic SNP removal skipped."
    echo ""
    VARIANTS_AFTER_MONO="N/A (step skipped)"
fi

##########################################################################
### DONE
##########################################################################

echo "[CRISP] Step 4 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""    local default="${2:-}"
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
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# Resolve input -- use Step 3 output if available, else Step 2
STEP3_FINAL=$(parse_param "STEP3_FINAL_PREFIX" "")
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")

if [ -n "${STEP3_FINAL}" ] && [ -f "${STEP3_FINAL}.bed" ]; then
    INPUT_PREFIX="${STEP3_FINAL}"
    echo "[CRISP] Using Step 3 output as input."
elif [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    INPUT_PREFIX="${CONVERTED_PREFIX}"
    echo "[CRISP] Step 3 output not found. Using Step 2 converted output."
else
    # Fall back to convention
    INPUT_PREFIX="${OUTPUT_DIR}/step3_callrate/${PROJECT_NAME}_mind${GENO}"
    if [ ! -f "${INPUT_PREFIX}.bed" ]; then
        INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
    fi
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
echo "        CALLRATE_MODE        : ${CALLRATE_MODE}"
echo "        GENO                 : ${GENO}"
echo "        VARIANT_CUSTOM_TIERS : ${VARIANT_CUSTOM_TIERS:-N/A}"
echo "        INPUT_PREFIX         : ${INPUT_PREFIX}"
echo "        OUTPUT_DIR           : ${OUTPUT_DIR}"
echo "        PROJECT_NAME         : ${PROJECT_NAME}"
echo "        KEEP_INTERMEDIATE    : ${KEEP_INTERMEDIATE}"
echo "        PLOT_ENGINE          : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE INPUT FILES
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

##########################################################################
### VALIDATE CUSTOM TIERS
##########################################################################

validate_tiers() {
    local tiers_str="$1"
    local prev=1.0
    local valid=1

    IFS=',' read -ra tiers <<< "${tiers_str}"

    if [ "${#tiers[@]}" -lt 2 ]; then
        echo "[CRISP] ERROR: CUSTOM mode requires at least 2 tiers."
        echo "[CRISP]        Got: ${tiers_str}"
        exit 1
    fi

    for tier in "${tiers[@]}"; do
        # Check numeric
        if ! [[ "${tier}" =~ ^[0-9]*\.?[0-9]+$ ]]; then
            echo "[CRISP] ERROR: Invalid tier value '${tier}' -- must be numeric between 0 and 1."
            exit 1
        fi
        # Check between 0 and 1
        if (( $(echo "${tier} <= 0 || ${tier} >= 1" | bc -l) )); then
            echo "[CRISP] ERROR: Tier '${tier}' out of range -- must be between 0 and 1 (exclusive)."
            exit 1
        fi
        # Check descending
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
### SIMPLE MODE
##########################################################################

run_simple() {

    echo "##########################################################################"
    echo "###            STEP 4: VARIANT CALL RATE (SIMPLE MODE)                ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Applying single --geno ${GENO} filter..."
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

    # Generate lmiss for plotting
    LMISS_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_simple_lmiss"
    generate_lmiss "${INPUT_PREFIX}" "${LMISS_PREFIX}" "${LOG_FILE}"
    LMISS_FILE="${LMISS_PREFIX}.lmiss"

    # Generate plots
    echo "[CRISP] Generating variant call rate plots..."
    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_snprate.R" \
            "SIMPLE" "${GENO}" "${STEP4_DIR}" "${LMISS_FILE}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Plot saved to: ${STEP4_DIR}/step4_snprate_simple.pdf"
    elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
        python3 "${SCRIPT_DIR}/scripts/plot_snprate.py" \
            "SIMPLE" "${GENO}" "${STEP4_DIR}" "${LMISS_FILE}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Plot saved to: ${STEP4_DIR}/step4_snprate_simple.pdf"
    fi
    echo ""

    # Write exclusion list -- variants in input not in output
    EXCLUSION_LIST="${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
    awk '{print $2}' "${INPUT_PREFIX}.bim" > /tmp/crisp_vars_before.txt
    awk '{print $2}' "${OUT_PREFIX}.bim"   > /tmp/crisp_vars_after.txt
    comm -23 <(sort /tmp/crisp_vars_before.txt) \
             <(sort /tmp/crisp_vars_after.txt) > "${EXCLUSION_LIST}"

    echo "[CRISP] Exclusion list written to: ${EXCLUSION_LIST}"
    echo ""

    write_report "SIMPLE" "${GENO}" \
        "${VARIANTS_BEFORE}" "${VARIANTS_AFTER}" "${VARIANTS_REMOVED}" \
        "" "" "" "" \
        "${OUT_PREFIX}"
}

##########################################################################
### CASCADE MODE
##########################################################################

run_cascade() {

    echo "##########################################################################"
    echo "###           STEP 4: VARIANT CALL RATE (CASCADE MODE)                ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Running cascade: --geno 0.25 >> 0.20 >> 0.10 >> 0.05"
    echo ""

    CASCADE_THRESHOLDS=("0.25" "0.20" "0.10" "0.05")
    run_multipass "CASCADE" "${CASCADE_THRESHOLDS[@]}"
}

##########################################################################
### CUSTOM MODE
##########################################################################

run_custom() {

    echo "##########################################################################"
    echo "###           STEP 4: VARIANT CALL RATE (CUSTOM MODE)                 ###"
    echo "##########################################################################"
    echo ""

    if [ -z "${VARIANT_CUSTOM_TIERS}" ]; then
        echo "[CRISP] ERROR: CALLRATE_MODE = CUSTOM but VARIANT_CUSTOM_TIERS is not set."
        echo "[CRISP]        Please add VARIANT_CUSTOM_TIERS to ${INSTRUCTION_FILE}"
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
### MULTI-PASS RUNNER (shared by CASCADE and CUSTOM)
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

    # Track all variants from input for exclusion list
    awk '{print $2}' "${INPUT_PREFIX}.bim" | sort > /tmp/crisp_vars_original.txt

    for i in "${!thresholds[@]}"; do

        PASS=$((i + 1))
        THRESHOLD="${thresholds[$i]}"
        OUT_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_geno${THRESHOLD}"
        LOG_FILE="${LOG_DIR}/step4_${mode,,}_pass${PASS}.log"

        echo "------------------------------------------------------------------"
        echo "[CRISP] Pass ${PASS} of ${n} -- --geno ${THRESHOLD}"
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

        # Generate lmiss for this pass input (pre-filter distribution)
        LMISS_PREFIX="${STEP4_DIR}/${PROJECT_NAME}_${mode,,}_pass${PASS}_geno${THRESHOLD}_lmiss"
        generate_lmiss "${CURRENT_INPUT}" "${LMISS_PREFIX}" "${LOG_FILE}"
        LMISS_FILES+=("${LMISS_PREFIX}.lmiss")

        # Update cumulative exclusion list
        awk '{print $2}' "${OUT_PREFIX}.bim" | sort > /tmp/crisp_vars_pass${PASS}.txt
        comm -23 /tmp/crisp_vars_original.txt \
                 /tmp/crisp_vars_pass${PASS}.txt > "${CUMULATIVE_EXCL}"

        # Clean intermediate if requested (keep final pass)
        if [ "${KEEP_INTERMEDIATE^^}" = "NO" ] && [ "${PASS}" -lt "${n}" ]; then
            for ext in .bed .bim .fam .log .nosex; do
                [ -f "${OUT_PREFIX}${ext}" ] && rm "${OUT_PREFIX}${ext}"
            done
            echo "[CRISP] Intermediate Pass ${PASS} files removed."
        fi

        CURRENT_INPUT="${OUT_PREFIX}"
    done

    VARIANTS_FINAL="${VARIANTS_PER_PASS[$((n-1))]}"

    echo "------------------------------------------------------------------"
    echo "[CRISP] ${mode} complete."
    echo ""
    echo "[CRISP] Variants before  : ${VARIANTS_BEFORE}"
    echo "[CRISP] Total removed    : ${TOTAL_REMOVED}"
    echo "[CRISP] Variants after   : ${VARIANTS_FINAL}"
    echo "[CRISP] Exclusion list   : ${CUMULATIVE_EXCL}"
    echo ""

    # Generate plots
    echo "[CRISP] Generating variant call rate plots..."
    if [ "${PLOT_ENGINE^^}" = "R" ]; then
        ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_snprate.R" \
            "${mode}" "${GENO}" "${STEP4_DIR}" \
            "${LMISS_FILES[@]}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Individual pass plots saved to  : ${STEP4_DIR}"
        echo "[OK]    Faceted comparison plot saved to: ${STEP4_DIR}"
    elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
        python3 "${SCRIPT_DIR}/scripts/plot_snprate.py" \
            "${mode}" "${GENO}" "${STEP4_DIR}" \
            "${LMISS_FILES[@]}" \
            >> "${LOG_DIR}/step4_plots.log" 2>&1
        echo "[OK]    Individual pass plots saved to  : ${STEP4_DIR}"
        echo "[OK]    Faceted comparison plot saved to: ${STEP4_DIR}"
    fi
    echo ""

    # Build pass summary string for report
    PASS_SUMMARY=""
    for i in "${!thresholds[@]}"; do
        PASS_SUMMARY="${PASS_SUMMARY}  Pass $((i+1)) (geno=${thresholds[$i]}) : ${REMOVED_PER_PASS[$i]} variants removed\n"
    done

    write_report "${mode}" "${GENO}" \
        "${VARIANTS_BEFORE}" "${VARIANTS_FINAL}" "${TOTAL_REMOVED}" \
        "${PASS_SUMMARY}" \
        "${CURRENT_INPUT}"
}

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
        echo "  CRISP -- STEP 4 VARIANT CALL RATE REPORT"
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
        echo "  VARIANT COUNTS"
        echo "  Variants before filtering : ${before}"
        echo "  Variants after filtering  : ${after}"
        echo "  Variants removed          : ${total_removed}"
        echo "  Samples (unchanged)       : ${SAMPLES}"
        if [ -n "${pass_summary}" ]; then
        echo "------------------------------------------------------------------"
        echo "  PASS SUMMARY"
        echo -e "${pass_summary}"
        fi
        echo "------------------------------------------------------------------"
        echo "  OUTPUT"
        echo "  Final BED prefix  : ${final_prefix}"
        echo "  Exclusion list    : ${OUTPUT_DIR}/${PROJECT_NAME}_exclusions_step4.txt"
        echo "  Plots             : ${STEP4_DIR}"
        echo "  Log               : ${LOG_DIR}"
        echo "------------------------------------------------------------------"
        echo "  PLOTS GENERATED"
        if [ "${mode}" = "SIMPLE" ]; then
            echo "  step4_snprate_simple.pdf"
        else
            echo "  step4_snprate_${mode,,}_pass[N].pdf (one per pass)"
            echo "  step4_snprate_${mode,,}_faceted.pdf"
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

echo "[CRISP] Step 4 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
