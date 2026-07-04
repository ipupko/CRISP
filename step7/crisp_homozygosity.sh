#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Heterozygosity and Homozygosity (Step 6)
### Version: 0.2.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Performs heterozygosity and homozygosity QC using PLINK output.
# Reuses the .sexcheck file from the sex check step for X chromosome data.
# The sex check step may be Step 3, 4 or 5 depending on STEP_ORDER.
#
# Steps:
#   6a: PLINK --het (autosomal heterozygosity)
#   6b: PLINK --homozyg (runs of homozygosity)
#   6c: X chromosome het check if no sexcheck file found
#   6d: Detection and plotting via plot_homozygosity.R or .py
#   6e: Count outliers
#   6f: Apply exclusions (optional)
#   6g: Write Step 6 report
#
# Outputs:
#   Homozygosity_Runs_vs_Zscore.pdf
#   report_homozygosity.txt
#   outliers_high_het.txt
#   outliers_excess_homo.txt
#   outliers_combined.txt
#
# Thresholds (overridable in crisp_instructions.txt):
#   HOM_Z_HIGH : Z-score upper threshold (excess homozygosity)
#   HOM_Z_LOW  : Z-score lower threshold (high heterozygosity)
#
# Usage:
#   bash crisp_homozygosity.sh
#   bash crisp_homozygosity.sh --config my_project.txt
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
            echo "Usage: bash crisp_homozygosity.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_homozygosity.sh [--config <file>]"
            exit 1
            ;;
    esac
done

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_homozygosity.sh --config <path>"
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
HOM_Z_HIGH=$(parse_param "HOM_Z_HIGH" "3")
HOM_Z_LOW=$(parse_param "HOM_Z_LOW" "-2")
HOMO_EXCLUDE=$(parse_param "HOMO_EXCLUDE" "YES")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# resolve input, use most recent clean output from sex check step,
# falling back through step4 and step2 if not found
STEP5_CLEAN=$(ls -t "${OUTPUT_DIR}/step5_sexcheck/${PROJECT_NAME}"*clean*.bed \
    2>/dev/null | head -1 | sed 's/\.bed$//')

if [ -n "${STEP5_CLEAN}" ] && [ -f "${STEP5_CLEAN}.bed" ]; then
    INPUT_PREFIX="${STEP5_CLEAN}"
    echo "[CRISP] Using sex check step clean output as input."
else
    INPUT_PREFIX=$(ls -t "${OUTPUT_DIR}/step4_snprate/${PROJECT_NAME}"*maf*.bed \
        2>/dev/null | head -1 | sed 's/\.bed$//')
    if [ -z "${INPUT_PREFIX}" ]; then
        INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
    fi
    echo "[CRISP] Using Step 4 output as input."
fi

# reuse sexcheck file from whichever step ran sex check
# searches across step3, step4 and step5 output dirs to support any STEP_ORDER
SEXCHECK_FILE=""
for step_dir in step5_sexcheck step4_snprate step3_callrate; do
    candidate="${OUTPUT_DIR}/${step_dir}/${PROJECT_NAME}_sexcheck.sexcheck"
    if [ -f "${candidate}" ]; then
        SEXCHECK_FILE="${candidate}"
        break
    fi
done

##########################################################################
### HEADER
##########################################################################

echo ""
echo "[CRISP] ${PREP_MSG}"
echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""
echo "[CRISP] Parameters:"
echo "        HOM_Z_HIGH   : ${HOM_Z_HIGH}"
echo "        HOM_Z_LOW    : ${HOM_Z_LOW}"
echo "        HOMO_EXCLUDE : ${HOMO_EXCLUDE}"
echo "        INPUT_PREFIX : ${INPUT_PREFIX}"
echo "        PLOT_ENGINE  : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE INPUT
##########################################################################

for ext in .bed .bim .fam; do
    if [ ! -f "${INPUT_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${INPUT_PREFIX}${ext}"
        echo "[CRISP]        Ensure the sex check step has been run first."
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

# check sexcheck file from whichever step produced it
if [ -z "${SEXCHECK_FILE}" ] || [ ! -f "${SEXCHECK_FILE}" ]; then
    echo "[CRISP] WARNING: No sexcheck file found in step3/4/5 output directories."
    echo "[CRISP]          Will run a fresh X chromosome check."
    RUN_X_HET=1
else
    echo "[OK]    Reusing sexcheck file: ${SEXCHECK_FILE}"
    RUN_X_HET=0
fi
echo ""

##########################################################################
### SET UP DIRECTORIES
##########################################################################

STEP6_DIR="${OUTPUT_DIR}/step6_homozygosity"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step6_homozygosity_report.txt"

mkdir -p "${STEP6_DIR}"
mkdir -p "${LOG_DIR}"

# per-run temp folder, avoids /tmp collisions on HPC shared nodes
CRISP_TMP="${OUTPUT_DIR}/.crisp_tmp_$$"
mkdir -p "${CRISP_TMP}"
trap "rm -rf ${CRISP_TMP}" EXIT

##########################################################################
### STEP 6a: AUTOSOMAL HETEROZYGOSITY
##########################################################################

echo "##########################################################################"
echo "###              STEP 6a: AUTOSOMAL HETEROZYGOSITY                    ###"
echo "##########################################################################"
echo ""

HET_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_het"
LOG_6A="${LOG_DIR}/step6a_het.log"

${PLINK1} \
    --bfile "${INPUT_PREFIX}" \
    --het \
    --out "${HET_PREFIX}" \
    >> "${LOG_6A}" 2>&1

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: PLINK --het failed. Check log: ${LOG_6A}"
    exit 1
fi

echo "[OK]    Heterozygosity file generated: ${HET_PREFIX}.het"
echo ""

##########################################################################
### STEP 6b: RUNS OF HOMOZYGOSITY
##########################################################################

echo "##########################################################################"
echo "###              STEP 6b: RUNS OF HOMOZYGOSITY                        ###"
echo "##########################################################################"
echo ""

HOH_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_homozyg"
LOG_6B="${LOG_DIR}/step6b_homozyg.log"

${PLINK1} \
    --bfile "${INPUT_PREFIX}" \
    --homozyg \
    --out "${HOH_PREFIX}" \
    >> "${LOG_6B}" 2>&1

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: PLINK --homozyg failed. Check log: ${LOG_6B}"
    exit 1
fi

echo "[OK]    ROH file generated: ${HOH_PREFIX}.hom.indiv"
echo ""

##########################################################################
### STEP 6c: X CHROMOSOME HET (FALLBACK IF STEP 5 FILE NOT FOUND)
##########################################################################

if [ "${RUN_X_HET}" -eq 1 ]; then

    echo "##########################################################################"
    echo "###           STEP 6c: X CHROMOSOME HET (FALLBACK)                   ###"
    echo "##########################################################################"
    echo ""

    XCHR_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_xchr_het"
    LOG_6C="${LOG_DIR}/step6c_xchr.log"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --chr X \
        --check-sex \
        --out "${XCHR_PREFIX}" \
        >> "${LOG_6C}" 2>&1

    SEXCHECK_FILE="${XCHR_PREFIX}.sexcheck"
    echo "[OK]    X chromosome sexcheck generated: ${SEXCHECK_FILE}"
    echo ""
fi

##########################################################################
### STEP 6d: DETECTION AND PLOTTING
##########################################################################

echo "##########################################################################"
echo "###              STEP 6d: DETECTION AND PLOTTING                      ###"
echo "##########################################################################"
echo ""

LOG_6D="${LOG_DIR}/step6d_plots.log"

if [ "${PLOT_ENGINE^^}" = "R" ]; then
    echo "[CRISP] Running plot_homozygosity.R..."
    ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_homozygosity.R" \
        "${STEP6_DIR}" \
        "${SEXCHECK_FILE}" \
        "${HET_PREFIX}.het" \
        "${HOH_PREFIX}.hom.indiv" \
        "${HOM_Z_HIGH}" \
        "${HOM_Z_LOW}" \
        >> "${LOG_6D}" 2>&1

elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
    echo "[CRISP] Running plot_homozygosity.py..."
    python3 "${SCRIPT_DIR}/scripts/plot_homozygosity.py" \
        "${STEP6_DIR}" \
        "${SEXCHECK_FILE}" \
        "${HET_PREFIX}.het" \
        "${HOH_PREFIX}.hom.indiv" \
        "${HOM_Z_HIGH}" \
        "${HOM_Z_LOW}" \
        >> "${LOG_6D}" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Plotting step failed. Check log: ${LOG_6D}"
    exit 1
fi

echo "[OK]    Detection and plotting complete."
echo ""

##########################################################################
### STEP 6e: COUNT OUTLIERS
##########################################################################

HET_EXCL="${STEP6_DIR}/outliers_high_het.txt"
HOMO_EXCL="${STEP6_DIR}/outliers_excess_homo.txt"
COMBINED_EXCL="${STEP6_DIR}/outliers_combined.txt"

N_HIGH_HET=0
N_EXCESS_HOMO=0
N_COMBINED=0

[ -f "${HET_EXCL}" ]      && N_HIGH_HET=$(wc -l    < "${HET_EXCL}"      | tr -d '[:space:]')
[ -f "${HOMO_EXCL}" ]     && N_EXCESS_HOMO=$(wc -l  < "${HOMO_EXCL}"     | tr -d '[:space:]')
[ -f "${COMBINED_EXCL}" ] && N_COMBINED=$(wc -l     < "${COMBINED_EXCL}" | tr -d '[:space:]')

echo "[CRISP] High heterozygosity outliers : ${N_HIGH_HET}"
echo "[CRISP] Excess homozygosity outliers : ${N_EXCESS_HOMO}"
echo "[CRISP] Combined outliers            : ${N_COMBINED}"
echo ""

##########################################################################
### STEP 6f: APPLY EXCLUSIONS (OPTIONAL)
##########################################################################

CLEAN_PREFIX="${INPUT_PREFIX}"
SAMPLES_AFTER="${SAMPLES}"

if [ "${HOMO_EXCLUDE^^}" = "YES" ] && \
   [ -f "${COMBINED_EXCL}" ] && [ -s "${COMBINED_EXCL}" ]; then

    echo "[CRISP] Applying exclusions..."

    CLEAN_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_step6_clean"
    LOG_EXCL="${LOG_DIR}/step6_exclusions.log"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --remove "${COMBINED_EXCL}" \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_EXCL}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Exclusion step failed. Check log: ${LOG_EXCL}"
        exit 1
    fi

    SAMPLES_AFTER=$(wc -l < "${CLEAN_PREFIX}.fam" | tr -d '[:space:]')
    echo "[OK]    Exclusions applied."
    echo "[CRISP] Samples before : ${SAMPLES}"
    echo "[CRISP] Samples removed: ${N_COMBINED}"
    echo "[CRISP] Samples after  : ${SAMPLES_AFTER}"
    echo ""

else
    echo "[CRISP] HOMO_EXCLUDE = NO, outliers retained in dataset."
    echo ""
fi

##########################################################################
### STEP 6g: WRITE STEP 6 REPORT
##########################################################################

{
    echo "=================================================================="
    echo "  CRISP: STEP 6 HETEROZYGOSITY AND HOMOZYGOSITY REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Input prefix : ${INPUT_PREFIX}"
    echo "------------------------------------------------------------------"
    echo "  PARAMETERS"
    echo "  HOM_Z_HIGH   : ${HOM_Z_HIGH}  (excess homozygosity threshold)"
    echo "  HOM_Z_LOW    : ${HOM_Z_LOW}  (high heterozygosity threshold)"
    echo "  HOMO_EXCLUDE : ${HOMO_EXCLUDE}"
    echo "------------------------------------------------------------------"
    echo "  OUTLIER COUNTS"
    echo "  High heterozygosity (Z < ${HOM_Z_LOW})  : ${N_HIGH_HET}"
    echo "  Excess homozygosity (Z > ${HOM_Z_HIGH}) : ${N_EXCESS_HOMO}"
    echo "  Combined outliers                        : ${N_COMBINED}"
    echo "------------------------------------------------------------------"
    echo "  SAMPLE COUNTS"
    echo "  Samples before exclusion : ${SAMPLES}"
    echo "  Samples after exclusion  : ${SAMPLES_AFTER}"
    echo "------------------------------------------------------------------"
    echo "  OUTPUT FILES"
    echo "  Plot         : ${STEP6_DIR}/Homozygosity_Runs_vs_Zscore.pdf"
    echo "  Full report  : ${STEP6_DIR}/report_homozygosity.txt"
    echo "  High het     : ${HET_EXCL}"
    echo "  Excess homo  : ${HOMO_EXCL}"
    echo "  Combined     : ${COMBINED_EXCL}"
    if [ "${HOMO_EXCLUDE^^}" = "YES" ]; then
    echo "  Clean prefix : ${CLEAN_PREFIX}"
    fi
    echo "  Log          : ${LOG_DIR}"
    echo "=================================================================="
    echo "  END OF REPORT"
    echo "=================================================================="
} | tee "${REPORT_FILE}"

echo ""
echo "[CRISP] Report written to: ${REPORT_FILE}"
echo ""
echo "[CRISP] Step 6 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
