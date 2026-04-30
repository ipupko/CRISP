#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: File Upload and Validation
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
##########################################################################
### DESCRIPTION
### This script reads crisp_instructions.txt and validates the input
### files specified by the user. It checks file existence, computes
### MD5 checksums, and reports basic dataset dimensions before any
### QC steps are performed.
##########################################################################

set -euo pipefail

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

INSTRUCTION_FILE="crisp_instructions.txt"

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Place crisp_instructions.txt in the same directory"
    echo "[CRISP]        as this script, or edit INSTRUCTION_FILE above."
    exit 1
fi

echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""

##########################################################################
### PARSE INSTRUCTION FILE
### Reads KEY = VALUE pairs. Ignores comment lines starting with #.
### Strips inline comments and leading/trailing whitespace.
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
INPUT_PATH=$(parse_param "INPUT_PATH")
OUTPUT_DIR=$(parse_param "OUTPUT_DIR" "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME" "project")

##########################################################################
### VALIDATE REQUIRED PARAMETERS
##########################################################################

echo "[CRISP] Parameters read from instruction file:"
echo "        INPUT_FORMAT  : ${INPUT_FORMAT}"
echo "        INPUT_PATH    : ${INPUT_PATH}"
echo "        OUTPUT_DIR    : ${OUTPUT_DIR}"
echo "        PROJECT_NAME  : ${PROJECT_NAME}"
echo ""

param_error=0

if [ -z "${INPUT_FORMAT}" ]; then
    echo "[CRISP] ERROR: INPUT_FORMAT is not set in ${INSTRUCTION_FILE}"
    param_error=1
fi

if [ -z "${INPUT_PATH}" ]; then
    echo "[CRISP] ERROR: INPUT_PATH is not set in ${INSTRUCTION_FILE}"
    param_error=1
fi

if [ "${param_error}" -eq 1 ]; then
    echo "[CRISP] Please edit ${INSTRUCTION_FILE} and re-run."
    exit 1
fi

##########################################################################
### RESOLVE INPUT FILES BASED ON FORMAT
### PLINK/PED : INPUT_PATH is a prefix. Appends expected extensions.
### VCF/BCF   : INPUT_PATH is a full file path.
### BGEN      : INPUT_PATH is a full file path. Also looks for .sample.
##########################################################################

echo "[CRISP] Resolving input files for format: ${INPUT_FORMAT}"
echo ""

declare -a FILES_TO_CHECK

case "${INPUT_FORMAT^^}" in
    PLINK)
        FILES_TO_CHECK=(
            "${INPUT_PATH}.bed"
            "${INPUT_PATH}.bim"
            "${INPUT_PATH}.fam"
        )
        ;;
    PED)
        FILES_TO_CHECK=(
            "${INPUT_PATH}.ped"
            "${INPUT_PATH}.map"
        )
        ;;
    VCF)
        if [ -f "${INPUT_PATH}" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}")
        elif [ -f "${INPUT_PATH}.vcf" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}.vcf")
        elif [ -f "${INPUT_PATH}.vcf.gz" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}.vcf.gz")
        else
            echo "[CRISP] ERROR: No VCF file found at: ${INPUT_PATH}"
            exit 1
        fi
        ;;
    BCF)
        if [ -f "${INPUT_PATH}" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}")
        elif [ -f "${INPUT_PATH}.bcf" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}.bcf")
        else
            echo "[CRISP] ERROR: No BCF file found at: ${INPUT_PATH}"
            exit 1
        fi
        ;;
    BGEN)
        if [ -f "${INPUT_PATH}" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}")
        elif [ -f "${INPUT_PATH}.bgen" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}.bgen")
        else
            echo "[CRISP] ERROR: No BGEN file found at: ${INPUT_PATH}"
            exit 1
        fi
        # Also check for companion .sample file
        SAMPLE_FILE="${INPUT_PATH%.bgen}.sample"
        if [ -f "${SAMPLE_FILE}" ]; then
            FILES_TO_CHECK+=("${SAMPLE_FILE}")
        else
            echo "[CRISP] WARNING: No .sample file found alongside BGEN."
            echo "                 Expected: ${SAMPLE_FILE}"
        fi
        ;;
    *)
        echo "[CRISP] ERROR: Unrecognised INPUT_FORMAT: ${INPUT_FORMAT}"
        echo "               Valid options: PLINK, PED, VCF, BCF, BGEN"
        exit 1
        ;;
esac

##########################################################################
### CHECK FILE EXISTENCE AND COMPUTE MD5 CHECKSUMS
##########################################################################

mkdir -p "${OUTPUT_DIR}"
MD5_REPORT="${OUTPUT_DIR}/${PROJECT_NAME}_input_md5.txt"
SUMMARY_REPORT="${OUTPUT_DIR}/${PROJECT_NAME}_input_summary.txt"

echo "[CRISP] Checking files and computing MD5 checksums..."
echo ""

missing=0

# Write MD5 report header
echo "# CRISP Input File MD5 Checksums" > "${MD5_REPORT}"
echo "# Project : ${PROJECT_NAME}" >> "${MD5_REPORT}"
echo "# Date    : $(date)" >> "${MD5_REPORT}"
echo "# Format  : ${INPUT_FORMAT}" >> "${MD5_REPORT}"
echo "" >> "${MD5_REPORT}"

for filepath in "${FILES_TO_CHECK[@]}"; do
    if [ ! -f "${filepath}" ]; then
        echo "[MISSING] ${filepath}"
        missing=$((missing + 1))
    else
        size=$(du -sh "${filepath}" | cut -f1)
        md5=$(md5sum "${filepath}" | awk '{print $1}')
        echo "[OK]      ${filepath}  (${size})  MD5: ${md5}"
        echo "${md5}  ${filepath}" >> "${MD5_REPORT}"
    fi
done

echo ""

if [ "${missing}" -gt 0 ]; then
    echo "[CRISP] ERROR: ${missing} expected file(s) not found."
    echo "[CRISP]        Check INPUT_PATH in ${INSTRUCTION_FILE} and re-run."
    exit 1
fi

echo "[CRISP] All input files found."
echo "[CRISP] MD5 checksums written to: ${MD5_REPORT}"
echo ""

##########################################################################
### REPORT DATASET DIMENSIONS
### Sample and variant counts extracted from PLINK files.
### For other formats a note is written to the summary.
##########################################################################

echo "[CRISP] Extracting dataset dimensions..."
echo ""

n_samples="N/A"
n_variants="N/A"

case "${INPUT_FORMAT^^}" in
    PLINK)
        n_samples=$(wc -l < "${INPUT_PATH}.fam" | tr -d '[:space:]')
        n_variants=$(wc -l < "${INPUT_PATH}.bim" | tr -d '[:space:]')
        ;;
    PED)
        n_samples=$(wc -l < "${INPUT_PATH}.ped" | tr -d '[:space:]')
        n_variants=$(wc -l < "${INPUT_PATH}.map" | tr -d '[:space:]')
        ;;
    VCF)
        n_samples=$(grep "^#CHROM" "${FILES_TO_CHECK[0]}" \
            | tr '\t' '\n' | tail -n +10 | wc -l | tr -d '[:space:]')
        n_variants=$(grep -v "^#" "${FILES_TO_CHECK[0]}" \
            | wc -l | tr -d '[:space:]')
        ;;
    BGEN)
        n_samples="See .sample file"
        n_variants="Requires bgenix"
        ;;
esac

echo "        Samples  : ${n_samples}"
echo "        Variants : ${n_variants}"
echo ""

##########################################################################
### WRITE SUMMARY REPORT
##########################################################################

{
    echo "=================================================================="
    echo "  CRISP -- INPUT SUMMARY REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Format       : ${INPUT_FORMAT}"
    echo "  Input path   : ${INPUT_PATH}"
    echo "  Output dir   : ${OUTPUT_DIR}"
    echo "------------------------------------------------------------------"
    echo "  FILES"
    for filepath in "${FILES_TO_CHECK[@]}"; do
        size=$(du -sh "${filepath}" | cut -f1)
        md5=$(md5sum "${filepath}" | awk '{print $1}')
        echo "    ${filepath}"
        echo "    Size : ${size}  |  MD5 : ${md5}"
        echo ""
    done
    echo "------------------------------------------------------------------"
    echo "  DIMENSIONS"
    echo "    Samples  : ${n_samples}"
    echo "    Variants : ${n_variants}"
    echo "=================================================================="
} > "${SUMMARY_REPORT}"

echo "[CRISP] Summary report written to: ${SUMMARY_REPORT}"
echo ""
echo "[CRISP] File upload and validation complete. Ready to proceed."
echo ""
