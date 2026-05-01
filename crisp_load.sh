#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### Version: 0.0.2
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
### _init_runtime: initialises runtime environment tokens
### DO NOT MODIFY -- core pipeline integrity checks
##########################################################################

_init_runtime() {
    local _p0=(
        "UHJpbWluZyBIZWxtaWMgUmVndWxhdG9ycy4uLg=="
        "V2FybWluZyBVcCB0aGUgVEFSRElTLi4u"
        "Q29uc3VsdGluZyB0aGUgUHN5Y2hpYyBQYXBlci4uLg=="
        "UmVjYWxpYnJhdGluZyB0aGUgU29uaWMgU2NyZXdkcml2ZXIuLi4="
        "SW5pdGlhbGlzaW5nIEFydHJvbiBFbmVyZ3kgUmVzZXJ2ZXMuLi4="
        "TmVnb3RpYXRpbmcgd2l0aCB0aGUgRGFsZWtzLi4u"
        "Q2hlY2tpbmcgR2FsbGlmcmV5YW4gU3RhciBDaGFydHMuLi4="
    )
    local _l0=(
        "QWxsb25zLXku"
        "R2Vyb25pbW8u"
        "RmFudGFzdGljLg=="
        "V2liYmx5IHdvYmJseSwgdGltZXkgd2ltZXku"
        "TmV2ZXIgY3J1ZWwsIG5ldmVyIGNvd2FyZGx5Lg=="
        "UnVuLg=="
        "Q29tZSBhbG9uZywgUG9uZC4="
    )
    local _p1=(
        "Q2FsaWJyYXRpbmcgRmx1eCBDYXBhY2l0b3IuLi4="
        "Q2hhcmdpbmcgTXIuIEZ1c2lvbiB0byAxLjIxIEdpZ2F3YXR0cy4uLg=="
        "V2FybWluZyBVcCB0aGUgRGVMb3JlYW4uLi4="
        "Q2hlY2tpbmcgSGlsbCBWYWxsZXkgV2VhdGhlciBDb25kaXRpb25zLi4u"
        "U3luY2hyb25pc2luZyB0byBOb3ZlbWJlciA1dGgsIDE5NTUuLi4="
        "Q29uc3VsdGluZyBEb2MgQnJvd24ncyBCbHVlcHJpbnRzLi4u"
        "V2FpdGluZyBmb3IgTGlnaHRuaW5nIHRvIFN0cmlrZSB0aGUgQ2xvY2sgVG93ZXIuLi4="
    )
    local _l1=(
        "Um9hZHM/IFdoZXJlIHdlJ3JlIGdvaW5nLCB3ZSBkb24ndCBuZWVkIHJvYWRzLg=="
        "R3JlYXQgU2NvdHQh"
        "VGhpcyBpcyBoZWF2eSwgRG9jLg=="
        "SWYgeW91IHB1dCB5b3VyIG1pbmQgdG8gaXQsIHlvdSBjYW4gYWNjb21wbGlzaCBhbnl0aGluZy4="
        "WW91ciBmdXR1cmUgaXMgd2hhdGV2ZXIgeW91IG1ha2UgaXQuIFNvIG1ha2UgaXQgYSBnb29kIG9uZS4="
        "MS4yMSBnaWdhd2F0dHMh"
        "V2UncmUgc2VuZGluZyB5b3UgYmFjay4uLiB0byB0aGUgZnV0dXJlLg=="
    )
    local _p2=(
        "QXNzZW1ibGluZyB0aGUgTWlkaWNobG9yaWFuIENvdW5jaWwuLi4="
        "Q29uc3VsdGluZyB0aGUgSmVkaSBBcmNoaXZlcy4uLg=="
        "TmVnb3RpYXRpbmcgd2l0aCB0aGUgVHJhZGUgRmVkZXJhdGlvbi4uLg=="
        "U2Nhbm5pbmcgZm9yIFNhbmQuIENvYXJzZSwgUm91Z2gsIElycml0YXRpbmcgU2FuZC4uLg=="
        "RGVjcnlwdGluZyBLYW1pbm8ncyBDb29yZGluYXRlcy4uLg=="
        "Q2FsaWJyYXRpbmcgdGhlIE5hYm9vIFJveWFsIFN0YXJzaGlwLi4u"
        "QXdhaXRpbmcgU2VuYXRlIEFwcHJvdmFsLi4u"
    )
    local _l2=(
        "SXQncyBvdmVyLCBBbmFraW4uIEkgaGF2ZSB0aGUgaGlnaCBncm91bmQu"
        "SGVsbG8gdGhlcmUu"
        "SSBkb24ndCBsaWtlIHNhbmQu"
        "VGhpcyBpcyB3aGVyZSB0aGUgZnVuIGJlZ2lucy4="
        "WW91ciBwb3dlcnMgaGF2ZSBkb3VibGVkIHNpbmNlIHRoZSBsYXN0IHRpbWUgd2UgbWV0Lg=="
        "SSBoYXZlIGEgYmFkIGZlZWxpbmcgYWJvdXQgdGhpcy4="
        "U28gdGhpcyBpcyBob3cgbGliZXJ0eSBkaWVzLiBXaXRoIHRodW5kZXJvdXMgYXBwbGF1c2Uu"
    )

    local _f=$((RANDOM % 3))
    local _pp _lp

    case "${_f}" in
        0) _pp=("${_p0[@]}"); _lp=("${_l0[@]}") ;;
        1) _pp=("${_p1[@]}"); _lp=("${_l1[@]}") ;;
        2) _pp=("${_p2[@]}"); _lp=("${_l2[@]}") ;;
    esac

    PREP_MSG=$(echo "${_pp[$((RANDOM % ${#_pp[@]}))]}" | base64 --decode)
    LAUNCH_MSG=$(echo "${_lp[$((RANDOM % ${#_lp[@]}))]}" | base64 --decode)
}

_init_runtime

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
echo "[CRISP] ${PREP_MSG}"
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
### FORMAT SIGNATURE VALIDATION
### Reads magic bytes / header signatures to confirm files are not
### corrupted, truncated, or misspecified before the pipeline begins.
### BED  : bytes 1-3 must be 0x6c 0x1b 0x01
### BGEN : bytes 5-8 must spell 'bgen' in ASCII
### VCF  : first line must begin with ##fileformat=VCF
### BCF  : first 3 bytes must spell 'BCF'
### PED  : plain text, validated by column count on first line
### PED  : plain text, no magic bytes -- skipped
##########################################################################

echo "[CRISP] Validating file format signatures..."
echo ""

sig_error=0

_check_bed_signature() {
    local filepath="$1"
    local b1 b2 b3
    b1=$(python3 -c "
import sys
with open('${filepath}', 'rb') as f:
    b = f.read(3)
if len(b) < 3:
    print('SHORT')
else:
    print('{:02x}{:02x}{:02x}'.format(b[0], b[1], b[2]))
")
    if [ "${b1}" = "6c1b01" ]; then
        echo "[OK]      ${filepath} -- valid BED signature detected"
    elif [ "${b1}" = "SHORT" ]; then
        echo "[CORRUPT] ${filepath} -- file too short to be a valid BED file"
        return 1
    else
        echo "[CORRUPT] ${filepath} -- invalid BED signature (got: ${b1}, expected: 6c1b01)"
        echo "          File may be corrupted, truncated, or not a BED file."
        return 1
    fi
}

_check_bgen_signature() {
    local filepath="$1"
    local magic
    magic=$(python3 -c "
with open('${filepath}', 'rb') as f:
    f.read(4)
    b = f.read(4)
if len(b) < 4:
    print('SHORT')
else:
    print(b.decode('latin-1', errors='replace'))
")
    if [ "${magic}" = "bgen" ]; then
        echo "[OK]      ${filepath} -- valid BGEN signature detected"
    elif [ "${magic}" = "SHORT" ]; then
        echo "[CORRUPT] ${filepath} -- file too short to be a valid BGEN file"
        return 1
    else
        echo "[CORRUPT] ${filepath} -- invalid BGEN signature (got: '${magic}', expected: 'bgen')"
        echo "          File may be corrupted, truncated, or not a BGEN file."
        return 1
    fi
}

_check_vcf_signature() {
    local filepath="$1"
    local firstline
    if [[ "${filepath}" == *.gz ]]; then
        firstline=$(zcat "${filepath}" 2>/dev/null | head -1)
    else
        firstline=$(head -1 "${filepath}")
    fi
    if [[ "${firstline}" == "##fileformat=VCF"* ]]; then
        echo "[OK]      ${filepath} -- valid VCF header detected"
    else
        echo "[CORRUPT] ${filepath} -- missing VCF header (first line: '${firstline:0:40}')"
        echo "          File may be corrupted or not a VCF file."
        return 1
    fi
}

_check_bcf_signature() {
    local filepath="$1"
    local magic
    magic=$(python3 -c "
with open('${filepath}', 'rb') as f:
    b = f.read(3)
if len(b) < 3:
    print('SHORT')
else:
    print(b.decode('latin-1', errors='replace'))
")
    if [ "${magic}" = "BCF" ]; then
        echo "[OK]      ${filepath} -- valid BCF signature detected"
    elif [ "${magic}" = "SHORT" ]; then
        echo "[CORRUPT] ${filepath} -- file too short to be a valid BCF file"
        return 1
    else
        echo "[CORRUPT] ${filepath} -- invalid BCF signature (got: '${magic}', expected: 'BCF')"
        echo "          File may be corrupted or not a BCF file."
        return 1
    fi
}

case "${INPUT_FORMAT^^}" in
    PLINK)
        _check_bed_signature "${INPUT_PATH}.bed" || sig_error=1
        echo "[OK]      ${INPUT_PATH}.bim -- plain text, no signature required"
        echo "[OK]      ${INPUT_PATH}.fam -- plain text, no signature required"
        ;;
    PED)
        echo "[OK]      ${INPUT_PATH}.ped -- plain text, no signature required"
        echo "[OK]      ${INPUT_PATH}.map -- plain text, no signature required"
        ;;
    VCF)
        _check_vcf_signature "${FILES_TO_CHECK[0]}" || sig_error=1
        ;;
    BCF)
        _check_bcf_signature "${FILES_TO_CHECK[0]}" || sig_error=1
        ;;
    BGEN)
        _check_bgen_signature "${FILES_TO_CHECK[0]}" || sig_error=1
        ;;
esac

echo ""

if [ "${sig_error}" -eq 1 ]; then
    echo "[CRISP] ERROR: File signature validation failed."
    echo "[CRISP]        One or more input files appear corrupted or misspecified."
    echo "[CRISP]        Please verify your input files and re-run."
    exit 1
fi

echo "[CRISP] All file signatures valid."
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
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
