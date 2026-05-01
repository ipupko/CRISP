#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Format Conversion
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
##########################################################################
### DESCRIPTION
### This script reads crisp_instructions.txt and converts input
### genotype files to PLINK BED/BIM/FAM format.
###
### Supported input formats:
###   PLINK  -- pass-through, creates working copy
###   PED    -- converted via PLINK 1.9 --make-bed
###   VCF    -- converted via PLINK 2 --vcf --make-bed
###   BCF    -- converted via PLINK 2 --bcf --make-bed
###   BGEN   -- converted via PLINK 2 --bgen --make-bed
###
### Chromosome cleaning applied post-conversion:
###   Retains chromosomes 1-22, X, Y, XY, MT only.
###   Non-standard chromosomes excluded and reported.
###
### Multiallelic variants excluded via --max-alleles 2
###   (VCF and BGEN only).
###
### Intermediate files retained or removed based on
###   KEEP_INTERMEDIATE in crisp_instructions.txt.
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
    echo "[CRISP]        Place crisp_instructions.txt in the same directory"
    echo "[CRISP]        as this script, or edit INSTRUCTION_FILE above."
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
INPUT_PATH=$(parse_param "INPUT_PATH")
OUTPUT_DIR=$(parse_param "OUTPUT_DIR" "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME" "project")
KEEP_INTERMEDIATE=$(parse_param "KEEP_INTERMEDIATE" "YES")
PLINK1=$(parse_param "PLINK1_PATH" "plink")
PLINK2=$(parse_param "PLINK2_PATH" "plink2")

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
echo "        INPUT_FORMAT      : ${INPUT_FORMAT}"
echo "        INPUT_PATH        : ${INPUT_PATH}"
echo "        OUTPUT_DIR        : ${OUTPUT_DIR}"
echo "        PROJECT_NAME      : ${PROJECT_NAME}"
echo "        KEEP_INTERMEDIATE : ${KEEP_INTERMEDIATE}"
echo "        PLINK1_PATH       : ${PLINK1}"
echo "        PLINK2_PATH       : ${PLINK2}"
echo ""

##########################################################################
### VALIDATE PARAMETERS
##########################################################################

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
### SET UP DIRECTORIES
##########################################################################

CONVERTED_DIR="${OUTPUT_DIR}/step2_converted"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_DIR="${OUTPUT_DIR}"

mkdir -p "${CONVERTED_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${REPORT_DIR}"

OUT_PREFIX="${CONVERTED_DIR}/${PROJECT_NAME}_converted"
CLEAN_PREFIX="${CONVERTED_DIR}/${PROJECT_NAME}_converted_chrclean"
LOG_FILE="${LOG_DIR}/step2_conversion.log"
REPORT_FILE="${REPORT_DIR}/${PROJECT_NAME}_step2_conversion_report.txt"
EXTRACT_LIST="${CONVERTED_DIR}/chr_extract_list.txt"

echo "[CRISP] Output directory : ${CONVERTED_DIR}"
echo "[CRISP] Log file         : ${LOG_FILE}"
echo ""

##########################################################################
### STEP 2A: FORMAT CONVERSION
##########################################################################

echo "##########################################################################"
echo "###                  STEP 2A: FORMAT CONVERSION                       ###"
echo "##########################################################################"
echo ""

case "${INPUT_FORMAT^^}" in

    PLINK)
        echo "[CRISP] Input is already BED/BIM/FAM. Creating working copy..."
        echo "[CRISP] Running: ${PLINK1} --bfile ${INPUT_PATH} --make-bed --allow-extra-chr --out ${OUT_PREFIX}"
        echo ""
        ${PLINK1} \
            --bfile "${INPUT_PATH}" \
            --make-bed \
            --allow-extra-chr \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: PLINK pass-through failed. Check log: ${LOG_FILE}"
            exit 1
        fi
        ;;

    PED)
        echo "[CRISP] Converting PED/MAP to BED/BIM/FAM using PLINK 1.9..."
        echo "[CRISP] Running: ${PLINK1} --ped ${INPUT_PATH}.ped --map ${INPUT_PATH}.map --make-bed --allow-extra-chr --out ${OUT_PREFIX}"
        echo ""
        ${PLINK1} \
            --ped "${INPUT_PATH}.ped" \
            --map "${INPUT_PATH}.map" \
            --make-bed \
            --allow-extra-chr \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: PED/MAP conversion failed. Check log: ${LOG_FILE}"
            exit 1
        fi
        ;;

    VCF)
        # Resolve VCF file path
        if [ -f "${INPUT_PATH}" ]; then
            VCF_FILE="${INPUT_PATH}"
        elif [ -f "${INPUT_PATH}.vcf.gz" ]; then
            VCF_FILE="${INPUT_PATH}.vcf.gz"
        elif [ -f "${INPUT_PATH}.vcf" ]; then
            VCF_FILE="${INPUT_PATH}.vcf"
        else
            echo "[CRISP] ERROR: No VCF file found at: ${INPUT_PATH}"
            exit 1
        fi

        echo "[CRISP] Converting VCF to BED/BIM/FAM using PLINK 2..."
        echo "[CRISP] Running: ${PLINK2} --vcf ${VCF_FILE} --make-bed --allow-extra-chr --max-alleles 2 --out ${OUT_PREFIX}"
        echo ""
        ${PLINK2} \
            --vcf "${VCF_FILE}" \
            --make-bed \
            --allow-extra-chr \
            --max-alleles 2 \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: VCF conversion failed. Check log: ${LOG_FILE}"
            exit 1
        fi
        ;;

    BCF)
        # Resolve BCF file path
        if [ -f "${INPUT_PATH}" ]; then
            BCF_FILE="${INPUT_PATH}"
        elif [ -f "${INPUT_PATH}.bcf" ]; then
            BCF_FILE="${INPUT_PATH}.bcf"
        else
            echo "[CRISP] ERROR: No BCF file found at: ${INPUT_PATH}"
            exit 1
        fi

        echo "[CRISP] Converting BCF to BED/BIM/FAM using PLINK 2..."
        echo "[CRISP] Running: ${PLINK2} --bcf ${BCF_FILE} --make-bed --allow-extra-chr --max-alleles 2 --out ${OUT_PREFIX}"
        echo ""
        ${PLINK2} \
            --bcf "${BCF_FILE}" \
            --make-bed \
            --allow-extra-chr \
            --max-alleles 2 \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: BCF conversion failed. Check log: ${LOG_FILE}"
            exit 1
        fi
        ;;

    BGEN)
        # Resolve BGEN file path
        if [ -f "${INPUT_PATH}" ]; then
            BGEN_FILE="${INPUT_PATH}"
        elif [ -f "${INPUT_PATH}.bgen" ]; then
            BGEN_FILE="${INPUT_PATH}.bgen"
        else
            echo "[CRISP] ERROR: No BGEN file found at: ${INPUT_PATH}"
            exit 1
        fi

        # Resolve companion .sample file
        SAMPLE_FILE="${BGEN_FILE%.bgen}.sample"
        if [ ! -f "${SAMPLE_FILE}" ]; then
            echo "[CRISP] ERROR: BGEN .sample file not found: ${SAMPLE_FILE}"
            exit 1
        fi

        echo "[CRISP] Converting BGEN to BED/BIM/FAM using PLINK 2..."
        echo "[CRISP] Running: ${PLINK2} --bgen ${BGEN_FILE} ref-first --sample ${SAMPLE_FILE} --make-bed --allow-extra-chr --max-alleles 2 --out ${OUT_PREFIX}"
        echo ""
        ${PLINK2} \
            --bgen "${BGEN_FILE}" ref-first \
            --sample "${SAMPLE_FILE}" \
            --make-bed \
            --allow-extra-chr \
            --max-alleles 2 \
            --out "${OUT_PREFIX}" \
            >> "${LOG_FILE}" 2>&1

        if [ $? -ne 0 ]; then
            echo "[CRISP] ERROR: BGEN conversion failed. Check log: ${LOG_FILE}"
            exit 1
        fi
        ;;

    *)
        echo "[CRISP] ERROR: Unrecognised INPUT_FORMAT: ${INPUT_FORMAT}"
        echo "[CRISP]        Valid options: PLINK, PED, VCF, BCF, BGEN"
        exit 1
        ;;
esac

echo "[OK]    Conversion complete."
echo ""

# Count variants before cleaning
VARIANTS_BEFORE=$(wc -l < "${OUT_PREFIX}.bim" | tr -d '[:space:]')
SAMPLES=$(wc -l < "${OUT_PREFIX}.fam" | tr -d '[:space:]')
echo "[CRISP] Variants after conversion : ${VARIANTS_BEFORE}"
echo "[CRISP] Samples                   : ${SAMPLES}"
echo ""

##########################################################################
### STEP 2B: CHROMOSOME CLEANING
### Retain only chromosomes 1-22, X, Y, XY, MT
### Exclude and report non-standard chromosome codes
##########################################################################

echo "##########################################################################"
echo "###               STEP 2B: CHROMOSOME CLEANING                        ###"
echo "##########################################################################"
echo ""

echo "[CRISP] Analysing chromosome composition..."
echo ""

# Extract unique chromosomes and count variants per chromosome
awk '{print $1}' "${OUT_PREFIX}.bim" | sort | uniq -c | sort -k2,2 > /tmp/crisp_chr_counts.txt

# Define valid chromosomes
VALID_CHR="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y XY MT chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrXY chrMT"

non_standard_found=0
excluded_total=0

while read -r count chrom; do
    is_valid=0
    for valid in ${VALID_CHR}; do
        if [ "${chrom}" = "${valid}" ]; then
            is_valid=1
            break
        fi
    done
    if [ "${is_valid}" -eq 0 ]; then
        echo "[CRISP WARNING] Non-standard chr '${chrom}' : ${count} variants -- will be excluded"
        non_standard_found=1
        excluded_total=$((excluded_total + count))
    fi
done < /tmp/crisp_chr_counts.txt

echo ""

if [ "${non_standard_found}" -eq 0 ]; then
    echo "[CRISP] No non-standard chromosomes found. Skipping chr cleaning step."
    CLEAN_PREFIX="${OUT_PREFIX}"
    VARIANTS_AFTER="${VARIANTS_BEFORE}"
else
    echo "[CRISP] Building variant extract list for standard chromosomes..."

    # Write extract list -- variant IDs on standard chromosomes only
    awk -v valid="${VALID_CHR}" '
    BEGIN { split(valid, arr, " "); for (i in arr) vset[arr[i]]=1 }
    { if ($1 in vset) print $2 }
    ' "${OUT_PREFIX}.bim" > "${EXTRACT_LIST}"

    VARIANTS_AFTER=$(wc -l < "${EXTRACT_LIST}" | tr -d '[:space:]')

    echo "[CRISP] Variants retained after chr cleaning : ${VARIANTS_AFTER}"
    echo "[CRISP] Variants excluded (non-standard chr) : ${excluded_total}"
    echo ""

    echo "[CRISP] Running: ${PLINK1} --bfile ${OUT_PREFIX} --extract ${EXTRACT_LIST} --make-bed --out ${CLEAN_PREFIX}"
    echo ""

    ${PLINK1} \
        --bfile "${OUT_PREFIX}" \
        --extract "${EXTRACT_LIST}" \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_FILE}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Chromosome cleaning failed. Check log: ${LOG_FILE}"
        exit 1
    fi

    echo "[OK]    Chromosome cleaning complete."
    echo ""

    # Clean up intermediate files if requested
    if [ "${KEEP_INTERMEDIATE^^}" = "NO" ]; then
        echo "[CRISP] KEEP_INTERMEDIATE = NO -- removing intermediate files..."
        for ext in .bed .bim .fam .log .nosex; do
            if [ -f "${OUT_PREFIX}${ext}" ]; then
                rm "${OUT_PREFIX}${ext}"
                echo "[CRISP]   Removed: ${OUT_PREFIX}${ext}"
            fi
        done
    else
        echo "[CRISP] KEEP_INTERMEDIATE = YES -- intermediate files retained."
    fi
    echo ""
fi

##########################################################################
### STEP 2C: WRITE CONVERSION REPORT
##########################################################################

echo "##########################################################################"
echo "###                  STEP 2C: CONVERSION REPORT                       ###"
echo "##########################################################################"
echo ""

{
    echo "=================================================================="
    echo "  CRISP -- STEP 2 CONVERSION REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Input format : ${INPUT_FORMAT}"
    echo "  Input path   : ${INPUT_PATH}"
    echo "------------------------------------------------------------------"
    echo "  CONVERSION"
    case "${INPUT_FORMAT^^}" in
        PLINK) echo "  Tool used    : PLINK 1.9 (pass-through)" ;;
        PED)   echo "  Tool used    : PLINK 1.9" ;;
        *)     echo "  Tool used    : PLINK 2" ;;
    esac
    echo "  Output dir   : ${CONVERTED_DIR}"
    echo "  Output prefix: $(basename "${CLEAN_PREFIX}")"
    echo "------------------------------------------------------------------"
    echo "  VARIANT COUNTS"
    echo "  Variants before chr cleaning : ${VARIANTS_BEFORE}"
    echo "  Variants after chr cleaning  : ${VARIANTS_AFTER}"
    echo "  Variants excluded            : ${excluded_total}"
    echo "  Samples                      : ${SAMPLES}"
    echo "------------------------------------------------------------------"
    if [ "${non_standard_found}" -eq 1 ]; then
        echo "  NON-STANDARD CHROMOSOMES EXCLUDED"
        while read -r count chrom; do
            is_valid=0
            for valid in ${VALID_CHR}; do
                if [ "${chrom}" = "${valid}" ]; then is_valid=1; break; fi
            done
            if [ "${is_valid}" -eq 0 ]; then
                echo "  chr '${chrom}' : ${count} variants excluded"
            fi
        done < /tmp/crisp_chr_counts.txt
        echo "------------------------------------------------------------------"
    fi
    echo "  OUTPUT FILES"
    for ext in .bed .bim .fam; do
        filepath="${CLEAN_PREFIX}${ext}"
        if [ -f "${filepath}" ]; then
            size=$(du -sh "${filepath}" | cut -f1)
            md5=$(md5sum "${filepath}" | awk '{print $1}')
            echo "  $(basename "${filepath}")"
            echo "  Size : ${size}  |  MD5 : ${md5}"
            echo ""
        fi
    done
    echo "------------------------------------------------------------------"
    echo "  INTERMEDIATE FILES"
    echo "  KEEP_INTERMEDIATE : ${KEEP_INTERMEDIATE}"
    echo "=================================================================="
    echo "  END OF REPORT"
    echo "=================================================================="
} | tee "${REPORT_FILE}"

echo ""
echo "[CRISP] Report written to: ${REPORT_FILE}"
echo ""
echo "[CRISP] Step 2 complete. Converted prefix: ${CLEAN_PREFIX}"
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
