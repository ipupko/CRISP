#!/bin/bash
##########################################################################
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.   
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88. 
# 888           888   .d88'  888  Y88bo.       888   .d88' 
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'  
# 888           888`88b.     888       `"Y88b  888         
# `88b    ooo   888  `88b.   888  oo     .d8P  888         
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o        
#                                                          
#                                                          
#                                                          
##########################################################################
### Comprehensive Robust Integrated SNP Processing
### CHUNK : File Upload and Validation
### Version: 0.3.0-alpha
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
###
### WHAT THIS SCRIPT DOES
### ----------------------
### This is the entry point of the CRISP pipeline. Before any quality
### control or analysis can begin, we need to know that the input data
### is actually there, readable, and what it claims to be.
###
### This script performs six tasks in order:
###
###   1. PARAMETERS   -- reads user settings from the instruction file
###   2. DECOMPRESS   -- unpacks compressed input files if needed
###   3. INPUT FILES  -- resolves and locates all expected input files
###   4. CHECKSUMS    -- computes MD5 fingerprints for every input file
###   5. SIGNATURES   -- reads file magic bytes to confirm format integrity
###   6. DIMENSIONS   -- counts samples and variants in the dataset
###   7. COMPRESS     -- optionally gzips uncompressed input (interactive)
###   8. METADATA     -- writes a structured record for COMPASS-AI
###   9. REPORT       -- writes a human-readable input summary
###
### WHY THIS MATTERS
### ----------------
### Corrupt or misspecified input files are the most common cause of
### silent pipeline failures. Catching problems at intake -- before
### any PLINK or other tool runs -- saves hours of debugging downstream.
### MD5 checksums also create a reproducibility record: you can prove
### exactly which files were processed and that they were not modified.
###
### USAGE
### -----
###   bash crisp_upload_chunk.sh
###   bash crisp_upload_chunk.sh --config my_project.txt
###
### OPTIONS
###   --config <file>   Path to instruction file
###                     Default: crisp_instructions.txt
###
### v0.3.0-alpha: adds COMPASS-AI metadata collection.
###   Project name anonymised as CR_PID_XXXXXX in the registry.
###   Opt out via METADATA_COLLECTION = NO in the instruction file.
##########################################################################

set -euo pipefail
# set -e  : exit immediately if any command returns a non-zero exit code
# set -u  : treat unset variables as errors (prevents silent empty-string bugs)
# set -o pipefail : if any command in a pipe fails, the whole pipe fails
#                   (without this, only the last command's exit code is checked)

##########################################################################
#  ___       _ _   _       _ _          
# |_ _|_ __ (_) |_(_) __ _| (_)___  ___ 
#  | || '_ \| | __| |/ _` | | / __|/ _ \
#  | || | | | | |_| | (_| | | \__ \  __/
# |___|_| |_|_|\__|_|\__,_|_|_|___/\___|
#                                       
##########################################################################
### Load the CRISP flavour library.
###
### _crisp_flavour.sh is a shared runtime library sourced by every CRISP
### script. It provides:
###   - _init_runtime() : sets PREP_MSG and LAUNCH_MSG (franchise quotes)
###   - _crisp_palette(): exports CRISP_PAL_* colour variables for plots
###
### SCRIPT_DIR resolves to the directory containing this script,
### regardless of where the user calls it from. This ensures the
### relative path to scripts/_crisp_flavour.sh always works.
##########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime
# PREP_MSG and LAUNCH_MSG are now set -- franchise-themed progress messages
# printed at the start and end of the script respectively.

##########################################################################
#  ____                                _                
# |  _ \ __ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
# | |_) / _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# |  __/ (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
#                                                       
##########################################################################
### Parse command-line arguments.
###
### The only supported flag is --config, which lets the user point to
### a custom instruction file. This is important for running multiple
### projects from the same directory, or for SLURM array jobs where
### each job uses a different config.
###
### Default behaviour (no flags) looks for crisp_instructions.txt in
### the current working directory.
##########################################################################

INSTRUCTION_FILE="crisp_instructions.txt"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            INSTRUCTION_FILE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: bash crisp_upload_chunk.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_upload_chunk.sh [--config <file>]"
            exit 1
            ;;
    esac
done

# Confirm the instruction file exists before attempting to read it.
# A clear error here prevents confusing downstream failures.
if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_upload_chunk.sh --config <path>"
    exit 1
fi

echo ""
echo "[CRISP] ${PREP_MSG}"
echo ""
echo "[CRISP] Reading instruction file: ${INSTRUCTION_FILE}"
echo ""

##########################################################################
### parse_param() -- instruction file parser
###
### Reads a single KEY = VALUE pair from the instruction file.
### Rules:
###   - Lines starting with # are comments and are ignored
###   - Inline comments (anything after #) are stripped
###   - Leading and trailing whitespace is stripped from both key and value
###   - If the key is not found, the default value (arg $2) is returned
###   - head -1 ensures only the first match is used if a key appears twice
###
### Usage: variable=$(parse_param "KEY" "default_value")
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

# ── Read all parameters from the instruction file ─────────────────────
# Required parameters (no default -- pipeline will abort if missing):
INPUT_FORMAT=$(parse_param "INPUT_FORMAT")
INPUT_PATH=$(parse_param "INPUT_PATH")

# Optional parameters with safe defaults:
OUTPUT_DIR=$(parse_param "OUTPUT_DIR"           "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME"       "project")
SCHEDULER=$(parse_param "SCHEDULER"             "NONE")
REFERENCE_ASSEMBLY=$(parse_param "REFERENCE_ASSEMBLY" "GRCh38")
COHORT_POPULATION=$(parse_param "COHORT_POPULATION"   "UNKNOWN")
COHORT_SUBPOPULATION=$(parse_param "COHORT_SUBPOPULATION" "NONE")
COHORT_ADMIXED=$(parse_param "COHORT_ADMIXED"   "NO")
METADATA_COLLECTION=$(parse_param "METADATA_COLLECTION" "YES")

# ── Echo all resolved parameters to the terminal ──────────────────────
# This makes it easy to confirm settings at a glance before the pipeline
# runs, and provides a paper trail in log files.
echo "[CRISP] Parameters read from instruction file:"
echo "        INPUT_FORMAT         : ${INPUT_FORMAT}"
echo "        INPUT_PATH           : ${INPUT_PATH}"
echo "        OUTPUT_DIR           : ${OUTPUT_DIR}"
echo "        PROJECT_NAME         : ${PROJECT_NAME}"
echo "        REFERENCE_ASSEMBLY   : ${REFERENCE_ASSEMBLY}"
echo "        COHORT_POPULATION    : ${COHORT_POPULATION}"
echo "        COHORT_SUBPOPULATION : ${COHORT_SUBPOPULATION}"
echo "        COHORT_ADMIXED       : ${COHORT_ADMIXED}"
echo "        SCHEDULER            : ${SCHEDULER}"
echo "        METADATA_COLLECTION  : ${METADATA_COLLECTION}"
echo ""

# ── Validate required parameters ──────────────────────────────────────
# INPUT_FORMAT and INPUT_PATH are mandatory. Report all missing parameters
# before aborting so the user can fix everything in one edit.
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
    echo "[CRISP]        Edit ${INSTRUCTION_FILE} and re-run."
    exit 1
fi

##########################################################################
#  ____                                                   
# |  _ \  ___  ___ ___  _ __ ___  _ __  _ __ ___  ___ ___ 
# | | | |/ _ \/ __/ _ \| '_ ` _ \| '_ \| '__/ _ \/ __/ __|
# | |_| |  __/ (_| (_) | | | | | | |_) | | |  __/\__ \__ \
# |____/ \___|\___\___/|_| |_| |_| .__/|_|  \___||___/___/
#                                |_|                      
##########################################################################
### Auto-detect and handle compressed input files.
###
### Many biobanks and data providers deliver genotype data as gzipped
### archives. CRISP detects compressed input automatically and
### decompresses it to a working directory before proceeding.
###
### Supported compressed formats:
###   PLINK : .bed.gz, .bim.gz, .fam.gz
###   PED   : .ped.gz, .map.gz
###   VCF   : .vcf.gz is handled natively in the resolve step below
###   BCF   : binary compressed by design -- no action needed
###   BGEN  : binary compressed by design -- no action needed
###
### If both compressed and uncompressed versions exist for PLINK/PED,
### the uncompressed version is used and a warning is printed.
### WORKING_PATH tracks the effective prefix after any decompression,
### so all downstream steps reference the correct location.
##########################################################################

INPUT_COMPRESSED=0
DECOMP_DIR="${OUTPUT_DIR}/step1_decompressed"
WORKING_PATH="${INPUT_PATH}"
# WORKING_PATH starts as INPUT_PATH and is updated to DECOMP_DIR/basename
# if decompression occurs. All downstream file references use WORKING_PATH,
# not INPUT_PATH, so the rest of the script is transparent to compression.

# ── Helper: decompress a single file ──────────────────────────────────
# gunzip -c writes to stdout without deleting the source file.
# This preserves the original compressed archive.
_decompress_file() {
    local src="$1"
    local dst="$2"
    echo "[CRISP] Decompressing: $(basename ${src}) -> $(basename ${dst})"
    gunzip -c "${src}" > "${dst}"
    echo "[OK]    Done."
}

case "${INPUT_FORMAT^^}" in
    PLINK)
        if [ -f "${INPUT_PATH}.bed.gz" ] && [ ! -f "${INPUT_PATH}.bed" ]; then
            # Compressed PLINK triplet found; uncompressed version absent.
            # Decompress all three files (.bed, .bim, .fam) to DECOMP_DIR.
            echo "[CRISP] Compressed PLINK input detected (.bed.gz)"
            echo "[CRISP] Decompressing to: ${DECOMP_DIR}"
            echo ""
            mkdir -p "${DECOMP_DIR}"
            BASENAME=$(basename "${INPUT_PATH}")
            _decompress_file "${INPUT_PATH}.bed.gz" "${DECOMP_DIR}/${BASENAME}.bed"
            # .bim and .fam may or may not be compressed independently
            for ext in bim fam; do
                if [ -f "${INPUT_PATH}.${ext}.gz" ]; then
                    _decompress_file \
                        "${INPUT_PATH}.${ext}.gz" \
                        "${DECOMP_DIR}/${BASENAME}.${ext}"
                elif [ -f "${INPUT_PATH}.${ext}" ]; then
                    cp "${INPUT_PATH}.${ext}" "${DECOMP_DIR}/${BASENAME}.${ext}"
                else
                    echo "[CRISP] ERROR: Cannot find ${INPUT_PATH}.${ext} or ${INPUT_PATH}.${ext}.gz"
                    exit 1
                fi
            done
            WORKING_PATH="${DECOMP_DIR}/${BASENAME}"
            INPUT_COMPRESSED=1
            echo ""
        elif [ -f "${INPUT_PATH}.bed.gz" ] && [ -f "${INPUT_PATH}.bed" ]; then
            # Both exist -- use uncompressed to avoid ambiguity.
            echo "[CRISP] WARNING: Both .bed and .bed.gz found."
            echo "                 Using uncompressed: ${INPUT_PATH}.bed"
            echo ""
        fi
        ;;
    PED)
        if [ -f "${INPUT_PATH}.ped.gz" ] && [ ! -f "${INPUT_PATH}.ped" ]; then
            echo "[CRISP] Compressed PED input detected (.ped.gz)"
            mkdir -p "${DECOMP_DIR}"
            BASENAME=$(basename "${INPUT_PATH}")
            _decompress_file "${INPUT_PATH}.ped.gz" "${DECOMP_DIR}/${BASENAME}.ped"
            if [ -f "${INPUT_PATH}.map.gz" ]; then
                _decompress_file "${INPUT_PATH}.map.gz" "${DECOMP_DIR}/${BASENAME}.map"
            else
                cp "${INPUT_PATH}.map" "${DECOMP_DIR}/${BASENAME}.map"
            fi
            WORKING_PATH="${DECOMP_DIR}/${BASENAME}"
            INPUT_COMPRESSED=1
            echo ""
        fi
        ;;
    VCF)
        # .vcf.gz is handled natively in the resolve step below.
        # bcftools and tabix read .vcf.gz directly; no decompression needed.
        ;;
    BCF|BGEN)
        # BCF and BGEN are binary formats that are inherently compressed.
        # No decompression step is required or appropriate.
        ;;
esac

##########################################################################
#  ___                   _     _____ _ _           
# |_ _|_ __  _ __  _   _| |_  |  ___(_) | ___  ___ 
#  | || '_ \| '_ \| | | | __| | |_  | | |/ _ \/ __|
#  | || | | | |_) | |_| | |_  |  _| | | |  __/\__ \
# |___|_| |_| .__/ \__,_|\__| |_|   |_|_|\___||___/
#           |_|                                    
##########################################################################
### Resolve the full list of expected input files.
###
### PLINK and PED formats use a filename prefix: INPUT_PATH points to
### the base name without extension, and we append the expected extensions.
### For example, INPUT_PATH = /data/cohort resolves to:
###   /data/cohort.bed, /data/cohort.bim, /data/cohort.fam
###
### VCF, BCF, and BGEN formats use a full file path. We try the path
### as given first, then try appending the standard extension if the
### user omitted it.
###
### For BGEN, we also look for the companion .sample file, which contains
### sample IDs and is required for most downstream tools.
###
### FILES_TO_CHECK is the master list used by all subsequent sections.
##########################################################################

echo "[CRISP] Resolving input files for format: ${INPUT_FORMAT}"
echo ""

declare -a FILES_TO_CHECK

case "${INPUT_FORMAT^^}" in
    PLINK)
        FILES_TO_CHECK=(
            "${WORKING_PATH}.bed"   # binary genotype data
            "${WORKING_PATH}.bim"   # variant information (one row per SNP)
            "${WORKING_PATH}.fam"   # sample information (one row per individual)
        )
        ;;
    PED)
        FILES_TO_CHECK=(
            "${WORKING_PATH}.ped"   # genotype data in text format
            "${WORKING_PATH}.map"   # variant positions
        )
        ;;
    VCF)
        # Try the path as given, then with .vcf, then with .vcf.gz
        if [ -f "${INPUT_PATH}" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}")
        elif [ -f "${INPUT_PATH}.vcf" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}.vcf")
        elif [ -f "${INPUT_PATH}.vcf.gz" ]; then
            FILES_TO_CHECK=("${INPUT_PATH}.vcf.gz")
            INPUT_COMPRESSED=1
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
        # The .sample file carries sample IDs for BGEN data.
        # Without it, downstream tools cannot link genotypes to individuals.
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
        echo "[CRISP]        Valid options: PLINK, PED, VCF, BCF, BGEN"
        exit 1
        ;;
esac

##########################################################################
#   ____ _               _                            
#  / ___| |__   ___  ___| | _____ _   _ _ __ ___  ___ 
# | |   | '_ \ / _ \/ __| |/ / __| | | | '_ ` _ \/ __|
# | |___| | | |  __/ (__|   <\__ \ |_| | | | | | \__ \
#  \____|_| |_|\___|\___|_|\_\___/\__,_|_| |_| |_|___/
#                                                     
##########################################################################
### Check file existence and compute MD5 checksums.
###
### WHY MD5?
### MD5 produces a 128-bit fingerprint of a file's contents. If a file
### is corrupted, truncated, or replaced between pipeline runs, the MD5
### will change. Recording these hashes at intake creates a reproducibility
### audit trail: you can always prove which exact files were processed.
###
### MD5 is not cryptographically secure, but for file integrity checking
### in a genomics pipeline it is entirely sufficient and very fast.
###
### Two outputs:
###   {project}_input_md5.txt    : one line per file: <md5>  <path>
###                                compatible with md5sum --check
###   {project}_input_summary.txt: human-readable report (written later)
##########################################################################

mkdir -p "${OUTPUT_DIR}"
MD5_REPORT="${OUTPUT_DIR}/${PROJECT_NAME}_input_md5.txt"
SUMMARY_REPORT="${OUTPUT_DIR}/${PROJECT_NAME}_input_summary.txt"

echo "[CRISP] Checking files and computing MD5 checksums..."
echo ""

missing=0

{
    echo "# CRISP Input File MD5 Checksums"
    echo "# Project : ${PROJECT_NAME}"
    echo "# Date    : $(date)"
    echo "# Format  : ${INPUT_FORMAT}"
    echo ""
} > "${MD5_REPORT}"
# Using a brace group with a single redirect is more efficient than
# multiple >> appends -- the file is opened once and written sequentially.

for filepath in "${FILES_TO_CHECK[@]}"; do
    if [ ! -f "${filepath}" ]; then
        echo "[MISSING]  ${filepath}"
        missing=$((missing + 1))
    else
        size=$(du -sh "${filepath}" | cut -f1)
        md5=$(md5sum "${filepath}" | awk '{print $1}')
        echo "[OK]       ${filepath}  (${size})  MD5: ${md5}"
        echo "${md5}  ${filepath}" >> "${MD5_REPORT}"
    fi
done

echo ""

if [ "${missing}" -gt 0 ]; then
    echo "[CRISP] ERROR: ${missing} expected file(s) not found."
    echo "[CRISP]        Check INPUT_PATH in ${INSTRUCTION_FILE}."
    exit 1
fi

echo "[CRISP] All input files found."
echo "[CRISP] MD5 checksums written to: ${MD5_REPORT}"
echo ""

##########################################################################
#  ____  _                   _                       
# / ___|(_) __ _ _ __   __ _| |_ _   _ _ __ ___  ___ 
# \___ \| |/ _` | '_ \ / _` | __| | | | '__/ _ \/ __|
#  ___) | | (_| | | | | (_| | |_| |_| | | |  __/\__ \
# |____/|_|\__, |_| |_|\__,_|\__|\__,_|_|  \___||___/
#          |___/                                     
##########################################################################
### Validate file format signatures (magic bytes).
###
### WHY MAGIC BYTES?
### File extensions can lie. A file named .bed might be a text file,
### a different binary format, or a truncated download. Reading the
### first few bytes of a file -- the "magic bytes" or file signature --
### is the only reliable way to confirm a file is what it claims to be.
###
### Expected signatures:
###   BED  : bytes 1-3 = 0x6c 0x1b 0x01  (PLINK BED magic)
###   BGEN : bytes 5-8 = ASCII "bgen"
###   VCF  : first line starts with "##fileformat=VCF"
###   BCF  : bytes 1-3 = ASCII "BCF"
###   PED  : plain text -- no magic bytes, signature check skipped
##########################################################################

echo "[CRISP] Validating file format signatures..."
echo ""

sig_error=0

# ── BED signature check ───────────────────────────────────────────────
# Python reads the raw bytes more reliably than bash for binary files.
# The walrus operator (:=) is not used here for bash 3 compatibility.
_check_bed_signature() {
    local filepath="$1"
    local sig
    sig=$(python3 -c "
with open('${filepath}', 'rb') as f:
    b = f.read(3)
print('SHORT' if len(b) < 3 else '{:02x}{:02x}{:02x}'.format(b[0], b[1], b[2]))
")
    if [ "${sig}" = "6c1b01" ]; then
        echo "[OK]       ${filepath}  (valid BED signature)"
    elif [ "${sig}" = "SHORT" ]; then
        echo "[CORRUPT]  ${filepath}  (file too short)"
        return 1
    else
        echo "[CORRUPT]  ${filepath}  (bad signature: ${sig}, expected 6c1b01)"
        return 1
    fi
}

# ── BGEN signature check ──────────────────────────────────────────────
# BGEN stores its magic string at bytes 5-8 (after a 4-byte offset block).
_check_bgen_signature() {
    local filepath="$1"
    local magic
    magic=$(python3 -c "
with open('${filepath}', 'rb') as f:
    f.read(4)
    b = f.read(4)
print('SHORT' if len(b) < 4 else b.decode('latin-1', errors='replace'))
")
    if [ "${magic}" = "bgen" ]; then
        echo "[OK]       ${filepath}  (valid BGEN signature)"
    elif [ "${magic}" = "SHORT" ]; then
        echo "[CORRUPT]  ${filepath}  (file too short)"
        return 1
    else
        echo "[CORRUPT]  ${filepath}  (bad signature: '${magic}', expected 'bgen')"
        return 1
    fi
}

# ── VCF signature check ───────────────────────────────────────────────
# VCF files must begin with a meta-information header line.
# The ##fileformat line is the mandatory first line per VCF spec (v4.1+).
_check_vcf_signature() {
    local filepath="$1"
    local firstline
    if [[ "${filepath}" == *.gz ]]; then
        firstline=$(zcat "${filepath}" 2>/dev/null | head -1)
    else
        firstline=$(head -1 "${filepath}")
    fi
    if [[ "${firstline}" == "##fileformat=VCF"* ]]; then
        echo "[OK]       ${filepath}  (valid VCF header)"
    else
        echo "[CORRUPT]  ${filepath}  (missing VCF header, got: '${firstline:0:40}')"
        return 1
    fi
}

# ── BCF signature check ───────────────────────────────────────────────
# BCF is the binary version of VCF. Its magic bytes are ASCII "BCF".
_check_bcf_signature() {
    local filepath="$1"
    local magic
    magic=$(python3 -c "
with open('${filepath}', 'rb') as f:
    b = f.read(3)
print('SHORT' if len(b) < 3 else b.decode('latin-1', errors='replace'))
")
    if [ "${magic}" = "BCF" ]; then
        echo "[OK]       ${filepath}  (valid BCF signature)"
    elif [ "${magic}" = "SHORT" ]; then
        echo "[CORRUPT]  ${filepath}  (file too short)"
        return 1
    else
        echo "[CORRUPT]  ${filepath}  (bad signature: '${magic}', expected 'BCF')"
        return 1
    fi
}

# ── Dispatch to the appropriate check ─────────────────────────────────
case "${INPUT_FORMAT^^}" in
    PLINK)
        _check_bed_signature "${WORKING_PATH}.bed" || sig_error=1
        echo "[OK]       ${WORKING_PATH}.bim  (plain text)"
        echo "[OK]       ${WORKING_PATH}.fam  (plain text)"
        ;;
    PED)
        echo "[OK]       ${WORKING_PATH}.ped  (plain text)"
        echo "[OK]       ${WORKING_PATH}.map  (plain text)"
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
    echo "[CRISP] ERROR: One or more files failed signature validation."
    echo "[CRISP]        Check your input files and re-run."
    exit 1
fi

echo "[CRISP] All file signatures valid."
echo ""

##########################################################################
#  ____  _                          _                 
# |  _ \(_)_ __ ___   ___ _ __  ___(_) ___  _ __  ___ 
# | | | | | '_ ` _ \ / _ \ '_ \/ __| |/ _ \| '_ \/ __|
# | |_| | | | | | | |  __/ | | \__ \ | (_) | | | \__ \
# |____/|_|_| |_| |_|\___|_| |_|___/_|\___/|_| |_|___/
#                                                     
##########################################################################
### Extract dataset dimensions.
###
### Knowing the sample and variant counts before QC is important for
### two reasons:
###   1. It sets expectations -- if counts look wrong (e.g. far fewer
###      samples than expected), that indicates a file problem.
###   2. It provides a baseline for the QC reports: "X samples removed
###      out of Y total" is more informative than just "X removed".
###
### Method by format:
###   PLINK : wc -l on .fam (samples) and .bim (variants) -- instant
###   PED   : wc -l on .ped and .map
###   VCF   : grep #CHROM for sample count; grep -cv for variant count
###   BGEN  : variant count requires bgenix; noted in report
##########################################################################

echo "[CRISP] Extracting dataset dimensions..."
echo ""

n_samples="N/A"
n_variants="N/A"

case "${INPUT_FORMAT^^}" in
    PLINK)
        n_samples=$(wc -l < "${WORKING_PATH}.fam" | tr -d '[:space:]')
        n_variants=$(wc -l < "${WORKING_PATH}.bim" | tr -d '[:space:]')
        ;;
    PED)
        n_samples=$(wc -l < "${WORKING_PATH}.ped" | tr -d '[:space:]')
        n_variants=$(wc -l < "${WORKING_PATH}.map" | tr -d '[:space:]')
        ;;
    VCF)
        # The #CHROM header line lists all sample IDs starting at column 10.
        # Splitting by tab and counting from col 10 onward gives sample count.
        if [[ "${FILES_TO_CHECK[0]}" == *.gz ]]; then
            n_samples=$(zcat "${FILES_TO_CHECK[0]}" | grep "^#CHROM" \
                | tr '\t' '\n' | tail -n +10 | wc -l | tr -d '[:space:]')
            n_variants=$(zcat "${FILES_TO_CHECK[0]}" | grep -cv "^#" || true)
        else
            n_samples=$(grep "^#CHROM" "${FILES_TO_CHECK[0]}" \
                | tr '\t' '\n' | tail -n +10 | wc -l | tr -d '[:space:]')
            n_variants=$(grep -cv "^#" "${FILES_TO_CHECK[0]}" || true)
        fi
        ;;
    BGEN)
        # Sample count comes from the .sample file (2 header lines + 1 per sample).
        # Variant count requires bgenix and is deferred to avoid a hard dependency here.
        n_samples="See .sample file"
        n_variants="Run bgenix to count"
        ;;
esac

echo "        Samples  : ${n_samples}"
echo "        Variants : ${n_variants}"
echo ""

##########################################################################
#   ____                                        
#  / ___|___  _ __ ___  _ __  _ __ ___  ___ ___ 
# | |   / _ \| '_ ` _ \| '_ \| '__/ _ \/ __/ __|
# | |__| (_) | | | | | | |_) | | |  __/\__ \__ \
#  \____\___/|_| |_| |_| .__/|_|  \___||___/___/
#                      |_|                      
##########################################################################
### Offer to compress uncompressed input files (interactive only).
###
### Large PLINK datasets (multi-GB .bed files) can be stored at roughly
### 3-4x compression with gzip at negligible quality loss. This step
### offers to compress the files in-place after validation is complete.
###
### This prompt is skipped automatically in two cases:
###   1. SLURM / non-interactive runs (SCHEDULER != NONE)
###      -- interactive prompts hang batch jobs
###   2. Input was already compressed
###      -- no point re-compressing a decompressed working copy
###
### Only applies to PLINK and PED. VCF.gz, BCF, BGEN are already compact.
##########################################################################

if [ "${INPUT_COMPRESSED}" -eq 0 ] && \
   [[ "${INPUT_FORMAT^^}" =~ ^(PLINK|PED)$ ]] && \
   [ "${SCHEDULER^^}" = "NONE" ]; then

    # Calculate total size of input files in GB for the prompt message.
    total_bytes=0
    for filepath in "${FILES_TO_CHECK[@]}"; do
        fsize=$(du -sb "${filepath}" | cut -f1)
        total_bytes=$((total_bytes + fsize))
    done
    total_gb=$(python3 -c "print(f'{${total_bytes}/1024**3:.1f}')")

    echo "[CRISP] Input files are uncompressed."
    echo "        Total size: ${total_gb}G"
    echo ""
    echo "        Compress input files with gzip to save storage?"
    echo "        Original files will be replaced with .gz versions."
    echo "        This does not affect pipeline execution."
    echo ""
    read -rp "        Compress? [y/N]: " compress_answer

    if [[ "${compress_answer,,}" == "y" || "${compress_answer,,}" == "yes" ]]; then
        echo ""
        echo "[CRISP] Compressing input files..."
        for filepath in "${FILES_TO_CHECK[@]}"; do
            echo "        gzip ${filepath}"
            gzip -v "${filepath}" 2>&1 | sed 's/^/        /'
        done
        echo ""
        echo "[OK]    Compression complete."
        echo "[CRISP] Note: Update INPUT_PATH in ${INSTRUCTION_FILE} if needed"
        echo "              for future runs referencing the .gz files."
    else
        echo ""
        echo "[CRISP] Skipping compression. Input files unchanged."
    fi
    echo ""

elif [ "${INPUT_COMPRESSED}" -eq 0 ] && \
     [[ "${INPUT_FORMAT^^}" =~ ^(PLINK|PED)$ ]] && \
     [ "${SCHEDULER^^}" != "NONE" ]; then
    echo "[CRISP] Running under ${SCHEDULER} -- compression prompt skipped."
    echo "[CRISP] To compress input files, run crisp_upload_chunk.sh locally"
    echo "        with SCHEDULER = NONE."
    echo ""
fi

##########################################################################
#  __  __      _            _       _        
# |  \/  | ___| |_ __ _  __| | __ _| |_ __ _ 
# | |\/| |/ _ \ __/ _` |/ _` |/ _` | __/ _` |
# | |  | |  __/ || (_| | (_| | (_| | || (_| |
# |_|  |_|\___|\__\__,_|\__,_|\__,_|\__\__,_|
#                                            
##########################################################################
### COMPASS-AI metadata collection.
###
### WHAT IS THIS?
### COMPASS-AI is the learned GWAS engine being developed as part of the
### Compass Genomics platform. To train effectively, it needs structured
### information about every dataset that CRISP has ever processed:
### format, ancestry, size, assembly, and provenance.
###
### This block captures those attributes as a structured JSON record
### immediately after intake, when all the relevant information is known.
###
### TWO OUTPUTS
### -----------
### 1. {project}_cohort_metadata.json  -- full record, stored locally.
###    Contains PROJECT_NAME and file paths. Never shared externally.
###
### 2. crisp_cohort_registry.jsonl     -- anonymised record, append-only.
###    PROJECT_NAME is replaced with a random CR_PID_XXXXXX identifier.
###    File paths are stripped. Safe for aggregation across cohorts.
###    Appends one line per run, so the file grows as CRISP is used.
###
### ANONYMISATION
### -------------
### The CR_PID (Pipeline ID) is a random 6-character hex string, e.g.
### CR_PID_3A7F2C. It is generated fresh each run using Python's
### random.choices() over the hex alphabet. There is no mapping stored
### between CR_PID and PROJECT_NAME anywhere in this script -- that
### mapping exists only in {project}_cohort_metadata.json locally.
###
### OPT OUT
### -------
### Set METADATA_COLLECTION = NO in the instruction file to skip this
### block entirely. No JSON, no JSONL, no CR_PID generated.
##########################################################################

if [ "${METADATA_COLLECTION^^}" = "YES" ]; then

    echo "[CRISP] Collecting COMPASS-AI metadata..."

    # Generate a random anonymised pipeline ID.
    # Format: CR_PID_ followed by 6 uppercase hex characters (0-9, A-F).
    # Example: CR_PID_3A7F2C
    # Python's random.choices is used because bash $RANDOM only gives 0-32767
    # and formatting as zero-padded hex is cumbersome in pure bash.
    CRISP_PID="CR_PID_$(python3 -c "
import random
print(''.join(random.choices('0123456789ABCDEF', k=6)))
")"

    # ── Build a JSON array of file records ────────────────────────────
    # Each entry captures extension, size in MB, and MD5 hash.
    # File paths are intentionally excluded from this array so that the
    # same structure can be safely written to the anonymised JSONL registry.
    file_json="["
    first_file=1
    for filepath in "${FILES_TO_CHECK[@]}"; do
        size_bytes=$(du -sb "${filepath}" | cut -f1)
        size_mb=$(python3 -c "print(round(${size_bytes}/1024/1024, 3))")
        md5=$(md5sum "${filepath}" | awk '{print $1}')
        ext="${filepath##*.}"   # extract extension: everything after last dot
        if [ "${first_file}" -eq 0 ]; then file_json="${file_json},"; fi
        file_json="${file_json}{\"extension\":\"${ext}\",\"size_mb\":${size_mb},\"md5\":\"${md5}\"}"
        first_file=0
    done
    file_json="${file_json}]"

    # ── Compute total input size across all files ─────────────────────
    total_bytes=0
    for filepath in "${FILES_TO_CHECK[@]}"; do
        fsize=$(du -sb "${filepath}" | cut -f1)
        total_bytes=$((total_bytes + fsize))
    done
    total_mb=$(python3 -c "print(round(${total_bytes}/1024/1024, 3))")

    # ── ISO 8601 timestamp in UTC ─────────────────────────────────────
    # UTC avoids timezone ambiguity when records are aggregated across
    # institutions and compute environments.
    RUN_TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    # ── Write full local record (includes project name) ───────────────
    METADATA_JSON="${OUTPUT_DIR}/${PROJECT_NAME}_cohort_metadata.json"
    python3 - <<PYEOF
import json

record = {
    "schema_version"      : "1.0",
    "crisp_version"       : "0.3.0-alpha",
    "crisp_pid"           : "${CRISP_PID}",
    "project_name"        : "${PROJECT_NAME}",
    "run_timestamp"       : "${RUN_TIMESTAMP}",
    "input_format"        : "${INPUT_FORMAT}",
    "reference_assembly"  : "${REFERENCE_ASSEMBLY}",
    "cohort_population"   : "${COHORT_POPULATION}",
    "cohort_subpopulation": "${COHORT_SUBPOPULATION}",
    "cohort_admixed"      : "${COHORT_ADMIXED}",
    "scheduler"           : "${SCHEDULER}",
    "input_compressed"    : bool(${INPUT_COMPRESSED}),
    "n_samples"           : "${n_samples}",
    "n_variants"          : "${n_variants}",
    "total_input_size_mb" : ${total_mb},
    "files"               : ${file_json},
}

with open("${METADATA_JSON}", "w") as f:
    json.dump(record, f, indent=2)

print(f"[CRISP] Cohort metadata written to : ${METADATA_JSON}")
PYEOF

    # ── Write anonymised registry entry (no project name, no paths) ───
    # This record is safe to aggregate. The CR_PID links it back to the
    # full local record for anyone who needs to trace provenance.
    REGISTRY_JSONL="${OUTPUT_DIR}/crisp_cohort_registry.jsonl"
    python3 - <<PYEOF
import json

record = {
    "schema_version"      : "1.0",
    "crisp_version"       : "0.3.0-alpha",
    "crisp_pid"           : "${CRISP_PID}",
    "run_timestamp"       : "${RUN_TIMESTAMP}",
    "input_format"        : "${INPUT_FORMAT}",
    "reference_assembly"  : "${REFERENCE_ASSEMBLY}",
    "cohort_population"   : "${COHORT_POPULATION}",
    "cohort_subpopulation": "${COHORT_SUBPOPULATION}",
    "cohort_admixed"      : "${COHORT_ADMIXED}",
    "scheduler"           : "${SCHEDULER}",
    "input_compressed"    : bool(${INPUT_COMPRESSED}),
    "n_samples"           : "${n_samples}",
    "n_variants"          : "${n_variants}",
    "total_input_size_mb" : ${total_mb},
    "n_files"             : len(${file_json}),
}

with open("${REGISTRY_JSONL}", "a") as f:
    f.write(json.dumps(record) + "\n")

print(f"[CRISP] Registry entry appended to  : ${REGISTRY_JSONL}")
print(f"[CRISP] Anonymised pipeline ID      : ${CRISP_PID}")
PYEOF

    echo ""

else
    echo "[CRISP] METADATA_COLLECTION = NO -- skipping COMPASS-AI metadata."
    echo ""
fi

##########################################################################
#  ____                       _   
# |  _ \ ___ _ __   ___  _ __| |_ 
# | |_) / _ \ '_ \ / _ \| '__| __|
# |  _ <  __/ |_) | (_) | |  | |_ 
# |_| \_\___| .__/ \___/|_|   \__|
#           |_|                   
##########################################################################
### Write the human-readable input summary report.
###
### This is a plain-text file intended for the researcher to review
### before proceeding with QC. It consolidates everything captured
### in this script: file inventory, MD5s, sizes, dimensions, and
### key pipeline parameters as read from the instruction file.
###
### The summary is also consumed by crisp_report.sh (Step 11) when
### assembling the final pipeline report.
##########################################################################

{
    echo "=================================================================="
    echo "  CRISP INPUT SUMMARY REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project           : ${PROJECT_NAME}"
    echo "  Date              : $(date)"
    echo "  CRISP version     : 0.3.0-alpha"
    echo "  Format            : ${INPUT_FORMAT}"
    echo "  Input path        : ${INPUT_PATH}"
    echo "  Output dir        : ${OUTPUT_DIR}"
    echo "  Reference assembly: ${REFERENCE_ASSEMBLY}"
    echo "  Cohort population : ${COHORT_POPULATION}"
    echo "  Input compressed  : $([ ${INPUT_COMPRESSED} -eq 1 ] && echo YES || echo NO)"
    if [ "${INPUT_COMPRESSED}" -eq 1 ]; then
    echo "  Decompressed to   : ${DECOMP_DIR}"
    fi
    if [ "${METADATA_COLLECTION^^}" = "YES" ]; then
    echo "  COMPASS-AI PID    : ${CRISP_PID}"
    fi
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
echo "[CRISP] Step 1 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
