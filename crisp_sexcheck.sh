#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Sex Check and Aneuploidy Detection (Step 5)
###
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Date Created : 16/07/2024
### Date Updated : 03/05/2026
##########################################################################
### DESCRIPTION
### Performs sex verification and chromosomal aneuploidy detection
### using PLINK F-statistics and Y chromosome variant counts.
###
### Steps:
###   5a: PLINK --check-sex (generates F-statistics per sample)
###   5b: PLINK Y chromosome count (variant missingness on chrY)
###   5c: Merge F-stat and Y-count data into single input for R
###   5d: Run plot_sexcheck.R for detection and plotting
###   5e: Write exclusion lists (mismatches and aneuploidies separate)
###   5f: Write Step 5 report
###
### Outputs:
###   Sex_check.pdf                      Three scatter plots
###   report_sex.mismatch_aneuploidies.txt  Detailed text report
###   id_list_sex_mismatch.txt           Sex mismatch exclusions
###   id_list_aneuploidies.txt           Aneuploidy exclusions
###
### Thresholds (all overridable in crisp_instructions.txt):
###   SEX_F_FEMALE_MAX   F threshold for female classification
###   SEX_F_MALE_MIN     F threshold for male classification
###   SEX_F_TURNER       F threshold for Turner (X0) detection
###   SEX_F_KLINEFELTER  F threshold for Klinefelter (XXY) detection
###   SEX_F_XXX          F threshold for triple-X detection
###   SEX_Y_USE_MEAN     Use mean Y count as threshold (YES/NO)
###   SEX_Y_MANUAL       Manual Y count threshold if not using mean
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
PLINK1=$(parse_param "PLINK1_PATH" "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# Sex check thresholds
SEX_F_FEMALE_MAX=$(parse_param "SEX_F_FEMALE_MAX"   "0.2")
SEX_F_MALE_MIN=$(parse_param "SEX_F_MALE_MIN"       "0.8")
SEX_F_TURNER=$(parse_param "SEX_F_TURNER"           "0.8")
SEX_F_KLINEFELTER=$(parse_param "SEX_F_KLINEFELTER" "0.8")
SEX_F_XXX=$(parse_param "SEX_F_XXX"                 "-0.15")
SEX_Y_USE_MEAN=$(parse_param "SEX_Y_USE_MEAN"       "YES")
SEX_Y_MANUAL=$(parse_param "SEX_Y_MANUAL"           "0")

SEXCHECK_EXCLUDE_MISMATCH=$(parse_param "SEXCHECK_EXCLUDE_MISMATCH" "YES")
SEXCHECK_EXCLUDE_ANEUPLOIDY=$(parse_param "SEXCHECK_EXCLUDE_ANEUPLOIDY" "YES")
SEXCHECK_AMEND_SEX=$(parse_param "SEXCHECK_AMEND_SEX" "NO")
STEP4_FINAL=$(parse_param "STEP4_FINAL_PREFIX" "")
CONVERTED_PREFIX=$(parse_param "CONVERTED_PREFIX" "")

if [ -n "${STEP4_FINAL}" ] && [ -f "${STEP4_FINAL}.bed" ]; then
    INPUT_PREFIX="${STEP4_FINAL}"
elif [ -n "${CONVERTED_PREFIX}" ] && [ -f "${CONVERTED_PREFIX}.bed" ]; then
    INPUT_PREFIX="${CONVERTED_PREFIX}"
else
    # Fall back to convention -- look for MAF filtered output first
    INPUT_PREFIX=$(ls -t "${OUTPUT_DIR}/step4_snprate/${PROJECT_NAME}"*maf*.bed \
        2>/dev/null | head -1 | sed 's/\.bed$//')
    if [ -z "${INPUT_PREFIX}" ]; then
        INPUT_PREFIX=$(ls -t "${OUTPUT_DIR}/step4_snprate/${PROJECT_NAME}"*.bed \
            2>/dev/null | head -1 | sed 's/\.bed$//')
    fi
    if [ -z "${INPUT_PREFIX}" ]; then
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
echo "        SEXCHECK_EXCLUDE_MISMATCH   : ${SEXCHECK_EXCLUDE_MISMATCH}"
echo "        SEXCHECK_EXCLUDE_ANEUPLOIDY : ${SEXCHECK_EXCLUDE_ANEUPLOIDY}"
echo "        SEXCHECK_AMEND_SEX          : ${SEXCHECK_AMEND_SEX}"
echo "        SEX_F_FEMALE_MAX  : ${SEX_F_FEMALE_MAX}"
echo "        SEX_F_MALE_MIN    : ${SEX_F_MALE_MIN}"
echo "        SEX_F_TURNER      : ${SEX_F_TURNER}"
echo "        SEX_F_KLINEFELTER : ${SEX_F_KLINEFELTER}"
echo "        SEX_F_XXX         : ${SEX_F_XXX}"
echo "        SEX_Y_USE_MEAN    : ${SEX_Y_USE_MEAN}"
echo "        SEX_Y_MANUAL      : ${SEX_Y_MANUAL}"
echo "        INPUT_PREFIX      : ${INPUT_PREFIX}"
echo "        OUTPUT_DIR        : ${OUTPUT_DIR}"
echo "        PROJECT_NAME      : ${PROJECT_NAME}"
echo "        PLOT_ENGINE       : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE MUTUAL EXCLUSION
### SEXCHECK_AMEND_SEX and SEXCHECK_EXCLUDE_MISMATCH cannot both be YES
##########################################################################

if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ] && \
   [ "${SEXCHECK_EXCLUDE_MISMATCH^^}" = "YES" ]; then
    echo "[CRISP] ERROR: SEXCHECK_AMEND_SEX = YES and SEXCHECK_EXCLUDE_MISMATCH = YES"
    echo "[CRISP]        cannot both be active at the same time."
    echo "[CRISP]        Please choose one approach in ${INSTRUCTION_FILE}:"
    echo "[CRISP]          To correct sex in FAM file : SEXCHECK_EXCLUDE_MISMATCH = NO"
    echo "[CRISP]          To exclude mismatches      : SEXCHECK_AMEND_SEX = NO"
    exit 1
fi

##########################################################################
### VALIDATE INPUT
##########################################################################

for ext in .bed .bim .fam; do
    if [ ! -f "${INPUT_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${INPUT_PREFIX}${ext}"
        echo "[CRISP]        Ensure Step 4 (crisp_snprate.sh) has been run first."
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

##########################################################################
### SET UP DIRECTORIES
##########################################################################

STEP5_DIR="${OUTPUT_DIR}/step5_sexcheck"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step5_sexcheck_report.txt"

mkdir -p "${STEP5_DIR}"
mkdir -p "${LOG_DIR}"

##########################################################################
### STEP 5a: PLINK --check-sex
### Generates F-statistic per sample using X chromosome variants
##########################################################################

echo "##########################################################################"
echo "###              STEP 5a: PLINK CHECK-SEX                             ###"
echo "##########################################################################"
echo ""

SEXCHECK_PREFIX="${STEP5_DIR}/${PROJECT_NAME}_sexcheck"
LOG_5A="${LOG_DIR}/step5a_checksex.log"

echo "[CRISP] Running: ${PLINK1} --bfile ${INPUT_PREFIX} --check-sex --out ${SEXCHECK_PREFIX}"
echo ""

${PLINK1} \
    --bfile "${INPUT_PREFIX}" \
    --check-sex \
    --out "${SEXCHECK_PREFIX}" \
    >> "${LOG_5A}" 2>&1

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: PLINK --check-sex failed. Check log: ${LOG_5A}"
    exit 1
fi

echo "[OK]    Sex check file generated: ${SEXCHECK_PREFIX}.sexcheck"
echo ""

##########################################################################
### STEP 5b: Y CHROMOSOME VARIANT COUNT
### Counts non-missing variants on chrY per sample
##########################################################################

echo "##########################################################################"
echo "###              STEP 5b: Y CHROMOSOME COUNT                          ###"
echo "##########################################################################"
echo ""

YCHR_PREFIX="${STEP5_DIR}/${PROJECT_NAME}_ychr"
LOG_5B="${LOG_DIR}/step5b_ychr.log"

echo "[CRISP] Running: ${PLINK1} --bfile ${INPUT_PREFIX} --chr Y --missing --out ${YCHR_PREFIX}"
echo ""

${PLINK1} \
    --bfile "${INPUT_PREFIX}" \
    --chr Y \
    --missing \
    --out "${YCHR_PREFIX}" \
    >> "${LOG_5B}" 2>&1

if [ $? -ne 0 ]; then
    echo "[CRISP] WARNING: Y chromosome count failed -- chrY may not be present."
    echo "[CRISP]          Creating empty Y count file."
    touch "${YCHR_PREFIX}.imiss"
fi

echo "[OK]    Y chromosome count complete."
echo ""

##########################################################################
### STEP 5c: MERGE F-STAT AND Y-COUNT
### Merges .sexcheck and Y chromosome .imiss into single file
### for the R plotting script
##########################################################################

echo "##########################################################################"
echo "###              STEP 5c: MERGE F-STAT AND Y-COUNT                    ###"
echo "##########################################################################"
echo ""

MERGED_FILE="${STEP5_DIR}/${PROJECT_NAME}_f_stat_y_count.txt"

python3 - << PYEOF
import pandas as pd
import os
import sys

sexcheck_file = "${SEXCHECK_PREFIX}.sexcheck"
ychr_file     = "${YCHR_PREFIX}.imiss"
output_file   = "${MERGED_FILE}"

# Load sexcheck
if not os.path.isfile(sexcheck_file):
    print(f"[CRISP] ERROR: .sexcheck file not found: {sexcheck_file}")
    sys.exit(1)

df_sex = pd.read_csv(sexcheck_file, sep=r'\s+')

# Load Y count -- count non-missing variants per sample
if os.path.isfile(ychr_file) and os.path.getsize(ychr_file) > 0:
    df_y = pd.read_csv(ychr_file, sep=r'\s+')
    # YCOUNT = total Y variants minus missing
    df_y['YCOUNT'] = df_y['N_GENO'] - df_y['N_MISS']
    df_y = df_y[['IID', 'YCOUNT']]
    df_merged = df_sex.merge(df_y, on='IID', how='left')
    df_merged['YCOUNT'] = df_merged['YCOUNT'].fillna(0).astype(int)
else:
    print("[CRISP] WARNING: No Y count data found. Setting YCOUNT to 0.")
    df_sex['YCOUNT'] = 0
    df_merged = df_sex

df_merged.to_csv(output_file, sep='\t', index=False)
print(f"[CRISP] Merged file written: {output_file}")
print(f"[CRISP] Samples in merged file: {len(df_merged)}")
PYEOF

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Merge step failed."
    exit 1
fi

echo ""

##########################################################################
### STEP 5d: RUN DETECTION AND PLOTTING
##########################################################################

echo "##########################################################################"
echo "###              STEP 5d: DETECTION AND PLOTTING                      ###"
echo "##########################################################################"
echo ""

LOG_5D="${LOG_DIR}/step5d_plots.log"

if [ "${PLOT_ENGINE^^}" = "R" ]; then
    echo "[CRISP] Running plot_sexcheck.R..."
    ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_sexcheck.R" \
        "${STEP5_DIR}" \
        "${STEP5_DIR}" \
        "${MERGED_FILE}" \
        "${SEX_F_FEMALE_MAX}" \
        "${SEX_F_MALE_MIN}" \
        "${SEX_F_TURNER}" \
        "${SEX_F_KLINEFELTER}" \
        "${SEX_F_XXX}" \
        "${SEX_Y_USE_MEAN}" \
        "${SEX_Y_MANUAL}" \
        >> "${LOG_5D}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: plot_sexcheck.R failed. Check log: ${LOG_5D}"
        exit 1
    fi

elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
    echo "[CRISP] Running plot_sexcheck.py..."
    python3 "${SCRIPT_DIR}/scripts/plot_sexcheck.py" \
        "${STEP5_DIR}" \
        "${STEP5_DIR}" \
        "${MERGED_FILE}" \
        "${SEX_F_FEMALE_MAX}" \
        "${SEX_F_MALE_MIN}" \
        "${SEX_F_TURNER}" \
        "${SEX_F_KLINEFELTER}" \
        "${SEX_F_XXX}" \
        "${SEX_Y_USE_MEAN}" \
        "${SEX_Y_MANUAL}" \
        >> "${LOG_5D}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: plot_sexcheck.py failed. Check log: ${LOG_5D}"
        exit 1
    fi
fi

echo "[OK]    Detection and plotting complete."
echo ""

##########################################################################
### STEP 5e: COUNT EXCLUSIONS
##########################################################################

MISMATCH_FILE="${STEP5_DIR}/id_list_sex_mismatch.txt"
ANEUPLOIDY_FILE="${STEP5_DIR}/id_list_aneuploidies.txt"

N_MISMATCH=0
N_ANEUPLOIDY=0

[ -f "${MISMATCH_FILE}" ]   && N_MISMATCH=$(wc -l < "${MISMATCH_FILE}" | tr -d '[:space:]')
[ -f "${ANEUPLOIDY_FILE}" ] && N_ANEUPLOIDY=$(wc -l < "${ANEUPLOIDY_FILE}" | tr -d '[:space:]')

echo "[CRISP] Sex mismatches flagged : ${N_MISMATCH}"
echo "[CRISP] Aneuploidies flagged   : ${N_ANEUPLOIDY}"
echo ""

##########################################################################
### STEP 5f: APPLY EXCLUSIONS (OPTIONAL)
### Removes sex mismatch and/or aneuploidy samples from dataset.
### Controlled by SEXCHECK_EXCLUDE_MISMATCH and
### SEXCHECK_EXCLUDE_ANEUPLOIDY in crisp_instructions.txt.
##########################################################################

EXCL_COMBINED="${STEP5_DIR}/${PROJECT_NAME}_step5_all_exclusions.txt"
> "${EXCL_COMBINED}"

CLEAN_PREFIX="${INPUT_PREFIX}"
N_EXCLUDED_MISMATCH=0
N_EXCLUDED_ANEUPLOIDY=0

if [ "${SEXCHECK_EXCLUDE_MISMATCH^^}" = "YES" ] && \
   [ -f "${MISMATCH_FILE}" ] && [ -s "${MISMATCH_FILE}" ]; then
    cat "${MISMATCH_FILE}" >> "${EXCL_COMBINED}"
    N_EXCLUDED_MISMATCH="${N_MISMATCH}"
    echo "[CRISP] Sex mismatch samples added to exclusion list: ${N_EXCLUDED_MISMATCH}"
else
    echo "[CRISP] SEXCHECK_EXCLUDE_MISMATCH = NO -- sex mismatch samples retained."
fi

if [ "${SEXCHECK_EXCLUDE_ANEUPLOIDY^^}" = "YES" ] && \
   [ -f "${ANEUPLOIDY_FILE}" ] && [ -s "${ANEUPLOIDY_FILE}" ]; then
    cat "${ANEUPLOIDY_FILE}" >> "${EXCL_COMBINED}"
    N_EXCLUDED_ANEUPLOIDY="${N_ANEUPLOIDY}"
    echo "[CRISP] Aneuploidy samples added to exclusion list  : ${N_EXCLUDED_ANEUPLOIDY}"
else
    echo "[CRISP] SEXCHECK_EXCLUDE_ANEUPLOIDY = NO -- aneuploidy samples retained."
fi

# Remove duplicates from combined exclusion list
if [ -s "${EXCL_COMBINED}" ]; then
    sort -u "${EXCL_COMBINED}" -o "${EXCL_COMBINED}"
    N_TOTAL_EXCL=$(wc -l < "${EXCL_COMBINED}" | tr -d '[:space:]')
    echo ""
    echo "[CRISP] Applying exclusions to dataset..."

    CLEAN_PREFIX="${STEP5_DIR}/${PROJECT_NAME}_step5_clean"
    LOG_EXCL="${LOG_DIR}/step5_exclusions.log"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --remove "${EXCL_COMBINED}" \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_EXCL}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Sample exclusion failed. Check log: ${LOG_EXCL}"
        exit 1
    fi

    SAMPLES_AFTER=$(wc -l < "${CLEAN_PREFIX}.fam" | tr -d '[:space:]')
    echo "[OK]    Exclusions applied."
    echo "[CRISP] Samples before exclusion : ${SAMPLES}"
    echo "[CRISP] Samples excluded         : ${N_TOTAL_EXCL}"
    echo "[CRISP] Samples after exclusion  : ${SAMPLES_AFTER}"
    echo ""
else
    echo ""
    echo "[CRISP] No samples to exclude. Dataset unchanged."
    echo ""
    SAMPLES_AFTER="${SAMPLES}"
fi

##########################################################################
### STEP 5g: SEX AMENDMENT IN FAM FILE (OPTIONAL)
### Updates PEDSEX in .fam file to match genotype-inferred sex.
### Only runs when SEXCHECK_AMEND_SEX = YES.
### Original FAM file is never overwritten.
### A WARNING is printed and logged for each amended sample.
##########################################################################

AMENDMENT_LOG="${OUTPUT_DIR}/${PROJECT_NAME}_sex_amendments.txt"
N_AMENDED=0

if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ]; then

    echo "##########################################################################"
    echo "###              STEP 5g: SEX AMENDMENT IN FAM FILE                   ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP WARNING] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "[CRISP WARNING] SEX AMENDMENT HAS BEEN REQUESTED."
    echo "[CRISP WARNING] The recorded sex in the FAM file will be updated to match"
    echo "[CRISP WARNING] the genotype-inferred sex for flagged samples."
    echo "[CRISP WARNING] Before proceeding, it is strongly recommended to consult"
    echo "[CRISP WARNING] the data provider for each amended sample to confirm"
    echo "[CRISP WARNING] whether the recorded or genotyped sex is correct."
    echo "[CRISP WARNING] Automatic amendment may not always be the right decision."
    echo "[CRISP WARNING] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""

    # Python handles the FAM amendment
    python3 - << PYEOF
import pandas as pd
import os
import sys

merged_file    = "${MERGED_FILE}"
input_fam      = "${INPUT_PREFIX}.fam"
amended_fam    = "${STEP5_DIR}/${PROJECT_NAME}_step5_amended.fam"
amendment_log  = "${AMENDMENT_LOG}"

df_sex  = pd.read_csv(merged_file, sep=r'\s+')
df_fam  = pd.read_csv(input_fam, sep=r'\s+',
                       header=None,
                       names=['FID','IID','PAT','MAT','SEX','PHENO'])

# Identify mismatches where SNPSEX is valid (1 or 2)
mismatches = df_sex[
    (df_sex['PEDSEX'] != df_sex['SNPSEX']) &
    (df_sex['SNPSEX'].isin([1,2]))
][['FID','IID','PEDSEX','SNPSEX']]

if len(mismatches) == 0:
    print("[CRISP] No sex mismatches with valid SNPSEX found. No amendments made.")
    open(amendment_log,'w').close()
    sys.exit(0)

# Write amendment log
with open(amendment_log, 'w') as log:
    log.write("=" * 64 + "\n")
    log.write("  CRISP -- SEX AMENDMENT LOG\n")
    log.write("  Comprehensive Robust Integrated SNP Processing\n")
    log.write("=" * 64 + "\n\n")
    log.write("  WARNING: These samples had their recorded sex updated\n")
    log.write("  to match the genotype-inferred sex (SNPSEX).\n")
    log.write("  Consult the data provider before proceeding.\n\n")
    log.write(f"  {'FID':<12} {'IID':<12} {'Original SEX':<15} {'Amended SEX':<12}\n")
    log.write("-" * 64 + "\n")
    for _, row in mismatches.iterrows():
        orig = {1:'Male',2:'Female',0:'Unknown'}.get(int(row['PEDSEX']),'Unknown')
        new  = {1:'Male',2:'Female'}.get(int(row['SNPSEX']),'Unknown')
        log.write(f"  {row['FID']:<12} {row['IID']:<12} {orig:<15} {new:<12}\n")
        print(f"[CRISP WARNING] Amending: {row['FID']} {row['IID']} "
              f"SEX {int(row['PEDSEX'])} -> {int(row['SNPSEX'])}")
    log.write("\n" + "=" * 64 + "\n")

# Apply amendments to FAM
amendment_map = {
    (str(row['FID']), str(row['IID'])): int(row['SNPSEX'])
    for _, row in mismatches.iterrows()
}

df_fam['SEX'] = df_fam.apply(
    lambda r: amendment_map.get((str(r['FID']), str(r['IID'])), r['SEX']),
    axis=1
)

df_fam.to_csv(amended_fam, sep=' ', header=False, index=False)
print(f"\n[CRISP] Amended FAM written : {amended_fam}")
print(f"[CRISP] Amendment log       : {amendment_log}")
print(f"[CRISP] Samples amended     : {len(mismatches)}")
PYEOF

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Sex amendment step failed."
        exit 1
    fi

    N_AMENDED=$(grep -c "Amending:" <<< "$(cat ${LOG_DIR}/step5_amendment.log 2>/dev/null)" || echo 0)

    echo ""
    echo "[CRISP WARNING] Amendment log written to: ${AMENDMENT_LOG}"
    echo "[CRISP WARNING] Please review before continuing the pipeline."
    echo ""
fi

##########################################################################
### STEP 5h: WRITE STEP 5 REPORT
##########################################################################

{
    echo "=================================================================="
    echo "  CRISP -- STEP 5 SEX CHECK REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Input prefix : ${INPUT_PREFIX}"
    echo "------------------------------------------------------------------"
    echo "  PARAMETERS"
    echo "  SEX_F_FEMALE_MAX            : ${SEX_F_FEMALE_MAX}"
    echo "  SEX_F_MALE_MIN              : ${SEX_F_MALE_MIN}"
    echo "  SEX_F_TURNER                : ${SEX_F_TURNER}"
    echo "  SEX_F_KLINEFELTER           : ${SEX_F_KLINEFELTER}"
    echo "  SEX_F_XXX                   : ${SEX_F_XXX}"
    echo "  SEX_Y_USE_MEAN              : ${SEX_Y_USE_MEAN}"
    if [ "${SEX_Y_USE_MEAN^^}" = "NO" ]; then
    echo "  SEX_Y_MANUAL                : ${SEX_Y_MANUAL}"
    fi
    echo "------------------------------------------------------------------"
    echo "  DETECTION RESULTS"
    echo "  Total samples               : ${SAMPLES}"
    echo "  Sex mismatches flagged      : ${N_MISMATCH}"
    echo "  Aneuploidies flagged        : ${N_ANEUPLOIDY}"
    echo "------------------------------------------------------------------"
    echo "  EXCLUSIONS"
    echo "  SEXCHECK_EXCLUDE_MISMATCH   : ${SEXCHECK_EXCLUDE_MISMATCH}"
    echo "  SEXCHECK_EXCLUDE_ANEUPLOIDY : ${SEXCHECK_EXCLUDE_ANEUPLOIDY}"
    if [ -s "${EXCL_COMBINED}" ]; then
    echo "  Samples excluded            : ${N_TOTAL_EXCL}"
    echo "  Samples after exclusion     : ${SAMPLES_AFTER}"
    echo "  Clean dataset prefix        : ${CLEAN_PREFIX}"
    else
    echo "  Samples excluded            : 0"
    fi
    echo "------------------------------------------------------------------"
    echo "  SEX AMENDMENT"
    echo "  SEXCHECK_AMEND_SEX          : ${SEXCHECK_AMEND_SEX}"
    if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ]; then
    echo ""
    echo "  WARNING: Sex amendment was applied. Please review:"
    echo "  ${AMENDMENT_LOG}"
    echo "  before continuing the pipeline."
    fi
    echo "------------------------------------------------------------------"
    echo "  OUTPUT FILES"
    echo "  Plots       : ${STEP5_DIR}/Sex_check.pdf"
    echo "  Full report : ${STEP5_DIR}/report_sex.mismatch_aneuploidies.txt"
    echo "  Mismatches  : ${MISMATCH_FILE}"
    echo "  Aneuploidies: ${ANEUPLOIDY_FILE}"
    if [ -s "${EXCL_COMBINED}" ]; then
    echo "  Excl. list  : ${EXCL_COMBINED}"
    fi
    if [ "${SEXCHECK_AMEND_SEX^^}" = "YES" ]; then
    echo "  Amendments  : ${AMENDMENT_LOG}"
    fi
    echo "  Log         : ${LOG_DIR}"
    echo "=================================================================="
    echo "  END OF REPORT"
    echo "=================================================================="
} | tee "${REPORT_FILE}"

echo ""
echo "[CRISP] Report written to: ${REPORT_FILE}"
echo ""
echo "[CRISP] Step 5 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
