#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Hardy-Weinberg Equilibrium Filtering (Step 8)
### Version: 0.1.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Performs Hardy-Weinberg Equilibrium filtering on genotype data.
# Variants with significant deviation from HWE are flagged and
# optionally removed.
#
# HWE deviation can indicate:
#   Genotyping error       : excess heterozygotes or homozygotes
#   Population structure   : deviation driven by ancestry mixing
#   True association signal: deviation in cases only (not controls)
#
# To avoid removing true association signals, CRISP supports
# control-only HWE filtering when phenotype data is available.
#
# Sub-steps:
#   8a: Compute HWE statistics via PLINK --hardy
#   8b: Apply HWE filter via PLINK --hwe
#   8c: Generate report and plots
#
# Outputs:
#   .hwe file                            HWE statistics per variant
#   HWE_distribution.pdf                 P-value distribution histogram
#   step8_hwe_report.txt                 Step 8 plain-text report
#   excluded_hwe_variants.txt            Variants excluded by HWE
#
# Parameters (all overridable in crisp_instructions.txt):
#   HWE_THRESHOLD  : p-value cutoff (default 1e-6)
#   HWE_MODE       : ALL | CONTROLS_ONLY
#   HWE_PHENO_FILE : phenotype file for CONTROLS_ONLY mode
#
# Usage:
#   bash crisp_hwe.sh
#   bash crisp_hwe.sh --config my_project.txt
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
            echo "Usage: bash crisp_hwe.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_hwe.sh [--config <file>]"
            exit 1
            ;;
    esac
done

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_hwe.sh --config <path>"
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

OUTPUT_DIR=$(parse_param "OUTPUT_DIR"         "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME"     "project")
KEEP_INTERMEDIATE=$(parse_param "KEEP_INTERMEDIATE" "YES")
HWE_THRESHOLD=$(parse_param "HWE_THRESHOLD"   "0.000001")
HWE_MODE=$(parse_param "HWE_MODE"             "ALL")
HWE_PHENO_FILE=$(parse_param "HWE_PHENO_FILE" "")
HWE_META=$(parse_param "HWE_META"             "NO")
HWE_STRATIFY=$(parse_param "HWE_STRATIFY"     "NO")
PLINK1=$(parse_param "PLINK1_PATH"            "plink")
RSCRIPT=$(parse_param "RSCRIPT_PATH"          "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE"       "R")

# resolve input from Step 7 clean output, fall back through earlier steps
INPUT_PREFIX=""
for step_dir in step7_relatedness step6_homozygosity \
                step5_sexcheck step4_snprate step3_callrate; do
    candidate=$(ls -t "${OUTPUT_DIR}/${step_dir}/${PROJECT_NAME}"*clean*.bed \
        2>/dev/null | head -1 | sed 's/\.bed$//')
    if [ -n "${candidate}" ] && [ -f "${candidate}.bed" ]; then
        INPUT_PREFIX="${candidate}"
        echo "[CRISP] Using ${step_dir} clean output as input."
        break
    fi
done

if [ -z "${INPUT_PREFIX}" ]; then
    INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
fi

##########################################################################
### VALIDATE CONTROLS-ONLY MODE
##########################################################################

if [ "${HWE_MODE^^}" = "CONTROLS_ONLY" ]; then
    if [ -z "${HWE_PHENO_FILE}" ]; then
        echo "[CRISP] ERROR: HWE_MODE = CONTROLS_ONLY requires HWE_PHENO_FILE to be set."
        echo "[CRISP]        Add HWE_PHENO_FILE = /path/to/pheno.txt to ${INSTRUCTION_FILE}"
        exit 1
    fi
    if [ ! -f "${HWE_PHENO_FILE}" ]; then
        echo "[CRISP] ERROR: HWE_PHENO_FILE not found: ${HWE_PHENO_FILE}"
        exit 1
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
echo "        HWE_THRESHOLD : ${HWE_THRESHOLD}"
echo "        HWE_MODE      : ${HWE_MODE}"
echo "        HWE_META      : ${HWE_META}"
echo "        HWE_STRATIFY  : ${HWE_STRATIFY}"
if [ "${HWE_MODE^^}" = "CONTROLS_ONLY" ] || [ "${HWE_STRATIFY^^}" = "YES" ]; then
echo "        HWE_PHENO_FILE: ${HWE_PHENO_FILE}"
fi
echo "        INPUT_PREFIX  : ${INPUT_PREFIX}"
echo "        PLOT_ENGINE   : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE INPUT
##########################################################################

for ext in .bed .bim .fam; do
    if [ ! -f "${INPUT_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${INPUT_PREFIX}${ext}"
        echo "[CRISP]        Ensure Step 7 (crisp_relatedness.sh) has been run first."
        exit 1
    fi
done

echo "[OK]    Input files verified."
echo ""

SAMPLES=$(wc -l < "${INPUT_PREFIX}.fam"  | tr -d '[:space:]')
VARIANTS=$(wc -l < "${INPUT_PREFIX}.bim" | tr -d '[:space:]')

echo "[CRISP] Samples  : ${SAMPLES}"
echo "[CRISP] Variants : ${VARIANTS}"
echo ""

##########################################################################
### SET UP DIRECTORIES
##########################################################################

STEP8_DIR="${OUTPUT_DIR}/step8_hwe"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step8_hwe_report.txt"

mkdir -p "${STEP8_DIR}"
mkdir -p "${LOG_DIR}"

CRISP_TMP="${OUTPUT_DIR}/.crisp_tmp_$$"
mkdir -p "${CRISP_TMP}"
trap "rm -rf ${CRISP_TMP}" EXIT

##########################################################################
### STEP 8a: COMPUTE HWE STATISTICS
##########################################################################

echo "##########################################################################"
echo "###              STEP 8a: COMPUTE HWE STATISTICS                      ###"
echo "##########################################################################"
echo ""

HWE_STATS_PREFIX="${STEP8_DIR}/${PROJECT_NAME}_hwe_stats"
LOG_8A="${LOG_DIR}/step8a_hardy.log"

if [ "${HWE_MODE^^}" = "CONTROLS_ONLY" ] && [ -n "${HWE_PHENO_FILE}" ]; then
    echo "[CRISP] Computing HWE statistics (controls only)..."
    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --pheno "${HWE_PHENO_FILE}" \
        --hardy \
        --out "${HWE_STATS_PREFIX}" \
        >> "${LOG_8A}" 2>&1
else
    echo "[CRISP] Computing HWE statistics (all samples)..."
    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --hardy \
        --out "${HWE_STATS_PREFIX}" \
        >> "${LOG_8A}" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: PLINK --hardy failed. Check log: ${LOG_8A}"
    exit 1
fi

echo "[OK]    HWE statistics computed: ${HWE_STATS_PREFIX}.hwe"
echo ""

##########################################################################
### STEP 8b: APPLY HWE FILTER
##########################################################################

echo "##########################################################################"
echo "###              STEP 8b: APPLY HWE FILTER                            ###"
echo "##########################################################################"
echo ""

CLEAN_PREFIX="${STEP8_DIR}/${PROJECT_NAME}_step8_clean"
LOG_8B="${LOG_DIR}/step8b_hwe_filter.log"

if [ "${HWE_MODE^^}" = "CONTROLS_ONLY" ] && [ -n "${HWE_PHENO_FILE}" ]; then
    echo "[CRISP] Applying HWE filter (controls only, p < ${HWE_THRESHOLD})..."
    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --pheno "${HWE_PHENO_FILE}" \
        --hwe "${HWE_THRESHOLD}" include-nonctrl \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_8B}" 2>&1
else
    echo "[CRISP] Applying HWE filter (all samples, p < ${HWE_THRESHOLD})..."
    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --hwe "${HWE_THRESHOLD}" \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_8B}" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: PLINK --hwe filter failed. Check log: ${LOG_8B}"
    exit 1
fi

echo "[OK]    HWE filter applied."
echo ""

# count excluded variants
VARIANTS_AFTER=$(wc -l < "${CLEAN_PREFIX}.bim" | tr -d '[:space:]')
N_EXCLUDED=$((VARIANTS - VARIANTS_AFTER))

echo "[CRISP] Variants before : ${VARIANTS}"
echo "[CRISP] Variants excluded: ${N_EXCLUDED}"
echo "[CRISP] Variants after  : ${VARIANTS_AFTER}"
echo ""

# write exclusion list
EXCL_FILE="${STEP8_DIR}/${PROJECT_NAME}_excluded_hwe_variants.txt"
python3 - << PYEOF
import pandas as pd
import os
import sys

hwe_file    = "${HWE_STATS_PREFIX}.hwe"
bim_before  = "${INPUT_PREFIX}.bim"
bim_after   = "${CLEAN_PREFIX}.bim"
excl_file   = "${EXCL_FILE}"
threshold   = float("${HWE_THRESHOLD}")

if not os.path.isfile(hwe_file):
    print("[CRISP] WARNING: HWE stats file not found. Skipping exclusion list.")
    sys.exit(0)

df_hwe = pd.read_csv(hwe_file, sep=r'\s+')

# identify variants in bim before but not in bim after
bim_cols = ['CHR','SNP','CM','POS','A1','A2']
df_before = pd.read_csv(bim_before, sep=r'\s+', header=None, names=bim_cols)
df_after  = pd.read_csv(bim_after,  sep=r'\s+', header=None, names=bim_cols)
excl_snps = set(df_before['SNP']) - set(df_after['SNP'])

df_excl = df_hwe[df_hwe['SNP'].isin(excl_snps)][['CHR','SNP','P']].copy()
df_excl.to_csv(excl_file, sep='\t', index=False)
print(f"[CRISP] Exclusion list written: {excl_file}")
print(f"[CRISP] Variants in exclusion list: {len(df_excl):,}")
PYEOF

echo ""

##########################################################################
### STEP 8c: STRATIFIED HWE AND META-ANALYSIS (OPTIONAL)
# Runs when HWE_META = YES
# Stratifies HWE p-values by case/control and/or ancestry group
# Combines across strata using Fisher's method
# Classifies each variant as: error, signal, stratification, or pass
##########################################################################

META_RESULTS=""

if [ "${HWE_META^^}" = "YES" ]; then

    echo "##########################################################################"
    echo "###              STEP 8c: STRATIFIED HWE AND META-ANALYSIS            ###"
    echo "##########################################################################"
    echo ""

    if [ -z "${HWE_PHENO_FILE}" ] || [ ! -f "${HWE_PHENO_FILE}" ]; then
        echo "[CRISP] WARNING: HWE_META = YES but HWE_PHENO_FILE not found."
        echo "[CRISP]          Stratification by phenotype skipped."
        echo "[CRISP]          Running meta-analysis on all samples only."
    fi

    META_RESULTS="${STEP8_DIR}/${PROJECT_NAME}_hwe_meta.txt"
    LOG_8C="${LOG_DIR}/step8c_meta.log"

    python3 - << PYEOF >> "${LOG_8C}" 2>&1
import pandas as pd
import numpy as np
from scipy import stats
import os
import sys

hwe_file    = "${HWE_STATS_PREFIX}.hwe"
pheno_file  = "${HWE_PHENO_FILE}"
meta_out    = "${META_RESULTS}"
threshold   = float("${HWE_THRESHOLD}")

if not os.path.isfile(hwe_file):
    print(f"[META] ERROR: HWE file not found: {hwe_file}")
    sys.exit(1)

df = pd.read_csv(hwe_file, sep=r'\s+')

# split into ALL, AFF (cases), UNAFF (controls) strata if present
strata = {}
for test in df['TEST'].unique():
    subset = df[df['TEST'] == test][['CHR','SNP','P']].copy()
    subset = subset[subset['P'].notna() & (subset['P'] > 0)]
    subset = subset.rename(columns={'P': f'P_{test}'})
    strata[test] = subset

# start from ALL stratum
if 'ALL' not in strata:
    print("[META] ERROR: No ALL rows in HWE file.")
    sys.exit(1)

df_meta = strata['ALL'].copy()

# merge in case/control strata if available
for test in ['AFF', 'UNAFF']:
    if test in strata:
        df_meta = df_meta.merge(strata[test], on=['CHR','SNP'], how='left')

# Fisher's combined p-value across available strata
p_cols = [c for c in df_meta.columns if c.startswith('P_')]
print(f"[META] Strata available: {p_cols}")

def fishers_combined(row, cols):
    ps = [row[c] for c in cols if pd.notna(row[c]) and row[c] > 0]
    if len(ps) == 0:
        return np.nan
    chi2 = -2 * sum(np.log(p) for p in ps)
    return stats.chi2.sf(chi2, df=2*len(ps))

df_meta['P_META'] = df_meta.apply(
    lambda r: fishers_combined(r, p_cols), axis=1
)

# classify each variant
def classify(row, threshold):
    p_all   = row.get('P_ALL',   1.0) or 1.0
    p_aff   = row.get('P_AFF',   1.0) or 1.0
    p_unaff = row.get('P_UNAFF', 1.0) or 1.0
    p_meta  = row.get('P_META',  1.0) or 1.0
    has_strata = 'P_AFF' in row.index and 'P_UNAFF' in row.index

    if p_meta >= threshold:
        return 'PASS'
    if not has_strata:
        return 'FAIL_ALL'
    if p_unaff < threshold and p_aff < threshold:
        return 'FAIL_ALL'          # likely genotyping error, remove
    if p_aff < threshold and p_unaff >= threshold:
        return 'FAIL_CASES_ONLY'   # potential true signal, flag and keep
    if p_unaff < threshold and p_aff >= threshold:
        return 'FAIL_CONTROLS_ONLY' # remove
    return 'FAIL_ALL'

df_meta['CLASS'] = df_meta.apply(
    lambda r: classify(r, threshold), axis=1
)

counts = df_meta['CLASS'].value_counts()
print(f"\n[META] Classification results:")
for cls, n in counts.items():
    print(f"       {cls:<25}: {n:,}")

df_meta.to_csv(meta_out, sep='\t', index=False)
print(f"\n[META] Meta-analysis results written: {meta_out}")
PYEOF

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: HWE meta-analysis failed. Check log: ${LOG_8C}"
        exit 1
    fi

    echo "[OK]    HWE meta-analysis complete."
    echo "[CRISP] Results: ${META_RESULTS}"
    echo ""

fi

##########################################################################
### STEP 8d: PLOTTING
##########################################################################

echo "##########################################################################"
echo "###              STEP 8c: PLOTTING                                    ###"
echo "##########################################################################"
echo ""

LOG_8C="${LOG_DIR}/step8d_plots.log"

if [ "${PLOT_ENGINE^^}" = "R" ]; then
    echo "[CRISP] Running plot_hwe.R..."
    ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_hwe.R" \
        "${STEP8_DIR}" \
        "${HWE_STATS_PREFIX}.hwe" \
        "${HWE_THRESHOLD}" \
        "${HWE_MODE}" \
        "${HWE_META}" \
        "${META_RESULTS}" \
        >> "${LOG_8C}" 2>&1

elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
    echo "[CRISP] Running plot_hwe.py..."
    python3 "${SCRIPT_DIR}/scripts/plot_hwe.py" \
        "${STEP8_DIR}" \
        "${HWE_STATS_PREFIX}.hwe" \
        "${HWE_THRESHOLD}" \
        "${HWE_MODE}" \
        "${HWE_META}" \
        "${META_RESULTS}" \
        >> "${LOG_8C}" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Plotting failed. Check log: ${LOG_8C}"
    exit 1
fi

echo "[OK]    Plots saved to: ${STEP8_DIR}"
echo ""

##########################################################################
### WRITE STEP 8 REPORT
##########################################################################

{
    echo "=================================================================="
    echo "  CRISP: STEP 8 HWE FILTERING REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Input prefix : ${INPUT_PREFIX}"
    echo "------------------------------------------------------------------"
    echo "  PARAMETERS"
    echo "  HWE_THRESHOLD : ${HWE_THRESHOLD}"
    echo "  HWE_MODE      : ${HWE_MODE}"
    echo "  HWE_META      : ${HWE_META}"
    if [ "${HWE_MODE^^}" = "CONTROLS_ONLY" ] || [ "${HWE_STRATIFY^^}" = "YES" ]; then
    echo "  HWE_PHENO_FILE: ${HWE_PHENO_FILE}"
    fi
    echo "------------------------------------------------------------------"
    echo "  VARIANT COUNTS"
    echo "  Variants before HWE filter : ${VARIANTS}"
    echo "  Variants excluded           : ${N_EXCLUDED}"
    echo "  Variants after HWE filter  : ${VARIANTS_AFTER}"
    echo "  Sample count (unchanged)   : ${SAMPLES}"
    if [ "${HWE_META^^}" = "YES" ] && [ -f "${META_RESULTS}" ]; then
    echo "------------------------------------------------------------------"
    echo "  META-ANALYSIS CLASSIFICATION"
    python3 -c "
import pandas as pd
df = pd.read_csv('${META_RESULTS}', sep='\t')
counts = df['CLASS'].value_counts()
labels = {
    'PASS'             : 'Pass HWE meta-analysis',
    'FAIL_ALL'         : 'Fail all strata (likely error)',
    'FAIL_CASES_ONLY'  : 'Fail cases only (potential signal)',
    'FAIL_CONTROLS_ONLY': 'Fail controls only',
}
for cls, n in counts.items():
    lbl = labels.get(cls, cls)
    print(f'  {lbl:<40} : {n:,}')
"
    fi
    echo "------------------------------------------------------------------"
    echo "  OUTPUT FILES"
    echo "  HWE stats     : ${HWE_STATS_PREFIX}.hwe"
    echo "  Clean dataset : ${CLEAN_PREFIX}"
    echo "  Excluded list : ${EXCL_FILE}"
    if [ "${HWE_META^^}" = "YES" ] && [ -n "${META_RESULTS}" ]; then
    echo "  Meta-analysis : ${META_RESULTS}"
    fi
    echo "  Plot          : ${STEP8_DIR}/HWE_distribution.pdf"
    echo "  QQ plot       : ${STEP8_DIR}/HWE_qq.pdf"
    echo "  Log           : ${LOG_DIR}"
    echo "=================================================================="
    echo "  END OF REPORT"
    echo "=================================================================="
} | tee "${REPORT_FILE}"

echo ""
echo "[CRISP] Report written to: ${REPORT_FILE}"
echo ""
echo "[CRISP] Step 8 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
