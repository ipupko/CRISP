#!/bin/bash
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### CHUNK: Relatedness (Step 7)
### Version: 0.2.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Performs pairwise relatedness estimation and removes related samples.
#
# Two methods available via RELATEDNESS_METHOD:
#
#   IBD (default)
#     PLINK 1.9 --genome: pairwise IBD estimation
#     Produces IBD0, IBD1, IBD2 and PI_HAT per pair
#     Best for single-ancestry cohorts up to ~50,000 samples
#     Requires LD pruning first for accurate estimates
#     Well understood, PI_HAT threshold widely used in literature
#
#   KING
#     PLINK 2 --make-king-table: kinship coefficient estimation
#     Robust to population stratification and admixed cohorts
#     Scales to biobank size without LD pruning
#     Recommended for large or multi-ancestry cohorts
#
# Sub-steps:
#   7a: LD pruning (IBD only)
#   7b: Relatedness estimation (IBD or KING)
#   7c: Identify pairs above threshold
#   7d: Select one sample per pair for exclusion (keep higher call rate)
#   7e: Apply exclusions
#   7f: Plotting
#   7g: Write Step 7 report
#
# Exclusion strategy:
#   ALL_RELATED       : remove all individuals in any related pair (default)
#                       simpler and safer for consanguineous cohorts
#                       equivalent to DIVERGE relatedness_removal.R approach
#   HIGHER_CALLRATE   : keep sample with better call rate per pair
#                       more conservative, retains more samples
#   .genome file (IBD) or .kin0 file (KING)
#   related_pairs.txt
#   exclusions_step7_relatedness.txt
#   Relatedness_IBD.pdf or Relatedness_KING.pdf
#   step7_relatedness_report.txt
#
# Usage:
#   bash crisp_relatedness.sh
#   bash crisp_relatedness.sh --config my_project.txt
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
            echo "Usage: bash crisp_relatedness.sh [--config <file>]"
            echo "  --config <file>   Instruction file path (default: crisp_instructions.txt)"
            exit 0
            ;;
        *)
            echo "[CRISP] ERROR: Unknown argument: $1"
            echo "[CRISP]        Usage: bash crisp_relatedness.sh [--config <file>]"
            exit 1
            ;;
    esac
done

##########################################################################
### LOCATE INSTRUCTION FILE
##########################################################################

if [ ! -f "${INSTRUCTION_FILE}" ]; then
    echo "[CRISP] ERROR: Instruction file not found: ${INSTRUCTION_FILE}"
    echo "[CRISP]        Usage: bash crisp_relatedness.sh --config <path>"
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

OUTPUT_DIR=$(parse_param "OUTPUT_DIR"           "./results")
PROJECT_NAME=$(parse_param "PROJECT_NAME"       "project")
KEEP_INTERMEDIATE=$(parse_param "KEEP_INTERMEDIATE" "YES")
RELATEDNESS_METHOD=$(parse_param "RELATEDNESS_METHOD" "IBD")
IBD_PIHAT_CUTOFF=$(parse_param "IBD_PIHAT_CUTOFF"     "0.125")
KING_CUTOFF=$(parse_param "KING_CUTOFF"               "0.0884")
RELATEDNESS_KEEP=$(parse_param "RELATEDNESS_KEEP"     "ALL_RELATED")
PLINK1=$(parse_param "PLINK1_PATH"   "plink")
PLINK2=$(parse_param "PLINK2_PATH"   "plink2")
RSCRIPT=$(parse_param "RSCRIPT_PATH" "Rscript")
PLOT_ENGINE=$(parse_param "PLOT_ENGINE" "R")

# LD pruning parameters for IBD
# Default r2=0.05 is tighter than 0.1 and more appropriate for
# South Asian and other cohorts with stronger LD structure
LD_WINDOW=$(parse_param "LD_WINDOW"     "50")
LD_STEP=$(parse_param "LD_STEP"         "5")
LD_R2=$(parse_param "LD_R2"             "0.05")

# resolve input from Step 6 output
STEP6_CLEAN=$(ls -t "${OUTPUT_DIR}/step6_homozygosity/${PROJECT_NAME}"*clean*.bed \
    2>/dev/null | head -1 | sed 's/\.bed$//')

if [ -n "${STEP6_CLEAN}" ] && [ -f "${STEP6_CLEAN}.bed" ]; then
    INPUT_PREFIX="${STEP6_CLEAN}"
    echo "[CRISP] Using Step 6 clean output as input."
else
    # fall back through earlier steps
    for step_dir in step5_sexcheck step4_snprate step3_callrate; do
        candidate=$(ls -t "${OUTPUT_DIR}/${step_dir}/${PROJECT_NAME}"*clean*.bed \
            2>/dev/null | head -1 | sed 's/\.bed$//')
        if [ -n "${candidate}" ] && [ -f "${candidate}.bed" ]; then
            INPUT_PREFIX="${candidate}"
            echo "[CRISP] Step 6 output not found. Using: ${INPUT_PREFIX}"
            break
        fi
    done
fi

if [ -z "${INPUT_PREFIX:-}" ]; then
    INPUT_PREFIX="${OUTPUT_DIR}/step2_converted/${PROJECT_NAME}_converted_chrclean"
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
echo "        RELATEDNESS_METHOD : ${RELATEDNESS_METHOD}"
if [ "${RELATEDNESS_METHOD^^}" = "IBD" ]; then
echo "        IBD_PIHAT_CUTOFF   : ${IBD_PIHAT_CUTOFF}"
echo "        LD_WINDOW / STEP   : ${LD_WINDOW} / ${LD_STEP} (r2=${LD_R2})"
else
echo "        KING_CUTOFF        : ${KING_CUTOFF}"
fi
echo "        RELATEDNESS_KEEP   : ${RELATEDNESS_KEEP}"
echo "        INPUT_PREFIX       : ${INPUT_PREFIX}"
echo "        PLOT_ENGINE        : ${PLOT_ENGINE}"
echo ""

##########################################################################
### VALIDATE INPUT
##########################################################################

for ext in .bed .bim .fam; do
    if [ ! -f "${INPUT_PREFIX}${ext}" ]; then
        echo "[CRISP] ERROR: Input file not found: ${INPUT_PREFIX}${ext}"
        echo "[CRISP]        Ensure Step 6 (crisp_homozygosity.sh) has been run first."
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

STEP7_DIR="${OUTPUT_DIR}/step7_relatedness"
LOG_DIR="${OUTPUT_DIR}/logs"
REPORT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_step7_relatedness_report.txt"

mkdir -p "${STEP7_DIR}"
mkdir -p "${LOG_DIR}"

# per-run temp folder, avoids /tmp collisions on HPC shared nodes
CRISP_TMP="${OUTPUT_DIR}/.crisp_tmp_$$"
mkdir -p "${CRISP_TMP}"
trap "rm -rf ${CRISP_TMP}" EXIT

##########################################################################
### STEP 7a: LD PRUNING (IBD ONLY)
##########################################################################

PRUNED_PREFIX="${INPUT_PREFIX}"

if [ "${RELATEDNESS_METHOD^^}" = "IBD" ]; then

    echo "##########################################################################"
    echo "###              STEP 7a: LD PRUNING (IBD preparation)               ###"
    echo "##########################################################################"
    echo ""
    echo "[CRISP] Pruning LD: window=${LD_WINDOW}, step=${LD_STEP}, r2=${LD_R2}..."
    echo ""

    PRUNE_PREFIX="${STEP7_DIR}/${PROJECT_NAME}_ldprune"
    LOG_7A="${LOG_DIR}/step7a_ldprune.log"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
        --allow-extra-chr \
        --out "${PRUNE_PREFIX}" \
        >> "${LOG_7A}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: LD pruning failed. Check log: ${LOG_7A}"
        exit 1
    fi

    PRUNED_PREFIX="${STEP7_DIR}/${PROJECT_NAME}_ldpruned"

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --extract "${PRUNE_PREFIX}.prune.in" \
        --allow-extra-chr \
        --make-bed \
        --out "${PRUNED_PREFIX}" \
        >> "${LOG_7A}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Variant extraction after LD pruning failed. Check log: ${LOG_7A}"
        exit 1
    fi

    PRUNED_VARIANTS=$(wc -l < "${PRUNED_PREFIX}.bim" | tr -d '[:space:]')
    echo "[OK]    LD pruning complete."
    echo "[CRISP] Variants before LD pruning : ${VARIANTS}"
    echo "[CRISP] Variants after LD pruning  : ${PRUNED_VARIANTS}"
    echo ""
fi

##########################################################################
### STEP 7b: RELATEDNESS ESTIMATION
##########################################################################

echo "##########################################################################"
echo "###              STEP 7b: RELATEDNESS ESTIMATION                      ###"
echo "##########################################################################"
echo ""

if [ "${RELATEDNESS_METHOD^^}" = "IBD" ]; then

    echo "[CRISP] Running PLINK --genome (IBD)..."
    echo ""

    GENOME_PREFIX="${STEP7_DIR}/${PROJECT_NAME}_ibd"
    LOG_7B="${LOG_DIR}/step7b_ibd.log"

    ${PLINK1} \
        --bfile "${PRUNED_PREFIX}" \
        --genome \
        --allow-extra-chr \
        --min "${IBD_PIHAT_CUTOFF}" \
        --out "${GENOME_PREFIX}" \
        >> "${LOG_7B}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: PLINK --genome failed. Check log: ${LOG_7B}"
        exit 1
    fi

    RELATEDNESS_FILE="${GENOME_PREFIX}.genome"
    echo "[OK]    IBD estimation complete."
    echo "[CRISP] Output: ${RELATEDNESS_FILE}"

elif [ "${RELATEDNESS_METHOD^^}" = "KING" ]; then

    echo "[CRISP] Running PLINK 2 --make-king-table (KING)..."
    echo ""

    KING_PREFIX="${STEP7_DIR}/${PROJECT_NAME}_king"
    LOG_7B="${LOG_DIR}/step7b_king.log"

    ${PLINK2} \
        --bfile "${INPUT_PREFIX}" \
        --make-king-table \
        --king-table-filter "${KING_CUTOFF}" \
        --out "${KING_PREFIX}" \
        >> "${LOG_7B}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: PLINK 2 --make-king-table failed. Check log: ${LOG_7B}"
        exit 1
    fi

    RELATEDNESS_FILE="${KING_PREFIX}.kin0"
    echo "[OK]    KING estimation complete."
    echo "[CRISP] Output: ${RELATEDNESS_FILE}"

fi

echo ""

##########################################################################
### STEP 7c: IDENTIFY RELATED PAIRS
##########################################################################

echo "##########################################################################"
echo "###              STEP 7c: IDENTIFY RELATED PAIRS                      ###"
echo "##########################################################################"
echo ""

PAIRS_FILE="${STEP7_DIR}/${PROJECT_NAME}_related_pairs.txt"

python3 - << PYEOF
import pandas as pd
import sys
import os

method         = "${RELATEDNESS_METHOD}".upper()
relatedness_f  = "${RELATEDNESS_FILE}"
pairs_file     = "${PAIRS_FILE}"
pihat_cutoff   = float("${IBD_PIHAT_CUTOFF}")
king_cutoff    = float("${KING_CUTOFF}")

if not os.path.isfile(relatedness_f):
    print(f"[CRISP] ERROR: Relatedness file not found: {relatedness_f}")
    sys.exit(1)

if method == "IBD":
    df = pd.read_csv(relatedness_f, sep=r'\s+')
    df_pairs = df[df['PI_HAT'] >= pihat_cutoff][
        ['FID1','IID1','FID2','IID2','Z0','Z1','Z2','PI_HAT']
    ].copy()
    df_pairs['METRIC'] = df_pairs['PI_HAT']
    df_pairs['METHOD'] = 'IBD'

    # classify relationship
    def classify_ibd(row):
        if row['PI_HAT'] >= 0.9:              return 'Duplicate/MZ twin'
        elif row['PI_HAT'] >= 0.4:            return 'First degree'
        elif row['PI_HAT'] >= 0.185:          return 'Second degree'
        else:                                  return 'Distant'
    df_pairs['RELATIONSHIP'] = df_pairs.apply(classify_ibd, axis=1)

elif method == "KING":
    df = pd.read_csv(relatedness_f, sep=r'\s+')
    df_pairs = df[df['KINSHIP'] >= king_cutoff][
        ['#FID1','ID1','FID2','ID2','KINSHIP']
    ].copy()
    df_pairs.columns = ['FID1','IID1','FID2','IID2','KINSHIP']
    df_pairs['METRIC'] = df_pairs['KINSHIP']
    df_pairs['METHOD'] = 'KING'

    # classify relationship
    def classify_king(row):
        if row['KINSHIP'] >= 0.354:   return 'Duplicate/MZ twin'
        elif row['KINSHIP'] >= 0.177: return 'First degree'
        elif row['KINSHIP'] >= 0.0884:return 'Second degree'
        else:                          return 'Distant'
    df_pairs['RELATIONSHIP'] = df_pairs.apply(classify_king, axis=1)

df_pairs.to_csv(pairs_file, sep='\t', index=False)
print(f"[CRISP] Related pairs identified : {len(df_pairs):,}")
print(f"[CRISP] Pairs file written       : {pairs_file}")

counts = df_pairs['RELATIONSHIP'].value_counts()
for rel, n in counts.items():
    print(f"         {rel:<25} : {n:,}")
PYEOF

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Pair identification step failed."
    exit 1
fi

echo ""

##########################################################################
### STEP 7d: SELECT SAMPLES FOR EXCLUSION
# Keep sample with higher call rate from each related pair
##########################################################################

echo "##########################################################################"
echo "###              STEP 7d: SELECT SAMPLES FOR EXCLUSION               ###"
echo "##########################################################################"
echo ""

EXCL_FILE="${STEP7_DIR}/${PROJECT_NAME}_exclusions_step7.txt"

python3 - << PYEOF
import pandas as pd
import sys
import os

pairs_file   = "${PAIRS_FILE}"
fam_file     = "${INPUT_PREFIX}.fam"
imiss_prefix = "${STEP7_DIR}/${PROJECT_NAME}_relatedness_imiss"
excl_file    = "${EXCL_FILE}"
keep_method  = "${RELATEDNESS_KEEP}".upper()
method       = "${RELATEDNESS_METHOD}".upper()

if os.path.isfile(pairs_file):
    df_pairs = pd.read_csv(pairs_file, sep='\t')
else:
    print("[CRISP] No related pairs file found. No exclusions needed.")
    open(excl_file, 'w').close()
    sys.exit(0)

if len(df_pairs) == 0:
    print("[CRISP] No related pairs above threshold. No exclusions needed.")
    open(excl_file, 'w').close()
    sys.exit(0)

if keep_method == "ALL_RELATED":
    # remove all individuals appearing in any related pair
    # simpler and safer for consanguineous cohorts with complex relatedness networks
    # equivalent to the DIVERGE relatedness_removal.R approach
    pairs_ids = pd.concat([
        df_pairs[['FID1','IID1']].rename(columns={'FID1':'FID','IID1':'IID'}),
        df_pairs[['FID2','IID2']].rename(columns={'FID2':'FID','IID2':'IID'})
    ]).drop_duplicates()
    to_exclude = set(zip(pairs_ids['FID'].astype(str),
                         pairs_ids['IID'].astype(str)))

elif keep_method == "HIGHER_CALLRATE":
    # keep sample with better call rate from each related pair
    if os.path.isfile(imiss_file):
        df_miss = pd.read_csv(imiss_file, sep=r'\s+')
        miss_map = dict(zip(df_miss['IID'], df_miss['F_MISS']))
    else:
        df_fam  = pd.read_csv(fam_file, sep=r'\s+', header=None,
                               names=['FID','IID','PAT','MAT','SEX','PHENO'])
        miss_map = {iid: 0.0 for iid in df_fam['IID']}

    to_exclude = set()
    for _, row in df_pairs.iterrows():
        iid1, iid2 = str(row['IID1']), str(row['IID2'])
        fid1, fid2 = str(row['FID1']), str(row['FID2'])
        if (fid1,iid1) in to_exclude or (fid2,iid2) in to_exclude:
            continue
        miss1 = miss_map.get(iid1, 0.0)
        miss2 = miss_map.get(iid2, 0.0)
        if miss1 >= miss2:
            to_exclude.add((fid1, iid1))
        else:
            to_exclude.add((fid2, iid2))
else:
    print(f"[CRISP] ERROR: Unknown RELATEDNESS_KEEP '{keep_method}'.")
    print("[CRISP]        Valid options: ALL_RELATED, HIGHER_CALLRATE")
    sys.exit(1)

with open(excl_file, 'w') as f:
    for fid, iid in sorted(to_exclude):
        f.write(f"{fid}\t{iid}\n")

print(f"[CRISP] Samples selected for exclusion : {len(to_exclude):,}")
print(f"[CRISP] Exclusion list written         : {excl_file}")
print(f"[CRISP] Strategy                       : {keep_method}")
PYEOF

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Exclusion selection step failed."
    exit 1
fi

echo ""

##########################################################################
### STEP 7e: APPLY EXCLUSIONS
##########################################################################

echo "##########################################################################"
echo "###              STEP 7e: APPLY EXCLUSIONS                            ###"
echo "##########################################################################"
echo ""

N_EXCL=0
SAMPLES_AFTER="${SAMPLES}"
CLEAN_PREFIX="${INPUT_PREFIX}"

if [ -f "${EXCL_FILE}" ] && [ -s "${EXCL_FILE}" ]; then

    N_EXCL=$(wc -l < "${EXCL_FILE}" | tr -d '[:space:]')
    CLEAN_PREFIX="${STEP7_DIR}/${PROJECT_NAME}_step7_clean"
    LOG_7E="${LOG_DIR}/step7e_exclusions.log"

    echo "[CRISP] Applying ${N_EXCL} exclusions..."

    ${PLINK1} \
        --bfile "${INPUT_PREFIX}" \
        --remove "${EXCL_FILE}" \
        --make-bed \
        --out "${CLEAN_PREFIX}" \
        >> "${LOG_7E}" 2>&1

    if [ $? -ne 0 ]; then
        echo "[CRISP] ERROR: Exclusion step failed. Check log: ${LOG_7E}"
        exit 1
    fi

    SAMPLES_AFTER=$(wc -l < "${CLEAN_PREFIX}.fam" | tr -d '[:space:]')

    echo "[OK]    Exclusions applied."
    echo "[CRISP] Samples before : ${SAMPLES}"
    echo "[CRISP] Samples removed: ${N_EXCL}"
    echo "[CRISP] Samples after  : ${SAMPLES_AFTER}"

else
    echo "[CRISP] No related samples to exclude. Dataset unchanged."
fi

echo ""

##########################################################################
### STEP 7f: PLOTTING
##########################################################################

echo "##########################################################################"
echo "###              STEP 7f: PLOTTING                                    ###"
echo "##########################################################################"
echo ""

LOG_7F="${LOG_DIR}/step7f_plots.log"

if [ "${PLOT_ENGINE^^}" = "R" ]; then
    echo "[CRISP] Running plot_relatedness.R..."
    ${RSCRIPT} "${SCRIPT_DIR}/scripts/plot_relatedness.R" \
        "${STEP7_DIR}" \
        "${RELATEDNESS_FILE}" \
        "${RELATEDNESS_METHOD}" \
        "${IBD_PIHAT_CUTOFF}" \
        "${KING_CUTOFF}" \
        >> "${LOG_7F}" 2>&1

elif [ "${PLOT_ENGINE^^}" = "PYTHON" ]; then
    echo "[CRISP] Running plot_relatedness.py..."
    python3 "${SCRIPT_DIR}/scripts/plot_relatedness.py" \
        "${STEP7_DIR}" \
        "${RELATEDNESS_FILE}" \
        "${RELATEDNESS_METHOD}" \
        "${IBD_PIHAT_CUTOFF}" \
        "${KING_CUTOFF}" \
        >> "${LOG_7F}" 2>&1
fi

if [ $? -ne 0 ]; then
    echo "[CRISP] ERROR: Plotting failed. Check log: ${LOG_7F}"
    exit 1
fi

echo "[OK]    Plots saved to: ${STEP7_DIR}"
echo ""

##########################################################################
### STEP 7g: WRITE REPORT
##########################################################################

N_PAIRS=$(wc -l < "${PAIRS_FILE}" 2>/dev/null | tr -d '[:space:]') || N_PAIRS=0
N_PAIRS=$((N_PAIRS - 1))  # subtract header
[ "${N_PAIRS}" -lt 0 ] && N_PAIRS=0

{
    echo "=================================================================="
    echo "  CRISP: STEP 7 RELATEDNESS REPORT"
    echo "  Comprehensive Robust Integrated SNP Processing"
    echo "=================================================================="
    echo "  Project      : ${PROJECT_NAME}"
    echo "  Date         : $(date)"
    echo "  Method       : ${RELATEDNESS_METHOD}"
    if [ "${RELATEDNESS_METHOD^^}" = "IBD" ]; then
    echo "  PI_HAT cutoff: ${IBD_PIHAT_CUTOFF}"
    echo "  LD pruning   : window=${LD_WINDOW}, step=${LD_STEP}, r2=${LD_R2}"
    else
    echo "  KING cutoff  : ${KING_CUTOFF}"
    fi
    echo "------------------------------------------------------------------"
    echo "  RELATEDNESS RESULTS"
    echo "  Related pairs found    : ${N_PAIRS}"
    echo "  Samples excluded       : ${N_EXCL}"
    echo "  Strategy               : ${RELATEDNESS_KEEP}"
    echo "------------------------------------------------------------------"
    echo "  SAMPLE COUNTS"
    echo "  Samples before         : ${SAMPLES}"
    echo "  Samples after          : ${SAMPLES_AFTER}"
    echo "------------------------------------------------------------------"
    echo "  OUTPUT FILES"
    echo "  Relatedness file : ${RELATEDNESS_FILE}"
    echo "  Related pairs    : ${PAIRS_FILE}"
    echo "  Exclusion list   : ${EXCL_FILE}"
    if [ "${RELATEDNESS_METHOD^^}" = "IBD" ]; then
    echo "  Plot             : ${STEP7_DIR}/Relatedness_IBD.pdf"
    else
    echo "  Plot             : ${STEP7_DIR}/Relatedness_KING.pdf"
    fi
    echo "  Log              : ${LOG_DIR}"
    echo "=================================================================="
    echo "  END OF REPORT"
    echo "=================================================================="
} | tee "${REPORT_FILE}"

echo ""
echo "[CRISP] Report written to: ${REPORT_FILE}"
echo ""
echo "[CRISP] Step 7 complete."
echo ""
echo "[CRISP] ${LAUNCH_MSG}"
echo ""
