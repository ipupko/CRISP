#!/usr/bin/env bash
# ##############################################################
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
# ##############################################################
# Script : crisp_pca_pass1.sh
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.4.0
# Purpose:
#   Step 6 — PCA Pass 1 (ancestry assignment, pre-cleaning).
#   This is a lightweight first-pass PCA run on the dataset
#   immediately after sex check (Step 5). The goal is to assign
#   broad ancestry labels to every sample so that Steps 7
#   (heterozygosity), 8 (relatedness), and 9 (HWE) can all
#   apply population-stratified thresholds — the key insight
#   being that per-ancestry thresholds are far more meaningful
#   than cohort-wide ones, especially for the Wahlund effect
#   in HWE testing.
#
#   This is NOT the final PCA. PCA Pass 2 (Step 10) runs on
#   the fully cleaned dataset and produces the definitive PC
#   covariates for GWAS. Pass 1 is ancestry-assignment only.
#
# Sub-steps:
#   6a — LD pruning           (PLINK 1.9, --indep-pairwise 50 5 0.05)
#   6b — PCA computation      (PLINK 2 primary, PLINK 1.9 fallback)
#   6c — MDS computation      (PLINK 1.9, --mds-plot, cross-validation)
#   6d — GMM clustering       (Python, sklearn) + SD check within clusters
#   6e — PCA/MDS cross-val    (cluster agreement rate)
#   6f — Ancestry assignment  (satellite file written)
#   6g — Plotting             (Plot 1 monochrome + Plot 2 coloured, both PCA & MDS)
#   6h — Report + JSON export
#
# Outputs:
#   step6_pca_pass1/
#     {PROJECT_NAME}_pruned.bed/.bim/.fam         LD-pruned dataset
#     {PROJECT_NAME}_pca.eigenvec                  PC scores (all samples)
#     {PROJECT_NAME}_pca.eigenval                  Variance explained
#     {PROJECT_NAME}_mds.mds                       MDS dimensions
#     {PROJECT_NAME}_ancestry.txt                  Satellite ancestry file
#     {PROJECT_NAME}_cluster_{N}.txt               Per-cluster sample ID lists
#     {PROJECT_NAME}_step6_pca_pass1_report.txt    Plain-text report
#     {PROJECT_NAME}_step6_pca_pass1.json          Machine-readable summary
#   plots/
#     {PROJECT_NAME}_pca_mono.pdf                  PCA Plot 1 — monochrome
#     {PROJECT_NAME}_pca_colour.pdf                PCA Plot 2 — coloured
#     {PROJECT_NAME}_mds_mono.pdf                  MDS Plot 1 — monochrome
#     {PROJECT_NAME}_mds_colour.pdf                MDS Plot 2 — coloured
#
# Satellite file columns:
#   FID  IID  SUPERPOP  SUBPOP  CONFIDENCE  SOURCE  OUTLIER
#   MEMBERSHIP_P  MDS_AGREE  MAF_POP_TAG  PC1..PCn  ASSIGN_PASS
#
# Parameters (all overridable in crisp_instructions.txt):
#   ASSIGN_POP           = NO    — write pop labels to satellite file
#   ASSIGN_SUB_POP       = NO    — detect sub-clusters within groups
#   RUN_1KG_PROJECTION   = NO    — merge with 1KG for named superpops
#   KG_PLINK_PREFIX      =       — path to 1KG PLINK prefix (if above YES)
#   COHORT_ADMIXED       = NO    — draw confidence ellipse per cluster
#   PCA_N_CLUSTERS       = AUTO  — number of GMM clusters (AUTO or integer)
#   PCA_OUTLIER_SD       = 6     — SD threshold for outlier flagging
#   PCA_ELLIPSE_SD       = 3     — SD radius for confidence ellipse
#   PCA_N_PCS            = 20    — number of PCs to compute
#   PCA_RUN_MDS_CROSSVAL = NO    — run MDS alongside PCA for cross-validation
#   MDS_N_DIMS           = 10    — MDS dimensions to compute
#
# FUTURE:
#   — xgen native PCA without PLINK decompression (Project Mimir)
#   — flashpca2 / smartpca / EIGENSOFT as alternative engines
#   — PLINK 2 --pca approx (randomised SVD) for biobank-scale N
#
# Usage:
#   bash crisp_pca_pass1.sh
#   bash crisp_pca_pass1.sh --config my_project.txt
#
##########################################################################

set -euo pipefail

##########################################################################
#  _____      _ _   _       _ _
# |_   _|    (_) | (_)     | (_)
#   | | _ __  _| |_ _  __ _| |_ ___  ___
#   | || '_ \| | __| |/ _` | | / __|/ _ \
#  _| || | | | | |_| | (_| | | \__ \  __/
#  \___/_| |_|_|\__|_|\__,_|_|_|___/\___|
##########################################################################
# WHY: Source _crisp_flavour.sh before anything else. It exports palette
# variables, logging functions, and _init_runtime which validates the
# environment. Without this the script has no colour-mode awareness
# and no standardised logging.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/_crisp_flavour.sh"
_init_runtime

##########################################################################
# ______
# | ___ \
# | |_/ /_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___
# |  __/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | | | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# \_|  \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
##########################################################################
# WHY: Parse --config flag first so all subsequent reads use the correct
# instruction file. Default is crisp_instructions.txt in the working dir.

INSTRUCTION_FILE="crisp_instructions.txt"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            INSTRUCTION_FILE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: bash crisp_pca_pass1.sh [--config <file>]"
            exit 0
            ;;
        *)
            crisp_error "Unknown argument: $1"
            exit 1
            ;;
    esac
done

##########################################################################
# ______                              _
# | ___ \                            | |
# | |_/ /_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___
# |  __/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | | | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# \_|  \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
##########################################################################
# WHY: Read every parameter from the instruction file with safe defaults.
# All parameters have sensible defaults so the script runs out-of-the-box
# on a minimal instruction file, while allowing full customisation.

crisp_log "Reading instruction file: ${INSTRUCTION_FILE}"

_param() {
    # Read a key = value pair from the instruction file.
    # Returns the default if the key is absent or commented out.
    local key="$1" default="$2"
    local val
    val=$(grep -E "^\s*${key}\s*=" "${INSTRUCTION_FILE}" 2>/dev/null \
          | tail -1 | sed 's/.*=\s*//' | sed 's/\s*#.*//' | tr -d '[:space:]')
    echo "${val:-${default}}"
}

# Core identifiers
PROJECT_NAME=$(_param "PROJECT_NAME"      "crisp_project")
INPUT_PREFIX=$(_param "INPUT_PREFIX"      "")
STEP6_INPUT=$(_param  "STEP6_INPUT"       "")   # override auto-detect if needed

# Directories
WORK_DIR=$(_param     "WORK_DIR"          ".")
STEP6_DIR="${WORK_DIR}/step6_pca_pass1"
PLOT_DIR="${STEP6_DIR}/plots"
LOG_DIR="${STEP6_DIR}/logs"

# Ancestry assignment flags
ASSIGN_POP=$(_param           "ASSIGN_POP"           "NO")
ASSIGN_SUB_POP=$(_param       "ASSIGN_SUB_POP"       "NO")
RUN_1KG_PROJECTION=$(_param   "RUN_1KG_PROJECTION"   "NO")
KG_PLINK_PREFIX=$(_param      "KG_PLINK_PREFIX"      "")
COHORT_ADMIXED=$(_param       "COHORT_ADMIXED"       "NO")

# PCA parameters
PCA_N_PCS=$(_param            "PCA_N_PCS"            "20")
PCA_LD_WINDOW=$(_param        "PCA_LD_WINDOW"        "50")
PCA_LD_STEP=$(_param          "PCA_LD_STEP"          "5")
PCA_LD_R2=$(_param            "PCA_LD_R2"            "0.05")
PCA_N_CLUSTERS=$(_param       "PCA_N_CLUSTERS"       "AUTO")
PCA_OUTLIER_SD=$(_param       "PCA_OUTLIER_SD"       "6")
PCA_ELLIPSE_SD=$(_param       "PCA_ELLIPSE_SD"       "3")

# MDS cross-validation
PCA_RUN_MDS_CROSSVAL=$(_param "PCA_RUN_MDS_CROSSVAL" "NO")
MDS_N_DIMS=$(_param           "MDS_N_DIMS"           "10")

# 1KG reference panel
KG_REF_POP=$(_param           "KG_REF_POP"           "AUTO")

# Step 4 MAF population tag (for satellite file concordance check)
MONO_REF_POP=$(_param         "MONO_REF_POP"         "")

# Plotting
PLOT_COLOUR_MODE=$(_param     "PLOT_COLOUR_MODE"     "STANDARD")
PLOT_BACKGROUND=$(_param      "PLOT_BACKGROUND"      "LIGHT")
PLOT_ENGINE=$(_param          "PLOT_ENGINE"          "BOTH")   # R | PYTHON | BOTH

# Misc
REFERENCE_ASSEMBLY=$(_param   "REFERENCE_ASSEMBLY"   "GRCh38")

##########################################################################
#  _   _       _ _     _       _
# | | | |     | (_)   | |     | |
# | | | | __ _| |_  __| | __ _| |_ ___
# | | | |/ _` | | |/ _` |/ _` | __/ _ \
# \ \_/ / (_| | | | (_| | (_| | ||  __/
#  \___/ \__,_|_|_|\__,_|\__,_|\__\___|
##########################################################################
# WHY: Fail fast with clear messages rather than letting PLINK or Python
# produce cryptic errors downstream. We check: instruction file exists,
# input PLINK dataset exists, PLINK is available, Python + sklearn
# are available for GMM clustering.

crisp_log "Validating environment..."

# Instruction file
[[ -f "${INSTRUCTION_FILE}" ]] \
    || crisp_die "Instruction file not found: ${INSTRUCTION_FILE}"

# Auto-detect input dataset from Step 5 output if not explicitly set
if [[ -z "${STEP6_INPUT}" ]]; then
    STEP5_DIR="${WORK_DIR}/step5_sexcheck"
    # Step 5 produces the clean dataset excluding sex-check failures
    CANDIDATE="${STEP5_DIR}/${PROJECT_NAME}_sexcheck_clean"
    if [[ -f "${CANDIDATE}.bed" ]]; then
        STEP6_INPUT="${CANDIDATE}"
        crisp_log "Auto-detected Step 5 output: ${STEP6_INPUT}"
    else
        # Fallback: use INPUT_PREFIX directly
        STEP6_INPUT="${WORK_DIR}/${INPUT_PREFIX}"
        crisp_warn "Step 5 output not found — falling back to INPUT_PREFIX: ${STEP6_INPUT}"
    fi
fi

[[ -f "${STEP6_INPUT}.bed" ]] \
    || crisp_die "Input .bed not found: ${STEP6_INPUT}.bed"
[[ -f "${STEP6_INPUT}.bim" ]] \
    || crisp_die "Input .bim not found: ${STEP6_INPUT}.bim"
[[ -f "${STEP6_INPUT}.fam" ]] \
    || crisp_die "Input .fam not found: ${STEP6_INPUT}.fam"

# PLINK availability
# WHY: We prefer PLINK 1.9 for LD pruning and MDS (--mds-plot not in PLINK 2).
# We prefer PLINK 2 for PCA (faster, more modern). We detect both and route
# accordingly, with clear warnings when falling back.
PLINK1=""
PLINK2=""
for cmd in plink plink1 plink19; do
    if command -v "${cmd}" &>/dev/null; then
        PLINK1="${cmd}"
        break
    fi
done
if command -v plink2 &>/dev/null; then
    PLINK2="plink2"
fi

[[ -n "${PLINK1}" ]] \
    || crisp_die "PLINK 1.9 not found. Required for LD pruning and MDS."

if [[ -z "${PLINK2}" ]]; then
    crisp_warn "PLINK 2 not found — PCA will use PLINK 1.9 fallback."
fi

# Python + sklearn for GMM clustering
python3 -c "import sklearn" 2>/dev/null \
    || crisp_die "Python sklearn not found. Required for GMM clustering. Install: pip install scikit-learn"

# 1KG reference check
if [[ "${RUN_1KG_PROJECTION^^}" == "YES" ]]; then
    [[ -n "${KG_PLINK_PREFIX}" ]] \
        || crisp_die "RUN_1KG_PROJECTION=YES but KG_PLINK_PREFIX is not set."
    [[ -f "${KG_PLINK_PREFIX}.bed" ]] \
        || crisp_die "1KG reference .bed not found: ${KG_PLINK_PREFIX}.bed"
fi

# Create output directories
mkdir -p "${STEP6_DIR}" "${PLOT_DIR}" "${LOG_DIR}"

# Initialise palette — exports CRISP_PAL_* for plot scripts
_crisp_palette "${PLOT_COLOUR_MODE}" "${PLOT_BACKGROUND}"

crisp_log "Validation complete."
crisp_log "Input   : ${STEP6_INPUT}"
crisp_log "Output  : ${STEP6_DIR}"
crisp_log "Assembly: ${REFERENCE_ASSEMBLY}"
crisp_log "Mode    : ASSIGN_POP=${ASSIGN_POP}  ASSIGN_SUB_POP=${ASSIGN_SUB_POP}  RUN_1KG_PROJECTION=${RUN_1KG_PROJECTION}"
crisp_log "Palette : ${PLOT_COLOUR_MODE} / ${PLOT_BACKGROUND}"

##########################################################################
#  _____ _                 ____
# /  ___| |               / ___|  __ _
# \ `--.| |_ ___ _ __    / /___  / _` |
#  `--. \ __/ _ \ '_ \   | ___ \| (_| |
# /\__/ / ||  __/ |_) |  | \_/ | \__,_|
# \____/ \__\___| .__/   \_____/
#               | |
#               |_|
##########################################################################
# 6a — LD PRUNING
# WHY: PCA on the full SNP set is dominated by regions of high LD (MHC,
# inversions, pericentromeric regions). Pruning to near-independent SNPs
# ensures PCs reflect population ancestry structure rather than local
# haplotype blocks. We use --indep-pairwise with the validated DIVERGE
# parameters: window=50 variants, step=5, r2=0.05. The --allow-extra-chr
# flag is non-negotiable: without it PLINK silently drops non-standard
# chromosomes which can include important variants.

crisp_log "6a — LD pruning (window=${PCA_LD_WINDOW}, step=${PCA_LD_STEP}, r2=${PCA_LD_R2})"

PRUNE_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_pruned"

${PLINK1} \
    --bfile "${STEP6_INPUT}" \
    --indep-pairwise "${PCA_LD_WINDOW}" "${PCA_LD_STEP}" "${PCA_LD_R2}" \
    --allow-extra-chr \
    --autosome \
    --out "${PRUNE_PREFIX}" \
    >> "${LOG_DIR}/step6a_prune.log" 2>&1

# WHY: --autosome restricts pruning to chromosomes 1-22. Sex chromosomes
# introduce sex-linked LD patterns that confound ancestry PCA; they are
# excluded from the pruned set used for ancestry inference. PC covariates
# for GWAS can include X-chr variants but ancestry assignment should not.

PRUNE_IN="${PRUNE_PREFIX}.prune.in"
N_PRUNED=$(wc -l < "${PRUNE_IN}" | tr -d ' ')
crisp_log "6a — Retained ${N_PRUNED} variants after LD pruning."

# Extract pruned variant set as a PLINK dataset (needed for MDS)
${PLINK1} \
    --bfile "${STEP6_INPUT}" \
    --extract "${PRUNE_IN}" \
    --allow-extra-chr \
    --make-bed \
    --out "${PRUNE_PREFIX}" \
    >> "${LOG_DIR}/step6a_extract.log" 2>&1

##########################################################################
#  _____ _                 ____ _
# /  ___| |               / ___| |__
# \ `--.| |_ ___ _ __    / /___| '_ \
#  `--. \ __/ _ \ '_ \   | ___ \ |_) |
# /\__/ / ||  __/ |_) |  | \_/ / |_) |
# \____/ \__\___| .__/   \_____/_.__/
#               | |
#               |_|
##########################################################################
# 6b — PCA COMPUTATION
# WHY: We prefer PLINK 2 for PCA — it uses a faster randomised SVD
# algorithm and produces cleaner output. PLINK 1.9 is the fallback.
# Both produce .eigenvec and .eigenval files in compatible formats.
# If 1KG projection is requested, we first merge the cohort with the
# 1KG reference panel on the overlapping pruned variant set, then run
# PCA on the combined dataset so cohort samples are projected into the
# 1KG PC space — giving named superpopulation coordinates (EUR, AFR, etc.)
# rather than anonymous group numbers.

PCA_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_pca"

if [[ "${RUN_1KG_PROJECTION^^}" == "YES" ]]; then
    crisp_log "6b — Merging with 1KG reference panel for named superpop projection..."

    # Extract overlapping variants between cohort pruned set and 1KG
    KG_MERGE_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_kg_merge"

    # Find SNPs present in both datasets via BIM intersection
    awk '{print $2}' "${PRUNE_PREFIX}.bim" | sort > "${STEP6_DIR}/cohort_snps.txt"
    awk '{print $2}' "${KG_PLINK_PREFIX}.bim" | sort > "${STEP6_DIR}/kg_snps.txt"
    comm -12 "${STEP6_DIR}/cohort_snps.txt" "${STEP6_DIR}/kg_snps.txt" \
        > "${STEP6_DIR}/shared_snps.txt"

    N_SHARED=$(wc -l < "${STEP6_DIR}/shared_snps.txt" | tr -d ' ')
    crisp_log "6b — ${N_SHARED} variants shared with 1KG reference."

    [[ "${N_SHARED}" -ge 1000 ]] \
        || crisp_warn "Fewer than 1000 shared variants with 1KG — projection may be unreliable."

    # Extract shared variants from both datasets then merge
    ${PLINK1} --bfile "${PRUNE_PREFIX}" \
        --extract "${STEP6_DIR}/shared_snps.txt" \
        --allow-extra-chr --make-bed \
        --out "${STEP6_DIR}/${PROJECT_NAME}_cohort_shared" \
        >> "${LOG_DIR}/step6b_shared.log" 2>&1

    ${PLINK1} --bfile "${KG_PLINK_PREFIX}" \
        --extract "${STEP6_DIR}/shared_snps.txt" \
        --allow-extra-chr --make-bed \
        --out "${STEP6_DIR}/${PROJECT_NAME}_kg_shared" \
        >> "${LOG_DIR}/step6b_kg_shared.log" 2>&1

    ${PLINK1} \
        --bfile "${STEP6_DIR}/${PROJECT_NAME}_cohort_shared" \
        --bmerge "${STEP6_DIR}/${PROJECT_NAME}_kg_shared" \
        --allow-extra-chr --make-bed \
        --out "${KG_MERGE_PREFIX}" \
        >> "${LOG_DIR}/step6b_merge.log" 2>&1

    PCA_INPUT="${KG_MERGE_PREFIX}"
    PCA_MODE="1KG_PROJECTION"
else
    PCA_INPUT="${PRUNE_PREFIX}"
    PCA_MODE="COHORT_ONLY"
fi

# Run PCA — PLINK 2 preferred, PLINK 1.9 fallback
if [[ -n "${PLINK2}" ]]; then
    crisp_log "6b — PCA via PLINK 2 (--pca ${PCA_N_PCS})"
    ${PLINK2} \
        --bfile "${PCA_INPUT}" \
        --pca "${PCA_N_PCS}" \
        --allow-extra-chr \
        --out "${PCA_PREFIX}" \
        >> "${LOG_DIR}/step6b_pca.log" 2>&1
    # WHY: PLINK 2 names output .eigenvec and .eigenval — same as PLINK 1.9.
    # No post-processing needed to unify formats.
else
    crisp_warn "6b — PCA via PLINK 1.9 fallback (--pca)"
    ${PLINK1} \
        --bfile "${PCA_INPUT}" \
        --pca "${PCA_N_PCS}" \
        --allow-extra-chr \
        --out "${PCA_PREFIX}" \
        >> "${LOG_DIR}/step6b_pca.log" 2>&1
fi

crisp_log "6b — PCA complete. Mode: ${PCA_MODE}"

##########################################################################
#  _____ _                 ____
# /  ___| |               / ___| ___
# \ `--.| |_ ___ _ __    / /___/ __|
#  `--. \ __/ _ \ '_ \   | ___ \ (__
# /\__/ / ||  __/ |_) |  | \_/ /\___|
# \____/ \__\___| .__/   \_____/
#               | |
#               |_|
##########################################################################
# 6c — MDS COMPUTATION
# WHY: MDS (Multidimensional Scaling) operates on a pairwise IBS distance
# matrix rather than the genotype covariance matrix used by PCA. For
# population structure the two methods produce very similar results, but
# running both and checking cluster assignment agreement is a valuable
# quality metric — high disagreement between PCA and MDS cluster labels
# signals that the structure may be complex or the clustering unstable.
# MDS is always run with PLINK 1.9 (--mds-plot); PLINK 2 dropped MDS.

MDS_PREFIX="${STEP6_DIR}/${PROJECT_NAME}_mds"

if [[ "${PCA_RUN_MDS_CROSSVAL^^}" == "YES" ]]; then
    crisp_log "6c — MDS computation (${MDS_N_DIMS} dimensions)"

    # Genome-wide IBS matrix then MDS
    ${PLINK1} \
        --bfile "${PRUNE_PREFIX}" \
        --read-genome "${STEP6_DIR}/${PROJECT_NAME}_ibs.genome" \
        --allow-extra-chr \
        --cluster --mds-plot "${MDS_N_DIMS}" \
        --out "${MDS_PREFIX}" \
        >> "${LOG_DIR}/step6c_mds.log" 2>&1 || {
        # WHY: --read-genome requires a pre-computed genome file. If absent,
        # compute IBS matrix first, then MDS in a second pass.
        crisp_warn "6c — No pre-computed IBS genome file. Computing IBS matrix first (slow for large N)."
        ${PLINK1} \
            --bfile "${PRUNE_PREFIX}" \
            --allow-extra-chr \
            --genome \
            --out "${STEP6_DIR}/${PROJECT_NAME}_ibs" \
            >> "${LOG_DIR}/step6c_ibs.log" 2>&1

        ${PLINK1} \
            --bfile "${PRUNE_PREFIX}" \
            --read-genome "${STEP6_DIR}/${PROJECT_NAME}_ibs.genome" \
            --allow-extra-chr \
            --cluster --mds-plot "${MDS_N_DIMS}" \
            --out "${MDS_PREFIX}" \
            >> "${LOG_DIR}/step6c_mds.log" 2>&1
    }

    crisp_log "6c — MDS complete."
    MDS_AVAILABLE="YES"
else
    crisp_log "6c — MDS skipped (PCA_RUN_MDS_CROSSVAL=NO)"
    MDS_AVAILABLE="NO"
fi

##########################################################################
#  _____ _                 ____     _
# /  ___| |               / ___|   | |
# \ `--.| |_ ___ _ __    / /___  __| |
#  `--. \ __/ _ \ '_ \   | ___ \/ _` |
# /\__/ / ||  __/ |_) |  | \_/ | (_| |
# \____/ \__\___| .__/   \_____/\__,_|
#               | |
#               |_|
##########################################################################
# 6d — GMM CLUSTERING + SD CHECK
# WHY: We use a two-method approach for cluster assignment:
#   1) GMM (Gaussian Mixture Model) assigns each sample to a cluster
#      and returns a posterior membership probability (MEMBERSHIP_P).
#      GMM is preferred over k-means because it handles elliptical
#      clusters and returns soft assignments — the probability value
#      is informative for admixed samples that sit between clusters.
#   2) SD check within each GMM cluster flags samples that are assigned
#      to a cluster but sit far from its centroid (> PCA_OUTLIER_SD SDs).
#      These are boundary cases and genuine outliers respectively.
# The combination answers two distinct questions:
#   GMM:  which group does this sample belong to?
#   SD:   is this sample a clean member of that group?

crisp_log "6d — GMM clustering + SD outlier check"

python3 "${SCRIPT_DIR}/scripts/crisp_gmm_cluster.py" \
    --eigenvec    "${PCA_PREFIX}.eigenvec" \
    --n-clusters  "${PCA_N_CLUSTERS}" \
    --outlier-sd  "${PCA_OUTLIER_SD}" \
    --pca-mode    "${PCA_MODE}" \
    --kg-prefix   "${KG_PLINK_PREFIX}" \
    --out-dir     "${STEP6_DIR}" \
    --project     "${PROJECT_NAME}" \
    2>&1 | tee "${LOG_DIR}/step6d_gmm.log"

# GMM script writes:
#   {PROJECT_NAME}_gmm_assignments.txt  (FID IID GROUP MEMBERSHIP_P OUTLIER)
#   {PROJECT_NAME}_cluster_{N}.txt      (FID IID — per-cluster ID lists)
#   {PROJECT_NAME}_gmm_summary.json     (cluster counts, K, BIC scores)

GMM_ASSIGN="${STEP6_DIR}/${PROJECT_NAME}_gmm_assignments.txt"
[[ -f "${GMM_ASSIGN}" ]] || crisp_die "GMM assignment file not produced: ${GMM_ASSIGN}"

N_CLUSTERS_DETECTED=$(python3 -c "
import json
with open('${STEP6_DIR}/${PROJECT_NAME}_gmm_summary.json') as f:
    d = json.load(f)
print(d['n_clusters'])
")
crisp_log "6d — Detected ${N_CLUSTERS_DETECTED} cluster(s)."

##########################################################################
#  _____ _                 ____
# /  ___| |               / ___| ___
# \ `--.| |_ ___ _ __    / /___/ _ \
#  `--. \ __/ _ \ '_ \   | ___ \  __/
# /\__/ / ||  __/ |_) |  | \_/ /\___|
# \____/ \__\___| .__/   \_____/
#               | |
#               |_|
##########################################################################
# 6e — PCA / MDS CROSS-VALIDATION
# WHY: If MDS was run, we compare cluster labels between PCA-GMM and
# MDS-GMM. High agreement (>90%) confirms the population structure is
# robust. Low agreement (<80%) is a warning that the structure is complex
# or unstable — the analyst should inspect both plots before proceeding.
# Agreement rate is written to the report and JSON, and recorded in the
# satellite file as MDS_AGREE (YES/NO per sample).

MDS_AGREE_FILE="${STEP6_DIR}/${PROJECT_NAME}_mds_agree.txt"

if [[ "${MDS_AVAILABLE}" == "YES" ]]; then
    crisp_log "6e — PCA/MDS cross-validation"

    python3 "${SCRIPT_DIR}/scripts/crisp_gmm_cluster.py" \
        --eigenvec    "${MDS_PREFIX}.mds" \
        --n-clusters  "${N_CLUSTERS_DETECTED}" \
        --outlier-sd  "${PCA_OUTLIER_SD}" \
        --pca-mode    "MDS" \
        --kg-prefix   "" \
        --out-dir     "${STEP6_DIR}" \
        --project     "${PROJECT_NAME}_mds" \
        2>&1 | tee "${LOG_DIR}/step6e_mds_gmm.log"

    # Compare PCA and MDS assignments per sample
    python3 "${SCRIPT_DIR}/scripts/crisp_crossval.py" \
        --pca-assign  "${GMM_ASSIGN}" \
        --mds-assign  "${STEP6_DIR}/${PROJECT_NAME}_mds_gmm_assignments.txt" \
        --out         "${MDS_AGREE_FILE}" \
        --project     "${PROJECT_NAME}" \
        2>&1 | tee "${LOG_DIR}/step6e_crossval.log"

    AGREE_RATE=$(python3 -c "
import json
with open('${STEP6_DIR}/${PROJECT_NAME}_crossval_summary.json') as f:
    d = json.load(f)
print(d['agreement_rate_pct'])
")
    crisp_log "6e — PCA/MDS cluster agreement: ${AGREE_RATE}%"

    if (( $(echo "${AGREE_RATE} < 80" | bc -l) )); then
        crisp_warn "PCA/MDS agreement below 80% (${AGREE_RATE}%). Population structure may be complex — inspect both plots."
    fi
else
    crisp_log "6e — Cross-validation skipped (MDS not run)."
    AGREE_RATE="NA"
fi

##########################################################################
#  _____ _                 ____  __
# /  ___| |               / ___|/ _|
# \ `--.| |_ ___ _ __    / /___| |_
#  `--. \ __/ _ \ '_ \   | ___ \  _|
# /\__/ / ||  __/ |_) |  | \_/ | |
# \____/ \__\___| .__/   \_____/_|
#               | |
#               |_|
##########################################################################
# 6f — ANCESTRY ASSIGNMENT & SATELLITE FILE
# WHY: The satellite file is the primary output of Step 6. It is a
# plain-text TSV with one row per sample, carrying all ancestry
# information produced by this step. Steps 7, 8, and 9 read from this
# file to apply population-stratified thresholds. It is designed to be
# updatable: Step 10 (PCA Pass 2) will overwrite CONFIDENCE from
# PROBABLE to CONFIDENT (or AMBIGUOUS) once the dataset is fully cleaned.
#
# ASSIGN_POP controls whether group labels are propagated to downstream
# steps — not whether the assignment computation runs. PCA always runs;
# the flag determines whether the satellite file is marked as active for
# downstream consumption.
#
# Confidence is always PROBABLE at this stage because the dataset still
# contains heterozygosity outliers, related samples, and HWE-failing
# variants — all of which can distort PC space slightly.

SATELLITE="${STEP6_DIR}/${PROJECT_NAME}_ancestry.txt"

crisp_log "6f — Writing satellite ancestry file"

python3 "${SCRIPT_DIR}/scripts/crisp_ancestry_satellite.py" \
    --gmm-assign      "${GMM_ASSIGN}" \
    --eigenvec        "${PCA_PREFIX}.eigenvec" \
    --mds-agree       "${MDS_AGREE_FILE}" \
    --maf-pop-tag     "${MONO_REF_POP}" \
    --assign-pop      "${ASSIGN_POP}" \
    --assign-sub-pop  "${ASSIGN_SUB_POP}" \
    --pca-mode        "${PCA_MODE}" \
    --n-pcs           "${PCA_N_PCS}" \
    --out             "${SATELLITE}" \
    --project         "${PROJECT_NAME}" \
    2>&1 | tee "${LOG_DIR}/step6f_satellite.log"

[[ -f "${SATELLITE}" ]] || crisp_die "Satellite file not produced: ${SATELLITE}"
N_SAMPLES=$(( $(wc -l < "${SATELLITE}") - 1 ))
N_OUTLIERS=$(awk -F'\t' 'NR>1 && $7=="YES" {n++} END {print n+0}' "${SATELLITE}")
crisp_log "6f — Satellite file written: ${N_SAMPLES} samples, ${N_OUTLIERS} outliers."

# Per-cluster PLINK datasets (one --keep call per detected group)
# WHY: Downstream steps (7, 8, 9) may need ancestry-stratified PLINK
# datasets for within-ancestry analyses. We produce them here so each
# downstream step can simply --bfile the relevant sub-cohort.
crisp_log "6f — Extracting per-cluster PLINK sub-datasets"

for i in $(seq 1 "${N_CLUSTERS_DETECTED}"); do
    CLUSTER_IDS="${STEP6_DIR}/${PROJECT_NAME}_cluster_${i}.txt"
    if [[ -f "${CLUSTER_IDS}" ]]; then
        ${PLINK1} \
            --bfile "${STEP6_INPUT}" \
            --keep  "${CLUSTER_IDS}" \
            --allow-extra-chr \
            --make-bed \
            --out   "${STEP6_DIR}/${PROJECT_NAME}_group${i}" \
            >> "${LOG_DIR}/step6f_cluster${i}.log" 2>&1
        crisp_log "6f — Group ${i} sub-dataset written."
    fi
done

##########################################################################
#  _____ _                 ____
# /  ___| |               / ___| __ _
# \ `--.| |_ ___ _ __    / /___  / _` |
#  `--. \ __/ _ \ '_ \   | ___ \| (_| |
# /\__/ / ||  __/ |_) |  | \_/ | \__, |
# \____/ \__\___| .__/   \_____/  __/ |
#               | |              |___/
#               |_|
##########################################################################
# 6g — PLOTTING
# WHY: Four plots are produced — monochrome and coloured versions for both
# PCA and MDS. The monochrome plot (Plot 1) shows cluster structure and
# ellipses without colour so the analyst can assess the geometry before
# any colour interpretation. The coloured plot (Plot 2) applies the CRISP
# palette with one colour per detected group. Outliers are always
# CRISP_PAL_FAIL regardless of mode.
#
# All six palette × background combinations (STANDARD/COLOURBLIND/NIGHT ×
# LIGHT/DARK) must render correctly. Colours are consumed exclusively from
# CRISP_PAL_* env vars — zero hardcoded hex values in either plot script.
#
# COHORT_ADMIXED=YES draws a confidence ellipse per cluster at
# PCA_ELLIPSE_SD standard deviations. The ellipse is SD-based for
# consistency with PCA_OUTLIER_SD. Mahalanobis ellipses deferred to v1.1.

_run_plots() {
    local ENGINE="$1"   # R or PYTHON
    local PLOT_TYPE="$2" # PCA or MDS
    local INPUT_FILE="$3"

    if [[ "${ENGINE}" == "R" ]]; then
        Rscript --vanilla \
            "${SCRIPT_DIR}/scripts/plot_pca_pass1.R" \
            --input       "${INPUT_FILE}" \
            --satellite   "${SATELLITE}" \
            --plot-type   "${PLOT_TYPE}" \
            --out-dir     "${PLOT_DIR}" \
            --project     "${PROJECT_NAME}" \
            --n-clusters  "${N_CLUSTERS_DETECTED}" \
            --ellipse-sd  "${PCA_ELLIPSE_SD}" \
            --admixed     "${COHORT_ADMIXED}" \
            >> "${LOG_DIR}/step6g_plot_${PLOT_TYPE,,}_r.log" 2>&1
    else
        python3 \
            "${SCRIPT_DIR}/scripts/plot_pca_pass1.py" \
            --input       "${INPUT_FILE}" \
            --satellite   "${SATELLITE}" \
            --plot-type   "${PLOT_TYPE}" \
            --out-dir     "${PLOT_DIR}" \
            --project     "${PROJECT_NAME}" \
            --n-clusters  "${N_CLUSTERS_DETECTED}" \
            --ellipse-sd  "${PCA_ELLIPSE_SD}" \
            --admixed     "${COHORT_ADMIXED}" \
            >> "${LOG_DIR}/step6g_plot_${PLOT_TYPE,,}_py.log" 2>&1
    fi
}

crisp_log "6g — Generating plots (engine: ${PLOT_ENGINE})"

# Export palette vars so plot scripts can read them
export CRISP_PAL_MODE="${PLOT_COLOUR_MODE}"
export CRISP_PAL_BACKGROUND="${PLOT_BACKGROUND}"

# PCA plots
if [[ "${PLOT_ENGINE}" == "R" || "${PLOT_ENGINE}" == "BOTH" ]]; then
    _run_plots "R" "PCA" "${PCA_PREFIX}.eigenvec"
fi
if [[ "${PLOT_ENGINE}" == "PYTHON" || "${PLOT_ENGINE}" == "BOTH" ]]; then
    _run_plots "PYTHON" "PCA" "${PCA_PREFIX}.eigenvec"
fi

# MDS plots
if [[ "${MDS_AVAILABLE}" == "YES" ]]; then
    if [[ "${PLOT_ENGINE}" == "R" || "${PLOT_ENGINE}" == "BOTH" ]]; then
        _run_plots "R" "MDS" "${MDS_PREFIX}.mds"
    fi
    if [[ "${PLOT_ENGINE}" == "PYTHON" || "${PLOT_ENGINE}" == "BOTH" ]]; then
        _run_plots "PYTHON" "MDS" "${MDS_PREFIX}.mds"
    fi
fi

crisp_log "6g — Plots written to ${PLOT_DIR}"

##########################################################################
#  _____ _                 ____ _
# /  ___| |               / ___| |__
# \ `--.| |_ ___ _ __    / /___| '_ \
#  `--. \ __/ _ \ '_ \   | ___ \ |_) |
# /\__/ / ||  __/ |_) |  | \_/ / |_) |
# \____/ \__\___| .__/   \_____/_.__/
#               | |
#               |_|
##########################################################################
# 6h — REPORT + JSON EXPORT
# WHY: Every CRISP step writes both a human-readable plain-text report
# and a machine-readable JSON summary. The JSON feeds COMPASS-AI and
# the Step 12 aggregated report. The plain-text report is the primary
# QC audit trail for the analyst.

crisp_log "6h — Writing report and JSON summary"

REPORT="${STEP6_DIR}/${PROJECT_NAME}_step6_pca_pass1_report.txt"
JSON="${STEP6_DIR}/${PROJECT_NAME}_step6_pca_pass1.json"

python3 - << PYEOF
import json, datetime, os

n_input   = int("$(wc -l < "${STEP6_INPUT}.fam" | tr -d ' ')")
n_pruned  = int("${N_PRUNED}")
n_clusters= int("${N_CLUSTERS_DETECTED}")
n_outliers= int("${N_OUTLIERS}")
n_samples = int("${N_SAMPLES}")
agree_rate= "${AGREE_RATE}"
pca_mode  = "${PCA_MODE}"
assign_pop= "${ASSIGN_POP}"
mds_run   = "${MDS_AVAILABLE}"
assembly  = "${REFERENCE_ASSEMBLY}"
project   = "${PROJECT_NAME}"
timestamp = datetime.datetime.utcnow().isoformat() + "Z"

# ── Plain-text report ──────────────────────────────────────────────────
lines = [
    "=" * 70,
    "  CRISP — Step 6 — PCA Pass 1 (Ancestry Assignment)",
    "=" * 70,
    f"  Project       : {project}",
    f"  Date          : {timestamp}",
    f"  Assembly      : {assembly}",
    f"  PCA mode      : {pca_mode}",
    f"  ASSIGN_POP    : {assign_pop}",
    "",
    "  Input",
    "-" * 70,
    f"  Samples (pre-Step 6) : {n_input}",
    f"  Variants after pruning: {n_pruned}",
    "",
    "  Clustering",
    "-" * 70,
    f"  Clusters detected    : {n_clusters}",
    f"  Outliers flagged     : {n_outliers}  ({n_outliers/n_samples*100:.1f}%)",
    f"  MDS cross-validation : {'YES — agreement ' + agree_rate + '%' if mds_run == 'YES' else 'NOT RUN'}",
    "",
    "  Satellite file",
    "-" * 70,
    f"  {project}_ancestry.txt",
    f"  Confidence tier      : PROBABLE (all assignments — PCA Pass 1)",
    f"  Columns              : FID IID SUPERPOP SUBPOP CONFIDENCE SOURCE",
    f"                         OUTLIER MEMBERSHIP_P MDS_AGREE MAF_POP_TAG",
    f"                         PC1..PC${os.environ.get('PCA_N_PCS','20')} ASSIGN_PASS",
    "",
    "  NOTE: Confidence will be upgraded to CONFIDENT by Step 10",
    "  (PCA Pass 2) after heterozygosity, relatedness, and HWE",
    "  filtering have been applied.",
    "=" * 70,
]

with open("${REPORT}", "w") as f:
    f.write("\n".join(lines) + "\n")

# ── JSON summary ───────────────────────────────────────────────────────
summary = {
    "schema_version"    : "1.0.0",
    "step"              : "6",
    "step_name"         : "pca_pass1",
    "project"           : project,
    "timestamp"         : timestamp,
    "reference_assembly": assembly,
    "pca_mode"          : pca_mode,
    "assign_pop"        : assign_pop,
    "n_samples_input"   : n_input,
    "n_variants_pruned" : n_pruned,
    "n_clusters"        : n_clusters,
    "n_outliers"        : n_outliers,
    "outlier_pct"       : round(n_outliers / n_samples * 100, 2) if n_samples > 0 else 0,
    "mds_crossval_run"  : mds_run == "YES",
    "mds_agreement_pct" : float(agree_rate) if agree_rate != "NA" else None,
    "confidence_tier"   : "PROBABLE",
    "satellite_file"    : "${SATELLITE}",
}

with open("${JSON}", "w") as f:
    json.dump(summary, f, indent=2)

print(f"Report : ${REPORT}")
print(f"JSON   : ${JSON}")
PYEOF

##########################################################################
# ______                      _
# | ___ \                    | |
# | |_/ /___ _ __   ___  _ __| |_
# |    // _ \ '_ \ / _ \| '__| __|
# | |\ \  __/ |_) | (_) | |  | |_
# \_| \_\___| .__/ \___/|_|   \__|
#           | |
#           |_|
##########################################################################

crisp_log ""
crisp_log "══════════════════════════════════════════════════════════════"
crisp_log "  Step 6 — PCA Pass 1 complete"
crisp_log "  Clusters   : ${N_CLUSTERS_DETECTED}"
crisp_log "  Outliers   : ${N_OUTLIERS}"
crisp_log "  MDS agree  : ${AGREE_RATE}%"
crisp_log "  Satellite  : ${SATELLITE}"
crisp_log "  → Steps 7, 8, 9 will read ancestry labels from satellite file"
crisp_log "══════════════════════════════════════════════════════════════"
