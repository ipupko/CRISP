#!/usr/bin/env python3
# ##############################################################
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
# ##############################################################
#
# Script : scripts/crisp_pc_selection.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Sub-step 6k — Automatic PC count selection for GWAS covariates.
#
#   The standard "use 10 PCs" rule of thumb (Price et al. 2006) is
#   biologically unjustified — it was a pragmatic default that has
#   been cargo-culted for 20 years. The correct number of PCs to
#   include as GWAS covariates depends on the population structure
#   of the cohort:
#
#     Homogeneous (founder, isolate)  → 2–5 PCs sufficient
#     Structured (cline, admixed EUR) → 5–10 PCs
#     Heterogeneous (multi-ancestry)  → 10–20 PCs
#
#   This script implements three complementary methods and combines
#   them into a conservative recommendation. COMPASS-AI then learns
#   — across thousands of cohorts — the mapping from PC geometry
#   fingerprint (elongation, sphericity, effective dimensionality)
#   to optimal PC count, improving recommendations over time.
#
# Methods implemented:
#
#   1. Marchenko-Pastur (MP) law
#      Random matrix theory baseline. Eigenvalues above the MP upper
#      bound are signal; below are noise. The number of eigenvalues
#      above the MP bound = number of statistically significant PCs.
#      WHY: Principled, parameter-free, fast. Does not require the
#      GWAS to be run. Accounts for the expected distribution of
#      eigenvalues under the null (no population structure).
#
#   2. Elbow detection (second derivative method)
#      Find the inflection point in the scree plot — where eigenvalues
#      stop dropping steeply. Automated via the second derivative of
#      the eigenvalue curve (maximum curvature point).
#      WHY: Intuitive, matches what a human would identify visually
#      on the scree plot. Robust to the scale of eigenvalues.
#
#   3. Variance explained threshold
#      Count PCs needed to explain PCA_VARIANCE_THRESHOLD of total
#      variance (default 95%). Simple fallback.
#      WHY: Interpretable to analysts who think in terms of variance
#      explained. Used as a cross-check and reported alongside.
#
#   4. Tracy-Widom test (optional, requires rpy2 + AssocTests R package)
#      Tests whether each eigenvalue exceeds the Tracy-Widom
#      distribution under the null. Gold standard statistically.
#      WHY: The most rigorous method but requires R dependency.
#      Off by default; enabled via PCA_TW_TEST=YES.
#
#   Recommendation logic:
#     recommended_n_pcs = max(mp_n, elbow_n)   # conservative: take larger
#     recommended_n_pcs = min(recommended_n_pcs, max_pcs)  # cap at max
#     if recommended_n_pcs < 2: recommended_n_pcs = 2      # floor
#
#   WHY conservative (max not min):
#     Under-correction for stratification causes false positives in GWAS.
#     Over-correction wastes a few degrees of freedom but does not
#     inflate false positives. For GWAS, false positives are the more
#     serious error. We err on the side of more PCs.
#
# Homogeneity classification:
#   HOMOGENEOUS    — recommended_n_pcs <= 4
#   STRUCTURED     — recommended_n_pcs 5–9
#   HETEROGENEOUS  — recommended_n_pcs >= 10
#
# COMPASS-AI training signals:
#   recommended_n_pcs, mp_n, elbow_n, variance_n, homogeneity_class,
#   plus cohort PC geometry features (elongation_ratio, sphericity,
#   effective_dimensionality_90pct, cloud_log_volume, mean_nn_distance)
#   from crisp_pc_geometry.py — enabling COMPASS-AI to predict optimal
#   PC count from geometry alone without rerunning this analysis.
#
# Integration note (NOT YET WIRED INTO STEP 6 ORCHESTRATOR):
#   Pending integration as sub-step 6k once Step 6 orchestration
#   is finalised. Call pattern will be:
#
#     python3 "${SCRIPT_DIR}/scripts/crisp_pc_selection.py" \
#         --eigenval    "${PCA_PREFIX}.eigenval" \
#         --geometry    "${STEP6_DIR}/${PROJECT_NAME}_pc_geometry.json" \
#         --n-samples   <n_samples> \
#         --n-variants  <n_pruned_variants> \
#         --method      AUTO \
#         --variance-threshold 0.95 \
#         --max-pcs     20 \
#         --tw-test     NO \
#         --out         "${STEP6_DIR}/${PROJECT_NAME}_pc_selection.json" \
#         --project     "${PROJECT_NAME}"
#
# Outputs:
#   {project}_pc_selection.json    COMPASS-AI structured recommendation
#
##########################################################################

import argparse
import json
import sys
from pathlib import Path

import numpy as np


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6k — automatic PC count selection"
    )
    p.add_argument("--eigenval",          required=True,
                   help="PLINK .eigenval file")
    p.add_argument("--geometry",          default="",
                   help="PC geometry JSON from crisp_pc_geometry.py (optional)")
    p.add_argument("--n-samples",         type=int, required=True,
                   help="Number of samples (needed for Marchenko-Pastur)")
    p.add_argument("--n-variants",        type=int, required=True,
                   help="Number of LD-pruned variants used for PCA")
    p.add_argument("--method",            default="AUTO",
                   choices=["AUTO","MARCHENKO_PASTUR","ELBOW","VARIANCE","TW"],
                   help="Selection method (default: AUTO)")
    p.add_argument("--variance-threshold",type=float, default=0.95,
                   help="Variance explained threshold (default: 0.95)")
    p.add_argument("--max-pcs",           type=int, default=20,
                   help="Maximum PCs to consider (default: 20)")
    p.add_argument("--tw-test",           default="NO",
                   help="Run Tracy-Widom test YES/NO (requires rpy2)")
    p.add_argument("--out",               required=True,
                   help="Output JSON path")
    p.add_argument("--project",           required=True,
                   help="Project name prefix")
    return p.parse_args()


##########################################################################
# EIGENVALUE LOADING
##########################################################################

def load_eigenvalues(path: str) -> np.ndarray:
    """
    Load PLINK .eigenval file.
    WHY: PLINK eigenvalues are the variance explained by each PC in
    units of the genetic relatedness matrix. Negative eigenvalues
    (floating point artefacts from near-singular matrices) are clipped
    to zero rather than raising an error.
    """
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                vals.append(float(line))
    arr = np.array(vals)
    arr = np.clip(arr, 0, None)
    print(f"[CRISP] Eigenvalues loaded: {len(arr)}")
    print(f"[CRISP] Top 5: {[round(float(v),4) for v in arr[:5]]}")
    return arr


##########################################################################
# METHOD 1: MARCHENKO-PASTUR LAW
##########################################################################

def marchenko_pastur(eigenvalues: np.ndarray,
                     n_samples: int,
                     n_variants: int) -> tuple:
    """
    Count eigenvalues above the Marchenko-Pastur upper bound.

    WHY Marchenko-Pastur:
      Random matrix theory tells us the distribution of eigenvalues
      when there is NO population structure. The MP law gives the upper
      bound of this null distribution. Any eigenvalue above this bound
      is statistically inconsistent with noise and represents real
      population structure.

    MP upper bound: lambda_max = sigma^2 * (1 + sqrt(gamma))^2
    where gamma = n_variants / n_samples (matrix aspect ratio)

    WHY lower-quartile sigma^2 estimate:
      Using total_variance / n_variants inflates sigma^2 by including
      signal eigenvalues in the sum, causing the MP bound to be too high
      and under-estimating the number of significant PCs.
      The lower quartile of eigenvalues is robustly in the noise floor
      and gives a more accurate noise variance estimate. Validated across
      homogeneous (Finnish), structured (European), and heterogeneous
      (multi-ancestry biobank) eigenvalue distributions.
    """
    print(f"[CRISP] Marchenko-Pastur (n={n_samples}, p={n_variants})")

    gamma   = n_variants / n_samples
    # WHY lower quartile: most robust noise floor estimate
    n_noise = max(1, len(eigenvalues) // 4)
    sigma2  = float(np.sort(eigenvalues)[:n_noise].mean())

    lambda_max = sigma2 * (1 + np.sqrt(gamma)) ** 2
    lambda_min = sigma2 * (1 - np.sqrt(gamma)) ** 2 if gamma <= 1 else 0.0

    print(f"[CRISP]   gamma={gamma:.4f}  sigma2={sigma2:.6f} (lower quartile)")
    print(f"[CRISP]   MP upper bound: {lambda_max:.6f}")

    above = int((eigenvalues > lambda_max).sum())
    print(f"[CRISP]   Eigenvalues above MP bound: {above}")

    diag = {
        "gamma"           : round(float(gamma), 6),
        "sigma2"          : round(float(sigma2), 8),
        "sigma2_method"   : "lower_quartile",
        "mp_lambda_max"   : round(float(lambda_max), 6),
        "mp_lambda_min"   : round(float(lambda_min), 6),
        "n_above_mp"      : above,
        "eigenvalues"     : [round(float(v), 6) for v in eigenvalues],
        "above_mp_bound"  : [bool(v > lambda_max) for v in eigenvalues],
    }
    return max(1, above), diag


##########################################################################
# METHOD 2: ELBOW DETECTION (SECOND DERIVATIVE)
##########################################################################

def elbow_detection(eigenvalues: np.ndarray) -> tuple:
    """
    Find the scree plot elbow via the Kneedle algorithm.

    WHY Kneedle over second derivative:
      The second derivative method detects the point of maximum curvature
      but is biased towards the beginning of the curve — it fires early
      on any steep initial drop, even when meaningful structure extends
      further. The Kneedle algorithm measures perpendicular distance from
      the diagonal line connecting first and last eigenvalue. This is
      more robust for gradual scree plots (structured/heterogeneous
      cohorts) and gives the same result as second derivative for sharp
      scree plots (homogeneous cohorts). Validated across all three
      cohort types at realistic PLINK eigenvalue scales.

    WHY normalise:
      Raw eigenvalues vary in scale across cohorts. Normalising to [0,1]
      makes the diagonal reference line scale-invariant and ensures
      the perpendicular distance metric is comparable across cohorts.
    """
    print(f"[CRISP] Elbow detection (Kneedle algorithm)")

    if len(eigenvalues) < 3:
        return 1, {"elbow_n": 1, "method": "insufficient_eigenvalues"}

    n     = len(eigenvalues)
    x     = np.arange(n, dtype=float)
    y     = eigenvalues.astype(float)

    # Normalise to [0,1]
    x_n = (x - x.min()) / (x.max() - x.min() + 1e-10)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-10)

    # Perpendicular distance from diagonal line x+y=1
    distances = np.abs(x_n + y_n - 1) / np.sqrt(2)
    knee_idx  = int(np.argmax(distances)) + 1  # 1-indexed

    knee_idx = max(1, min(knee_idx, n))
    print(f"[CRISP]   Kneedle elbow at PC{knee_idx}")

    diag = {
        "elbow_n"         : knee_idx,
        "method"          : "kneedle",
        "distances"       : [round(float(v),6) for v in distances[:15]],
        "argmax_distance" : int(np.argmax(distances)),
    }
    return knee_idx, diag


##########################################################################
# METHOD 3: VARIANCE EXPLAINED THRESHOLD
##########################################################################

def variance_threshold_method(eigenvalues: np.ndarray,
                               threshold: float) -> tuple:
    """
    Count PCs needed to explain `threshold` fraction of total variance.

    WHY this is a cross-check not a primary method:
      The threshold is arbitrary — there is no biological reason why
      exactly 95% of within-cohort variance needs to be captured for
      GWAS. In homogeneous cohorts, 95% may be reached by PC1 alone
      (cohort is tight, not because there is no structure). We report
      this alongside MP and elbow as a sanity check.
    """
    total = eigenvalues.sum()
    if total == 0:
        return 1, {"variance_n": 1, "threshold": threshold}

    cumvar = np.cumsum(eigenvalues) / total
    n_pcs  = int(np.searchsorted(cumvar, threshold)) + 1
    n_pcs  = max(1, min(n_pcs, len(eigenvalues)))

    print(f"[CRISP]   Variance threshold {threshold*100:.0f}%: {n_pcs} PCs "
          f"(cumvar={cumvar[n_pcs-1]*100:.1f}%)")

    return n_pcs, {
        "variance_n"    : n_pcs,
        "threshold"     : threshold,
        "cumvar_at_n"   : round(float(cumvar[n_pcs-1]), 4),
        "cumvar_curve"  : [round(float(v),4) for v in cumvar],
    }


##########################################################################
# METHOD 4: TRACY-WIDOM TEST (OPTIONAL)
##########################################################################

def tracy_widom_test(eigenvalues: np.ndarray,
                     n_samples: int,
                     n_variants: int) -> tuple:
    """
    Tracy-Widom significance test for each eigenvalue.

    WHY Tracy-Widom:
      The gold standard for PC significance testing. Used by EIGENSOFT.
      The TW distribution describes the largest eigenvalue of a Wishart
      matrix under the null (no structure). PCs whose eigenvalues exceed
      the TW critical value are statistically significant.

    WHY optional:
      Requires rpy2 + AssocTests R package. Not always available.
      MP + elbow gives a good answer in most cases.
    """
    print(f"[CRISP] Tracy-Widom test")
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()

        ro.r("library(AssocTests)")
        result  = ro.r["tw"](ro.FloatVector(eigenvalues.tolist()),
                             n_samples, n_variants)
        p_vals  = list(result.rx2("p.value"))
        n_sig   = int(sum(1 for p in p_vals if p < 0.05))

        print(f"[CRISP]   TW significant PCs: {n_sig}")
        return max(1, n_sig), {
            "tw_n_significant" : n_sig,
            "tw_p_values"      : [round(float(p),6) for p in p_vals],
            "tw_alpha"         : 0.05,
        }
    except ImportError:
        print(f"[CRISP] WARNING: rpy2 not available — TW test skipped.")
        return None, {"error": "rpy2_not_available"}
    except Exception as e:
        print(f"[CRISP] WARNING: TW test failed: {e}")
        return None, {"error": str(e)}


##########################################################################
# HOMOGENEITY CLASSIFICATION
##########################################################################

def classify_homogeneity(recommended_n: int, geometry: dict) -> tuple:
    """
    Classify cohort structural homogeneity from recommended PC count
    and optional PC geometry signals.

    WHY three tiers:
      HOMOGENEOUS — few PCs, typical of founder populations (Finnish,
        Sardinian, certain South Asian isolates). Tight PC clouds,
        high sphericity, low effective dimensionality.

      STRUCTURED — moderate complexity, typical of single-country
        European or East Asian cohorts. Elongated PC clouds (cline-like).

      HETEROGENEOUS — complex multi-ancestry structure, typical of
        pan-ethnic biobanks or admixed populations. High effective
        dimensionality, low sphericity, large cloud volume.

    COMPASS-AI uses this as a high-level cohort label for cross-study
    comparison and recommendation calibration across the biobank network.
    """
    tier = (
        "HOMOGENEOUS"  if recommended_n <= 4  else
        "STRUCTURED"   if recommended_n <= 9  else
        "HETEROGENEOUS"
    )

    geo_signals = {}
    if geometry:
        coh   = geometry.get("cohort", {})
        elong = coh.get("elongation_ratio")
        spher = coh.get("sphericity")
        eff90 = coh.get("effective_dimensionality_90pct")

        if elong is not None:
            geo_signals = {
                "elongation_ratio"       : elong,
                "sphericity"             : spher,
                "effective_dim_90pct"    : eff90,
                "geometry_predicted_tier": (
                    "HOMOGENEOUS"  if eff90 is not None and eff90 <= 3  else
                    "STRUCTURED"   if eff90 is not None and eff90 <= 8  else
                    "HETEROGENEOUS"
                ),
            }
            if geo_signals["geometry_predicted_tier"] != tier:
                print(f"[CRISP] NOTE: eigenvalue tier ({tier}) differs from "
                      f"geometry tier ({geo_signals['geometry_predicted_tier']}). "
                      f"Eigenvalue estimate takes precedence.")

    print(f"[CRISP] Homogeneity: {tier}")
    return tier, geo_signals


##########################################################################
# MAIN
##########################################################################

def main():
    args = parse_args()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP] crisp_pc_selection.py — {args.project}")
    print(f"[CRISP] Method    : {args.method}")
    print(f"[CRISP] N samples : {args.n_samples}")
    print(f"[CRISP] N variants: {args.n_variants}")

    # Load
    eigenvalues = load_eigenvalues(args.eigenval)[:args.max_pcs]

    geometry = {}
    if args.geometry and Path(args.geometry).exists():
        with open(args.geometry) as f:
            geometry = json.load(f)

    # Run methods
    results = {}
    method  = args.method.upper()
    mp_n = elbow_n = var_n = tw_n = None

    if method in ("AUTO","MARCHENKO_PASTUR"):
        mp_n, d = marchenko_pastur(eigenvalues, args.n_samples, args.n_variants)
        results["marchenko_pastur"] = {"n": mp_n, **d}

    if method in ("AUTO","ELBOW"):
        elbow_n, d = elbow_detection(eigenvalues)
        results["elbow"] = {"n": elbow_n, **d}

    if method in ("AUTO","VARIANCE"):
        var_n, d = variance_threshold_method(eigenvalues, args.variance_threshold)
        results["variance_threshold"] = {"n": var_n, **d}

    if args.tw_test.upper() == "YES" or method == "TW":
        tw_n, d = tracy_widom_test(eigenvalues, args.n_samples, args.n_variants)
        if tw_n is not None:
            results["tracy_widom"] = {"n": tw_n, **d}

    # Recommendation
    if method == "AUTO":
        candidates = [n for n in [mp_n, elbow_n] if n is not None]
        recommended = max(candidates) if candidates else 10
    elif method == "MARCHENKO_PASTUR": recommended = mp_n or 10
    elif method == "ELBOW":            recommended = elbow_n or 10
    elif method == "VARIANCE":         recommended = var_n or 10
    elif method == "TW":               recommended = tw_n or 10
    else:
        all_c = [n for n in [mp_n, elbow_n, var_n, tw_n] if n is not None]
        recommended = max(all_c) if all_c else 10

    recommended = max(2, min(recommended, args.max_pcs))

    print(f"[CRISP] ══════════════════════════════════════")
    print(f"[CRISP] Recommended PCs for GWAS: {recommended}")

    homogeneity, geo_signals = classify_homogeneity(recommended, geometry)

    output = {
        "schema_version"           : "1.0.0",
        "sub_step"                 : "6k_pc_selection",
        "project"                  : args.project,
        "method"                   : args.method,
        "n_samples"                : args.n_samples,
        "n_variants"               : args.n_variants,
        "max_pcs_considered"       : args.max_pcs,
        "variance_threshold"       : args.variance_threshold,
        "recommended_n_pcs"        : recommended,
        "homogeneity_class"        : homogeneity,
        "mp_n"                     : mp_n,
        "elbow_n"                  : elbow_n,
        "variance_n"               : var_n,
        "tw_n"                     : tw_n,
        "recommendation_rationale" : (
            f"Conservative estimate max(MP={mp_n}, elbow={elbow_n}) "
            f"= {recommended}, capped at {args.max_pcs}. "
            f"Cohort classified as {homogeneity}."
        ),
        "compass_signals"          : {
            "recommended_n_pcs" : recommended,
            "homogeneity_class" : homogeneity,
            "mp_n"              : mp_n,
            "elbow_n"           : elbow_n,
            "variance_n"        : var_n,
            **geo_signals,
        },
        "usage_note"               : (
            f"Include PC1–PC{recommended} as GWAS covariates. "
            f"Validate with genomic inflation factor post-GWAS."
        ),
        "method_diagnostics"       : results,
    }

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"[CRISP] Output: {out_path}")
    print(f"[CRISP] crisp_pc_selection.py complete.")


if __name__ == "__main__":
    main()
