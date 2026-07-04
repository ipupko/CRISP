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
# Script : scripts/crisp_ancestry_satellite.py
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Assemble the Step 6 satellite ancestry file — the primary output
#   of PCA Pass 1 and the single source of truth for ancestry labels
#   consumed by Steps 7, 8, and 9.
#
#   The satellite file is designed as the ancestry equivalent of the
#   FAM file: one row per sample, plain-text TSV, versioned by
#   schema_version in the companion JSON, and updatable by downstream
#   steps (Step 10 PCA Pass 2 upgrades CONFIDENCE from PROBABLE to
#   CONFIDENT or AMBIGUOUS).
#
# Column schema (v1.0.0):
#   FID          — PLINK family ID
#   IID          — PLINK individual ID (stable per-sample identifier)
#   SUPERPOP     — Assigned population (e.g. EUR, GROUP_1)
#   SUBPOP       — Sub-population if ASSIGN_SUB_POP=YES, else —
#   CONFIDENCE   — Always PROBABLE at Step 6 (upgraded by Step 10)
#   SOURCE       — How the assignment was made (PCA_PASS1_COHORT or
#                  PCA_PASS1_1KG)
#   OUTLIER      — YES if sample > PCA_OUTLIER_SD from cluster centroid
#   MEMBERSHIP_P — GMM posterior probability of assigned cluster
#   MDS_AGREE    — YES/NO if PCA and MDS assignments agree (NA if MDS
#                  not run)
#   MAF_POP_TAG  — Population tag from Step 4 gnomAD comparison, or —
#   PC1..PCn     — Raw PC scores from .eigenvec
#   ASSIGN_PASS  — 1 = assignment active and propagated downstream
#                  0 = PCA ran but ASSIGN_POP=NO; labels not propagated
#
# WHY ASSIGN_PASS:
#   ASSIGN_POP controls whether population labels are used downstream,
#   not whether PCA runs. ASSIGN_PASS=0 means the satellite file is
#   written (for reference and auditability) but Steps 7/8/9 treat all
#   samples as a single unstratified group. ASSIGN_PASS=1 means Steps
#   7/8/9 read SUPERPOP and apply stratified thresholds.
#
# Usage:
#   python3 crisp_ancestry_satellite.py \
#       --gmm-assign    sim_project_gmm_assignments.txt \
#       --eigenvec      sim_project_pca.eigenvec \
#       --mds-agree     sim_project_mds_agree.txt \
#       --maf-pop-tag   EUR \
#       --assign-pop    YES \
#       --assign-sub-pop NO \
#       --pca-mode      COHORT_ONLY \
#       --n-pcs         20 \
#       --out           sim_project_ancestry.txt \
#       --project       sim_project
#
##########################################################################

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd


##########################################################################
# ARGUMENT PARSING
##########################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="CRISP Step 6 — assemble satellite ancestry file"
    )
    p.add_argument("--gmm-assign",     required=True,
                   help="GMM assignment file from crisp_gmm_cluster.py")
    p.add_argument("--eigenvec",       required=True,
                   help="PLINK .eigenvec file with raw PC scores")
    p.add_argument("--mds-agree",      default="",
                   help="MDS agreement file from crisp_crossval.py (optional)")
    p.add_argument("--maf-pop-tag",    default="",
                   help="Population tag from Step 4 MAF comparison (or empty)")
    p.add_argument("--assign-pop",     default="NO",
                   help="ASSIGN_POP flag from instruction file (YES/NO)")
    p.add_argument("--assign-sub-pop", default="NO",
                   help="ASSIGN_SUB_POP flag from instruction file (YES/NO)")
    p.add_argument("--pca-mode",       default="COHORT_ONLY",
                   help="PCA computation mode (COHORT_ONLY or 1KG_PROJECTION)")
    p.add_argument("--n-pcs",          type=int, default=20,
                   help="Number of PCs to include in satellite file (default: 20)")
    p.add_argument("--out",            required=True,
                   help="Output satellite file path")
    p.add_argument("--project",        required=True,
                   help="Project name for JSON summary")
    return p.parse_args()


##########################################################################
# FILE READERS
##########################################################################

def read_eigenvec(path: str, n_pcs: int) -> pd.DataFrame:
    """
    Read PLINK .eigenvec and return FID, IID, PC1..PCn.
    Handles both PLINK 1.9 (no header) and PLINK 2 (#FID header).
    WHY: PC scores are included in the satellite file so downstream
    steps and analysts can read all ancestry information from a single
    file without needing to re-read the .eigenvec separately.
    """
    with open(path) as f:
        first = f.readline().strip()

    if first.startswith("#"):
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = [c.lstrip("#") for c in df.columns]
    else:
        df = pd.read_csv(path, sep=r"\s+", header=None)
        n_cols = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i+1}" for i in range(n_cols)]

    # Keep only requested PCs
    pc_cols = [f"PC{i+1}" for i in range(n_pcs) if f"PC{i+1}" in df.columns]
    return df[["FID", "IID"] + pc_cols]


##########################################################################
# MAIN
##########################################################################

def main():
    args = parse_args()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[CRISP] crisp_ancestry_satellite.py — {args.project}")
    print(f"[CRISP] GMM assignments : {args.gmm_assign}")
    print(f"[CRISP] Eigenvec        : {args.eigenvec}")
    print(f"[CRISP] ASSIGN_POP      : {args.assign_pop}")
    print(f"[CRISP] ASSIGN_SUB_POP  : {args.assign_sub_pop}")
    print(f"[CRISP] PCA mode        : {args.pca_mode}")

    # ASSIGN_PASS — 1 if ASSIGN_POP=YES, 0 otherwise
    assign_pass = 1 if args.assign_pop.upper() == "YES" else 0

    # Source label for the satellite file
    source = (
        "PCA_PASS1_1KG"    if args.pca_mode == "1KG_PROJECTION"
        else "PCA_PASS1_COHORT"
    )

    # ── Read GMM assignments ──────────────────────────────────────────
    gmm = pd.read_csv(args.gmm_assign, sep="\t")
    required = ["FID", "IID", "GROUP", "MEMBERSHIP_P", "OUTLIER"]
    missing = [c for c in required if c not in gmm.columns]
    if missing:
        print(f"[CRISP] ERROR: GMM assignment file missing columns: {missing}",
              file=sys.stderr)
        sys.exit(1)

    print(f"[CRISP] Samples from GMM: {len(gmm)}")

    # ── Read PC scores ────────────────────────────────────────────────
    eigenvec = read_eigenvec(args.eigenvec, args.n_pcs)
    pc_cols  = [c for c in eigenvec.columns if c.startswith("PC")]

    print(f"[CRISP] PC columns      : {len(pc_cols)} (PC1..PC{len(pc_cols)})")

    # ── Read MDS agreement ────────────────────────────────────────────
    # WHY: MDS agreement is optional — only present when PCA_RUN_MDS_CROSSVAL=YES.
    # If absent, MDS_AGREE column is filled with NA for all samples.
    if args.mds_agree and os.path.exists(args.mds_agree):
        mds_agree = pd.read_csv(args.mds_agree, sep="\t")
        print(f"[CRISP] MDS agreement   : {len(mds_agree)} samples")
    else:
        mds_agree = None
        print(f"[CRISP] MDS agreement   : not available (MDS not run)")

    # ── Assemble satellite DataFrame ──────────────────────────────────
    # Start from GMM assignments — this is the canonical sample list
    sat = gmm[["FID", "IID", "GROUP", "MEMBERSHIP_P", "OUTLIER"]].copy()

    # SUPERPOP — the GMM group label (named superpop or GROUP_N)
    sat = sat.rename(columns={"GROUP": "SUPERPOP"})

    # SUBPOP — populated only when ASSIGN_SUB_POP=YES
    # WHY: Sub-population detection (e.g. KPP vs Punjab within SAS)
    # requires a second GMM pass within each cluster. At v1.0 with
    # ASSIGN_SUB_POP=NO, we write — as a placeholder. The column is
    # present so downstream tools can rely on the schema without
    # checking whether the flag was set.
    sat["SUBPOP"] = "—"

    # CONFIDENCE — always PROBABLE at Step 6 (upgraded by Step 10)
    sat["CONFIDENCE"] = "PROBABLE"

    # SOURCE — how the assignment was made
    sat["SOURCE"] = source

    # Reorder to match schema: FID IID SUPERPOP SUBPOP CONFIDENCE SOURCE OUTLIER
    sat = sat[["FID", "IID", "SUPERPOP", "SUBPOP",
               "CONFIDENCE", "SOURCE", "OUTLIER", "MEMBERSHIP_P"]]

    # MDS_AGREE — per-sample cross-validation result
    if mds_agree is not None:
        sat = sat.merge(
            mds_agree[["IID", "MDS_AGREE"]],
            on="IID", how="left"
        )
        sat["MDS_AGREE"] = sat["MDS_AGREE"].fillna("NA")
    else:
        sat["MDS_AGREE"] = "NA"

    # MAF_POP_TAG — Step 4 population tag for concordance reference
    # WHY: Preserving the Step 4 DB-derived population tag alongside the
    # Step 6 PCA-derived assignment creates a two-point evidence chain.
    # Analysts can check whether the reference DB and PCA agree on
    # population identity, which adds confidence when they match and
    # flags discordance when they don't.
    maf_tag = args.maf_pop_tag.strip() if args.maf_pop_tag else "—"
    sat["MAF_POP_TAG"] = maf_tag

    # PC scores — merge from eigenvec on IID
    sat = sat.merge(eigenvec[["IID"] + pc_cols], on="IID", how="left")

    # ASSIGN_PASS — flag controlling downstream propagation
    sat["ASSIGN_PASS"] = assign_pass

    # ── Round PC columns to 6 decimal places ─────────────────────────
    # WHY: Raw PC scores have 15+ decimal places of floating-point noise.
    # Rounding to 6 DP retains all meaningful precision while keeping
    # the file human-readable and reducing file size.
    for pc in pc_cols:
        sat[pc] = sat[pc].round(6)

    # ── Round MEMBERSHIP_P ────────────────────────────────────────────
    sat["MEMBERSHIP_P"] = sat["MEMBERSHIP_P"].round(4)

    # ── Write satellite file ──────────────────────────────────────────
    sat.to_csv(out_path, sep="\t", index=False)
    print(f"[CRISP] Satellite file written: {out_path}")
    print(f"[CRISP] Rows    : {len(sat)}")
    print(f"[CRISP] Columns : {list(sat.columns)}")

    # ── Summary statistics ────────────────────────────────────────────
    n_total    = len(sat)
    n_outliers = (sat["OUTLIER"] == "YES").sum()
    n_assigned = n_total  # all samples receive an assignment
    superpop_counts = sat["SUPERPOP"].value_counts().to_dict()

    if mds_agree is not None:
        n_mds_agree    = (sat["MDS_AGREE"] == "YES").sum()
        n_mds_disagree = (sat["MDS_AGREE"] == "NO").sum()
    else:
        n_mds_agree = n_mds_disagree = 0

    print(f"[CRISP] Population breakdown:")
    for pop, n in sorted(superpop_counts.items()):
        print(f"[CRISP]   {pop}: {n} samples")

    # ── Write companion JSON ──────────────────────────────────────────
    # WHY: The companion JSON carries the schema_version so that future
    # pipeline versions can detect format changes and handle them
    # gracefully, even if the satellite file format has evolved.
    summary = {
        "schema_version"    : "1.0.0",
        "project"           : args.project,
        "satellite_file"    : str(out_path),
        "pca_mode"          : args.pca_mode,
        "assign_pop"        : args.assign_pop.upper(),
        "assign_sub_pop"    : args.assign_sub_pop.upper(),
        "assign_pass"       : assign_pass,
        "confidence_tier"   : "PROBABLE",
        "source"            : source,
        "n_samples"         : n_total,
        "n_outliers"        : int(n_outliers),
        "n_pc_columns"      : len(pc_cols),
        "maf_pop_tag"       : maf_tag,
        "mds_crossval_run"  : mds_agree is not None,
        "n_mds_agree"       : int(n_mds_agree),
        "n_mds_disagree"    : int(n_mds_disagree),
        "superpop_counts"   : {k: int(v) for k, v in superpop_counts.items()},
        "columns"           : list(sat.columns),
    }

    json_path = out_path.parent / f"{args.project}_ancestry_satellite.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[CRISP] Companion JSON   : {json_path}")
    print(f"[CRISP] crisp_ancestry_satellite.py complete.")


if __name__ == "__main__":
    main()
