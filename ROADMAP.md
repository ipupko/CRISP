# CRISP Roadmap
**Comprehensive Robust Integrated SNP Processing**

This document outlines the planned development trajectory for CRISP.
It is a living document and will be updated as the project evolves.
Community feedback and contributions are welcome at every stage.

---

## Current Status: Active Development (pre v1.0)

The following chunks are built and tested:

| Chunk | Step | Status |
|-------|------|--------|
| `crisp_upload_chunk.sh` | Step 1: File validation and MD5 | Complete |
| `crisp_convert.sh` | Step 2: Format conversion | Complete |
| `crisp_callrate.sh` | Step 3: Sample call rate | Complete |
| `crisp_snprate.sh` | Step 4: Variant call rate, monomorphic, MAF | Complete |
| `crisp_sexcheck.sh` | Step 5: Sex check and aneuploidy detection | Complete |
| `crisp_homozygosity.sh` | Step 6: Heterozygosity and homozygosity | Complete |
| `crisp_relatedness.sh` | Step 7: Relatedness and duplicates | In progress |
| `crisp_hwe.sh` | Step 8: HWE filtering | Pending |
| `crisp_pca.sh` | Step 9: PCA and ancestry | Pending |
| `crisp_amend.sh` | Step 10: Amendments and clean dataset | Pending |
| `crisp_report.sh` | Step 11: End-to-end summary report | Pending |

**Plotting engine parity (R and Python) achieved for Steps 1 through 6.**

---

## Version 1.0: Full Pipeline

**Target:** First fully working end-to-end pipeline release.

All steps complete and tested. Every step produces:
- A structured plain-text report
- PDF plots via R (default) or Python (optional)
- A machine-readable JSON export
- An exclusion list for the amendments step

Additional v1.0 deliverables:
- Master script wiring all steps in sequence
- Step order flexibility for Steps 3, 4 and 5
- Full instruction file reference documentation
- Test dataset for validating a fresh install
- Contribution guidelines

---

## Version 1.1: Distribution and Packaging

**Target:** Make CRISP installable and citable.

- Bioconda recipe submission for `conda install -c bioconda crisp`
- JOSS (Journal of Open Source Software) paper submission
- GitHub Releases with versioned archives
- Automated test suite for CI
- Snakemake wrapper for HPC workflow integration
- Nextflow support as optional alternative

---

## Version 1.2: LaTeX QC Report and HWE Meta-analysis

**Target:** Publication-ready output and advanced variant QC.

**LaTeX QC Report**

A full compilable LaTeX document generated at pipeline completion,
covering every QC step with formatted tables, embedded figures, and
summary statistics. Designed for direct inclusion in methods sections
or supplementary materials.

```text
LATEX_REPORT    = YES
LATEX_ENGINE    = pdflatex    # Options: pdflatex, xelatex
```

**HWE Meta-analysis**

Rather than a simple p-value threshold, CRISP will compute HWE
statistics stratified by case/control status and ancestry group,
then meta-analyse across strata using Fisher's method. This prevents
removal of genuine association signals that show HWE deviation only
in cases, and prevents spurious failures driven by population
structure.

```text
HWE_META        = YES         # Default: NO
HWE_STRATIFY    = YES         # Stratify by phenotype
HWE_PHENO_FILE  = /path/to/pheno.txt
```

Outputs include a Manhattan-style plot of HWE p-values across
chromosomes and a QQ plot for each stratum.

---

## Version 1.3: Notification System

**Target:** Let the pipeline tell you when it is done.

Email notification via sendmail or Gmail SMTP, phone notification
via Pushover, and Slack webhook support. All optional and configurable
via the instruction file.

```text
NOTIFY_EMAIL    = your@email.com
NOTIFY_PUSHOVER = YES
NOTIFY_SLACK    = https://hooks.slack.com/...
```

---

## Version 1.4: Methodology Preparation Modules

**Target:** Extend CRISP beyond QC into analysis-ready dataset preparation.

Each module prepares a clean CRISP output dataset for a specific
downstream analysis methodology. Developed in consultation with
domain experts.

Planned modules:

| Module | Purpose |
|--------|---------|
| `crisp_prepare_mr.sh` | Mendelian Randomisation: instrument selection, clumping, palindromic SNP handling, exposure and outcome harmonisation |
| `crisp_prepare_prs.sh` | Polygenic Risk Scoring: LD reference panel alignment, score file formatting |
| `crisp_prepare_finemapping.sh` | Fine-mapping: credible set preparation, LD matrix generation |
| `crisp_prepare_coloc.sh` | Colocalization: regional summary stat formatting, MAF alignment |
| `crisp_prepare_gwas.sh` | GWAS: covariates, kinship matrix, population stratification |

---

## Version 1.5: Multi-cohort Mode

**Target:** Process multiple datasets through a single QC run.

Batch processing of multiple datasets using identical QC parameters,
with a combined cross-cohort summary report for consortium-level
comparison. Directly useful for meta-analysis projects where
harmonised QC across cohorts is essential.

```text
COHORT_MODE     = YES
COHORT_LIST     = /path/to/cohort_list.txt
```

---

## Version 1.6: Array Manifest Integration

**Target:** Smarter handling of non-standard probe names.

Support for Illumina array manifest files to rescue non-standard
probe names and assign correct genomic coordinates. Planned to
include rule-based coordinate rescue followed by ML-assisted
resolution for ambiguous cases, with per-array-type model
persistence across runs.

```text
ARRAY_MANIFEST  = /path/to/manifest.csv
USE_ML_RESCUE   = YES
```

---

## Version 2.0: CRISP-py

**Target:** A fully Python-native reimplementation of CRISP.

Separate repository: `ipupko/CRISP-py`

CRISP-py replaces the bash orchestration chunks and R plotting scripts
with pure Python, giving a single-language codebase that is easier to
package, test, and distribute via PyPI. The instruction file format
will remain identical to CRISP so existing users can switch without
changing their configuration.

Key differences from CRISP v1.x:
- Single language: Python throughout
- Installable via pip: `pip install crisp-py`
- Native unit test suite
- Identical instruction file format for seamless migration
- R plotting remains available as an optional backend

---

## Version 2.1: GUI

**Target:** A graphical interface for users who prefer not to edit
plain-text files.

Likely implemented as a Shiny app given the existing R ecosystem,
doubling as an interactive QC report viewer where plots from each
step are rendered dynamically rather than as static PDFs. A desktop
option via PyQt or tkinter is also possible for CRISP-py users.

The instruction file acts as the interface contract throughout --
the GUI simply reads and writes crisp_instructions.txt so the
underlying pipeline never needs to change.

---

## Version 3.0: Cloud Native

**Target:** Native support for major cloud platforms.

AWS, GCP and Azure job submission alongside the existing SLURM
integration. Dataset streaming from cloud storage without requiring
local copies of large files. Cost estimation before run submission.

```text
CLOUD_PROVIDER  = AWS         # Options: AWS, GCP, AZURE, NONE
CLOUD_BUCKET    = s3://my-bucket/data
```

---

## Future Considerations

The following ideas are under consideration with no timeline assigned.

**Containerisation**
Docker and Singularity images for fully reproducible environments
without any manual dependency installation.

**Automated reference panel integration**
Direct download and alignment of 1000 Genomes or gnomAD reference
panels for ancestry inference in PCA.

**Interactive QC dashboard**
A browser-based dashboard showing all QC metrics across steps in
a single view, with the ability to adjust thresholds interactively
and rerun affected steps.

**Citation tracking**
Automatic generation of a methods paragraph citing the correct
versions of PLINK, R, and CRISP used in each run for direct
inclusion in manuscript methods sections.

---

*Last updated: May 2026*
*Developed by Igor Pupko: https://github.com/ipupko/CRISP*
