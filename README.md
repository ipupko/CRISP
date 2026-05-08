<img width="1408" height="768" alt="Firefly_Logo for CRISP (Comprehensive Robust Integrated SNP Processing), a genomic data quali 726701" src="https://github.com/user-attachments/assets/f7312a7e-9584-42c7-b902-4351c3cfa930" />



# CRISP
Comprehensive Robust Integrated SNP Processing
!CURRENTLY UNDER ACTIVE DEVELOPMENT!


CRISP is a modular, instruction-file-driven pipeline for quality control and processing of large-scale genotyping data. It was developed with one goal in mind: to make robust, reproducible genotype QC accessible to anyone working with SNP array data, regardless of their computational background.

Whether you are processing a few hundred samples on a local Linux server or running a biobank-scale dataset on an HPC cluster, CRISP handles the heavy lifting through a single plain-text instruction file. You specify what you want. CRISP does the rest.

---

## Table of Contents

- [Background](#background)
- [What CRISP Does](#what-crisp-does)
- [Supported Input Formats](#supported-input-formats)
- [Requirements](#requirements)
- [Getting Started](#getting-started)
- [Output Structure](#output-structure)
- [A Note on Thresholds](#a-note-on-thresholds)
- [HPC and SLURM](#hpc-and-slurm)
- [Questions and Contributions](#questions-and-contributions)

---

## Background

Anyone who has worked with large genotyping datasets knows that quality control is rarely straightforward. Call rate filtering, sex verification, heterozygosity outliers, chromosomal aneuploidies, relatedness -- each step has its own quirks, and getting them right matters enormously for downstream analyses. Mistakes at the QC stage propagate silently into association results, PRS models, and everything that follows.

CRISP grew out of a collection of standalone R scripts developed and refined through real genotyping projects. Those scripts worked well, but running them required manual intervention at every step, careful bookkeeping of intermediate files, and a reasonable amount of bioinformatics experience just to get started. CRISP brings all of those steps together into a single coherent pipeline, with consistent logging, structured outputs, and sensible defaults at every stage.

> If you are new to genotype QC, CRISP is designed to guide you through the process. If you are experienced, it is designed to get out of your way.

---

## What CRISP Does

CRISP processes genotyping data through a series of independent, chainable QC steps. Each step can be enabled or disabled individually. All thresholds have tested defaults and can be overridden in the instruction file.

| Step | Name | Description |
|:----:|------|-------------|
| 1 | File Validation | MD5 checksum generation, input file summary, instruction file parsing |
| 2 | Format Conversion | Optional conversion of PED/MAP, VCF, BCF, or BGEN to BED/BIM/FAM |
| 3 | Sample Call Rate | Sample-level missingness filtering (`--mind`) |
| 4 | Variant Call Rate | Variant-level missingness filtering (`--geno`), tiered reporting |
| 5 | Sex Check | F-statistic and Y-count sex verification, mismatch detection |
| 6 | Aneuploidy Detection | Turner (X0), Klinefelter (XXY), and triple-X (XXX) identification |
| 7 | Heterozygosity and Homozygosity | Z-score outlier detection, runs of homozygosity |
| 8 | Relatedness | IBD-based duplicate and relative detection (KING available as alternative) |
| 9 | Variant QC | MAF filtering, HWE testing, monomorphic SNP removal |
| 10 | PCA | LD pruning, principal component analysis, ancestry visualisation |
| 11 | Amendments | Application of all exclusion lists, clean final dataset production |
| 12 | Summary Report | End-to-end QC summary across all steps |

---

## Supported Input Formats

CRISP accepts the four most common genotyping data formats used in human genetics research:

| Format | Notes |
|--------|-------|
| `BED/BIM/FAM` | Binary PLINK format. Most common for processed array data |
| `PED/MAP` |     Text PLINK format. Automatically converted to BED in Step 2 |
| `VCF / BCF` |   Standard variant call format, including gzipped VCF |
| `BGEN` |        Oxford format, commonly used for UK Biobank data |

All formats are converted to BED/BIM/FAM internally before QC begins.

---

## Requirements

CRISP is written in Python (orchestration), R (QC plots and reports), and Bash (master script). It runs on any Unix/Linux/macOS system. Optional SLURM support is available for HPC environments.

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | >= 1.90b | Core genotype processing |
| [PLINK 2](https://www.cog-genomics.org/plink/2.0/) | >= 2.0 | Format conversion, optional KING relatedness |
| Python | >= 3.8 | Pipeline orchestration |
| R | >= 4.0 | QC plots and reports |

### R Packages

```r
install.packages(c("tidyverse", "data.table", "ggplot2", "scales"))
```

### Python Packages

Only required if `PLOT_ENGINE = PYTHON`:

```bash
pip install matplotlib seaborn pandas
```

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/ipupko/CRISP.git
cd CRISP
```

### 2. Set up your instruction file

```bash
cp pipeline_instructions.txt my_project.txt
```

Open `my_project.txt` and set at minimum:

```text
INPUT_FORMAT    = PLINK
INPUT_PATH      = /path/to/your/data
OUTPUT_DIR      = ./results
PROJECT_NAME    = my_project
```

Everything else has a sensible default. Adjust thresholds as needed for your dataset.

### 3. Run the pipeline

```bash
bash master_pipeline.sh --config my_project.txt
```

That is it. CRISP will work through each enabled step in sequence, writing logs and structured outputs as it goes.

---

## Output Structure

Every run produces a structured output directory:

```
results/
├── my_project_step1_input_summary.txt     # Input file inventory and MD5 checksums
├── my_project_parsed_config.json          # Full record of parameters used
├── step3_sample_callrate/
│   ├── *.irem                             # Excluded sample IDs
│   └── *.imiss                            # Per-sample missingness report
├── step4_variant_callrate/
│   ├── SNP_excluded.txt                   # Variant exclusion counts by tier
│   └── callrate_histogram.pdf             # Missingness distribution plot
├── step5_sex_check/
│   ├── Sex_check.pdf                      # F-statistic vs Y-count plots
│   └── report_sex.mismatch_aneuploidies.txt
├── step6_aneuploidy/
│   └── id_list_aneuploidies.txt           # Samples flagged for exclusion
├── step7_homozygosity/
│   ├── Homozygosity_Runs_vs_Zscore.pdf
│   └── outliers_high_heterozygosity_iids.txt
├── step8_relatedness/
│   ├── relatedness_report.txt
│   └── id_list_related.txt
└── final/
    ├── my_project_clean.bed/.bim/.fam     # Final QC-passing dataset
    └── my_project_QC_summary.txt          # End-to-end summary report
```

---

## A Note on Thresholds

The default thresholds in CRISP reflect commonly used values in the human genetics literature and are appropriate for most standard SNP array datasets. That said, no single set of thresholds is universally correct. The right values depend on your array platform, sample ancestry, sample size, and downstream analysis goals.

If you are unsure where to start, run CRISP with the defaults first and review the diagnostic plots before making decisions. The plots are designed to help you understand your data, not just to flag exclusions automatically.

---

## HPC and SLURM

For users running on a cluster, set `SCHEDULER = SLURM` in the instruction file. CRISP will submit each step as a dependent SLURM job, with memory and time allocations configurable at the top of `master_pipeline.sh`. Logs for each job are written to the `logs/` directory.

---

## Questions and Contributions

CRISP is actively under development. If you encounter an issue, have a suggestion, or want to contribute, please open an [issue](https://github.com/ipupko/CRISP/issues).

Feedback from people actually using the pipeline on real datasets is particularly valuable.

---

<p align="center">
  Developed by <a href="https://github.com/ipupko">Igor Pupko</a>
</p>
