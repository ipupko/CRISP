<img width="1408" height="768" alt="Firefly_Logo for CRISP (Comprehensive Robust Integrated SNP Processing), a genomic data quali 726701" src="https://github.com/user-attachments/assets/f7312a7e-9584-42c7-b902-4351c3cfa930" />

# CRISP 
### Comprehensive Robust Integrated SNP Processing !!! Under active Development!!!
*Part of the [Compass Genomics](https://github.com/ipupko) suite*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.2.0-green.svg)]()
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20HPC%20%7C%20SLURM-orange.svg)]()
[![Language](https://img.shields.io/badge/language-Bash%20%7C%20R%20%7C%20Python-lightgrey.svg)]()

---

CRISP is a genotype quality control pipeline for large-scale SNP array data. It is built around a single plain-text instruction file, runs on anything from a laptop to a 1000-core HPC cluster, and produces structured reports and publication-ready plots at every step.

Most QC pipelines are a collection of bash scripts someone wrote for their own dataset. CRISP is what happens when you decide to do it properly.

---

## Why CRISP

Genotype QC is not glamorous but it matters enormously. A poorly QC'd dataset introduces systematic biases that no downstream analysis can fully recover from. CRISP was written because the existing options were either too opaque, too rigid, or too tied to a specific cohort's conventions.

A few things that are genuinely different here:

- **One instruction file controls everything.** No Python environments to configure, no Snakemake DAGs to write. If a wet lab collaborator can edit a text file, they can run CRISP.
- **Every step is auditable.** MD5 checksums at the start, structured reports at every step, JSON exports throughout. You can hand the output folder to a reviewer and they can reconstruct exactly what was done.
- **Three call rate modes.** SIMPLE (single threshold), CASCADE (progressive tightening), and CUSTOM (user-defined tiers). Most pipelines only offer SIMPLE.
- **Dual plotting engine.** R and Python produce identical outputs. No more arguments about ggplot2 vs matplotlib.
- **Biobank scale tested.** Explicitly designed and tested at UK Biobank scale -- 487k samples, 784k variants.
- **AI-assisted threshold recommendation coming in v1.3.** CRISP will analyse your missingness distribution and suggest thresholds based on where natural breaks occur, rather than applying generic defaults.

---

## Pipeline steps

| Step | Name | Description |
|------|------|-------------|
| 1 | File validation | MD5 checksums, file signature checks, dataset dimensions |
| 2 | Format conversion | PLINK, PED, VCF, BCF and BGEN all supported |
| 3 | Sample call rate | Per-sample missingness filtering, SIMPLE/CASCADE/CUSTOM modes |
| 4 | Variant call rate | Per-variant missingness, monomorphic removal, MAF filtering |
| 5 | Sex check | F-statistic and Y chromosome count, mismatch and aneuploidy detection |
| 6 | Homozygosity | Heterozygosity outliers, runs of homozygosity, Z-score filtering |
| 7 | Relatedness | IBD (default) or KING at biobank scale |
| 8 | HWE | Hardy-Weinberg equilibrium filtering with optional meta-analysis |
| 9 | PCA | LD pruning, principal components, ancestry visualisation |
| 10 | Amendments | Applies all exclusion lists, produces final clean dataset |
| 11 | Summary report | End-to-end QC summary across all steps |

Steps 3, 4 and 5 are interchangeable via `STEP_ORDER` in the instruction file. Different groups apply these in different orders and CRISP does not enforce a preference.

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| PLINK 1.9 | >= 1.90b | Core QC operations |
| PLINK 2 | >= 2.0 | Format conversion, KING relatedness |
| R | >= 4.0 | Default plotting engine |
| Python | >= 3.8 | Alternative plotting engine, Step 2 orchestration |

R packages: `ggplot2`, `data.table`, `scales`

Python packages (only needed if `PLOT_ENGINE = PYTHON`): `matplotlib`, `pandas`, `numpy`

CRISP checks all tool paths at startup and reports clearly before running anything. Missing optional tools generate warnings, missing required tools cause an immediate exit.

---

## Quick start

```bash
git clone https://github.com/ipupko/CRISP.git
cd CRISP

cp pipeline_instructions.txt my_project.txt
```

Open `my_project.txt` and set at minimum:

```text
INPUT_FORMAT        = PLINK
INPUT_PATH          = /path/to/your/data
OUTPUT_DIR          = ./results
PROJECT_NAME        = my_project
COHORT_POPULATION   = EUR
```

Then run:

```bash
bash crisp_upload_chunk.sh --config my_project.txt
bash crisp_convert.sh      --config my_project.txt
bash crisp_callrate.sh     --config my_project.txt
bash crisp_snprate.sh      --config my_project.txt
bash crisp_sexcheck.sh     --config my_project.txt
# ... and so on
```

Or via the master script:

```bash
bash crisp_master.sh --config my_project.txt
```

On SLURM:

```bash
sbatch crisp_master.sh --config my_project.txt
```

---

## Instruction file

Everything is controlled through `crisp_instructions.txt`. The file ships with sensible defaults - only `INPUT_FORMAT`, `INPUT_PATH` and `COHORT_POPULATION` are compulsory.

```text
# --- INPUT -------------------------------------------
INPUT_FORMAT        = PLINK     # PLINK | PED | VCF | BCF | BGEN
INPUT_PATH          = /data/ukb/ukb_qc_2024
OUTPUT_DIR          = ./results
PROJECT_NAME        = ukb_qc_2024

# --- COHORT (REQUIRED) -------------------------------
COHORT_POPULATION   = EUR       # UNKNOWN | EUR | AFR | EAS | SAS | AMR | MID | CSA
COHORT_SUBPOPULATION= NONE      # e.g. FIN, HUN, PAK, ETH

# --- CALL RATE ----------------------------------------
CALLRATE_MODE       = CASCADE   # SIMPLE | CASCADE | CUSTOM
MIND                = 0.05
GENO                = 0.05

# --- VARIANT QC --------------------------------------
MAF                 = 0.01
REMOVE_MONO         = YES
FILTER_MAF          = YES
HWE                 = 1e-6

# --- SEX CHECK ---------------------------------------
SEXCHECK_EXCLUDE_MISMATCH   = YES
SEXCHECK_EXCLUDE_ANEUPLOIDY = YES

# --- OUTPUT ------------------------------------------
PLOT_ENGINE         = R         # R | PYTHON
EXPORT_FORMAT       = PLINK     # PLINK | VCF | BGEN | GEN | ALL
SCHEDULER           = NONE      # NONE | SLURM
```

See `pipeline_instructions.txt` for the full refernce including tool paths, phasing parameters, downstream tool preparation flags, and advanced options.

---

## Output structure

```
results/
├── ukb_qc_2024_input_md5.txt
├── ukb_qc_2024_input_summary.txt
├── ukb_qc_2024_parsed_config.json
├── logs/
├── step2_converted/
├── step3_callrate/
│   ├── ukb_qc_2024_exclusions_step3.txt
│   └── step3_callrate_cascade_faceted.pdf
├── step4_snprate/
│   ├── ukb_qc_2024_exclusions_step4.txt
│   ├── ukb_qc_2024_exclusions_step4_monomorphic.txt
│   ├── ukb_qc_2024_exclusions_step4_maf.txt
│   └── step4_snprate_diagnostic_summary.pdf
├── step5_sexcheck/
│   ├── ukb_qc_2024_exclusions_sex_mismatch.txt
│   ├── ukb_qc_2024_exclusions_aneuploidies.txt
│   └── step5_sexcheck.pdf
├── step6_homozygosity/
│   └── step6_homozygosity_outliers.pdf
└── final/
    ├── ukb_qc_2024_clean.bed
    ├── ukb_qc_2024_clean.bim
    ├── ukb_qc_2024_clean.fam
    └── ukb_qc_2024_QC_summary.txt
```

---

## Compass Genomics suite

CRISP is part of a broader set of tools for end-to-end genomic data processing:

| Tool | Purpose | Status |
|------|---------|--------|
| **CRISP** | Genotype QC | Active development |
| **CRISP-py** | Python-native reimplementation of CRISP | Planned v2.0 |
| **CAIRN** | Phenotype and metabolomics imputation | Active development |
| **CODA** | Post-analysis QC and result visualisation | Planned post CRISP-py |

All tools share the same instuction file philosophy. If you know how to run CRISP you will know how to run CAIRN.

---

## Roadmap highlights

- **v1.1** -- Bioconda, JOSS paper, Snakemake/Nextflow wrappers, SHAPEIT2 phasing, multi-allelic support
- **v1.3** -- AI-assisted threshold recommendation, automated methods paragraph generation
- **v1.5** -- TOPMed and HRC pre-imputation preparation modules
- **v2.0** -- CRISP-py with native format conversion engine (no PLINK dependency for Step 2)

Full roadmap in [ROADMAP.md](ROADMAP.md).

---

## Citation

If you use CRISP in your research please cite this repository until the JOSS paper is published:

```
Pupko I. (2026) CRISP: Comprehensive Robust Integrated SNP Processing.
Compass Genomics. https://github.com/ipupko/CRISP
```


---

## Contributing

Issues and feature requests are welcome - please open one before submitting a pull request so we can discuss the proposed change first. If you are using CRISP on a dataset and hit something unexpected, a minimal reproducible example is enormously helpful.

---

## License

MIT. See [LICENSE](LICENSE).
