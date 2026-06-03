# CRISP Roadmap
**Comprehensive Robust Integrated SNP Processing**
*Part of the Compass Genomics suite*

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
| `crisp_snprate.sh` | Step 4: Variant call rate, monomorphic, MAF | Under review |
| `crisp_sexcheck.sh` | Step 5: Sex check and aneuploidy detection | Under review |
| `crisp_homozygosity.sh` | Step 6: Heterozygosity and homozygosity | Under review |
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
- Step order flexibility for Steps 3, 4 and 5 via STEP_ORDER parameter
- Full instruction file reference documentation
- Test dataset for validating a fresh install
- Contribution guidelines
- `crisp_check_software.sh` -- startup tool validation script
- Downstream analysis preparation module (`crisp_prepare_output.sh`)

**Software validation at startup**

Before running any QC steps, `crisp_master.sh` calls
`crisp_check_software.sh` which validates every tool path specified
in the instruction file and prints a clean status table. Missing
required tools cause an immediate exit. Missing optional tools
(phasing, liftover) generate warnings only.

```
[CRISP] Software check:
        plink        : OK   (v1.90b6.21)
        plink2       : OK   (v2.00a3.7)
        Rscript      : OK   (v4.3.1)
        python3      : OK   (v3.11.4)
        bcftools     : OK   (v1.17)
        shapeit2     : MISSING (required for phasing at v1.1)
        liftOver     : MISSING (required for build conversion at v1.1)
```

**Downstream analysis preparation**

CRISP prepares output files for the downstream tools specified via
`PREPARE_*` flags in the instruction file.

| Flag | Tool | What CRISP prepares |
|------|------|---------------------|
| PREPARE_REGENIE | REGENIE | Annotation files, phenotype and covariate files |
| PREPARE_MTAG | MTAG | Harmonised summary statistics per trait |
| PREPARE_SCOPA | SCOPA | Covariance matrix and phenotype files |
| PREPARE_SNPTEST | SNPTEST | GEN/SAMPLE files per chromosome |
| PREPARE_SAIGE | SAIGE | Sparse GRM and phenotype files |
| PREPARE_BOLT_LMM | BOLT-LMM | LD score files and phenotype files |
| PREPARE_PLINK2 | PLINK2 | PGEN/PVAR/PSAM format output |
| PREPARE_METAL | METAL | Per-cohort formatted summary statistics |

**Step order flexibility**

Steps 3 (sample call rate), 4 (variant call rate) and 5 (sex check)
are interchangeable. Different research groups apply these in different
orders depending on their QC philosophy. CRISP does not enforce a
fixed order for these three steps.

```text
STEP_ORDER = 3,4,5    # Default: sample call rate, variant QC, sex check
STEP_ORDER = 5,3,4    # Sex check first, then sample call rate, then variant QC
STEP_ORDER = 4,5,3    # Variant QC first, then sex check, then sample call rate
```

Each step resolves its input dynamically from the most recently
produced clean dataset rather than assuming a fixed predecessor.
Steps 6 (homozygosity), 7 (relatedness), 8 (HWE), 9 (PCA),
10 (amendments) and 11 (report) remain fixed in that order.

**Cohort population declaration (compulsory from v1.0)**

Every CRISP project must declare the cohort ancestry in the instruction
file. CRISP records this in all reports from v1.0 onwards. The field
is informational only in v1.0. Population-aware QC adjustments,
reference panel selection and threshold calibration are applied
automatically from v1.3 onwards.

```text
COHORT_POPULATION    = UNKNOWN    # Required. Default: UNKNOWN
COHORT_SUBPOPULATION = NONE       # Optional subgroup
COHORT_ADMIXED       = NO         # YES if mixed ancestry
```

Broad categories: UNKNOWN, EUR, AFR, AMR, EAS, SAS, MID, CSA

Subpopulation codes:
- EUR : FIN, HUN, ITA, SWE, GBR, IBS, TSI, CEU, NFE
- AFR : YRI, LWK, GWD, MSL, ESN, ASW, ACB, ETH, MLI
- SAS : PAK, LKA, BEB, STU, ITU, PJL, GIH
- EAS : CHN, JPT, KHV, CDX, CHS, CHB
- AMR : MXL, PUR, CLM, PEL

**Export step (Step 11b): Output format conversion**

After the final amended clean dataset is produced, CRISP offers
conversion to the user's preferred output format via
`crisp_export.sh`. Controlled by `EXPORT_FORMAT` in the instruction
file.

```text
EXPORT_FORMAT = PLINK    # Options: PLINK, VCF, BGEN, GEN, ALL
EXPORT_BUILD  = GRCh37   # Options: GRCh37, GRCh38 (liftover applied if needed)
```

| Format | Tool | Output |
|--------|------|--------|
| PLINK | PLINK 1.9 | .bed .bim .fam |
| VCF | PLINK 2 | per-chromosome .vcf.gz + .tbi index |
| BGEN | PLINK 2 | per-chromosome .bgen + .sample |
| GEN | PLINK 2 | per-chromosome .gen + .sample (Oxford format for SNPTEST/BGENIE) |
| ALL | all of the above | all formats produced simultaneously |

VCF and BGEN outputs are bgzipped and tabix-indexed automatically.
Per-chromosome split is optional via `EXPORT_SPLIT_CHR = YES`.
GEN format is particularly useful for SNPTEST and BGENIE association
analysis pipelines.

---

## Version 1.1: Distribution and Packaging

**Target:** Make CRISP installable and citable.

- Bioconda recipe submission for `conda install -c bioconda crisp`
- JOSS (Journal of Open Source Software) paper submission
- GitHub Releases with versioned archives
- Automated test suite for CI
- Snakemake wrapper for HPC workflow integration
- Nextflow support as optional alternative
- Multi-allelic variant support in Step 2 via MULTIALLELIC_MODE parameter
- Chromosome splitting for large chromosomes prior to phasing
- SHAPEIT2 phasing integration with chunk stitching via --ligate

**Reference Panel Population**

CRISP uses population-matched reference panels for phasing (SHAPEIT2)
and pre-imputation QC. Using a mismatched reference population
silently degrades phasing accuracy and imputation quality, particularly
for rare variants.

Default is European (EUR). Adjust `REF_POPULATION` to match your
cohort ancestry.

```text
REF_POPULATION = EUR    # Default: EUR
                        # Options: EUR, AFR, AMR, EAS, SAS, CSA, MID, ALL
                        #
                        # EUR : European (default)
                        # AFR : African
                        # AMR : Admixed American
                        # EAS : East Asian
                        # SAS : South Asian (Indian subcontinent)
                        # CSA : Central and South Asian
                        # MID : Middle Eastern
                        # ALL : All populations combined (use with caution)
```

Reference panel paths are set per population. CRISP ships with a
reference panel configuration file (`crisp_ref_panels.txt`) mapping
each population code to the correct panel files on the user's system.

```text
# crisp_ref_panels.txt (example)
EUR_GENETIC_MAP   = /ref/shapeit/genetic_map_hg19_EUR
EUR_REF_PANEL     = /ref/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
AFR_GENETIC_MAP   = /ref/shapeit/genetic_map_hg19_AFR
AFR_REF_PANEL     = /ref/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
EAS_GENETIC_MAP   = /ref/shapeit/genetic_map_hg19_EAS
EAS_REF_PANEL     = /ref/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
```

A warning is printed if `REF_POPULATION` does not match the ancestry
inferred from the PCA step (Step 9). This is a warning not an error
since mixed-ancestry cohorts are common and intentional use of
combined panels is valid.

```text
MULTIALLELIC_MODE = EXCLUDE    # Default: EXCLUDE (current behaviour)
                               # KEEP   : retain multi-allelic sites as-is
                               # SPLIT  : split into biallelic records via PLINK 2
```

---

## Version 1.2: LaTeX QC Report, HWE Meta-analysis and Standalone KING

**Target:** Publication-ready output, advanced variant QC, and improved relatedness estimation.

**Standalone KING (replacing PLINK 2 built-in)**

CRISP v1.0 uses PLINK 2's built-in `--make-king-table` for KING
kinship estimation. From v1.2, this is replaced with the standalone
KING software by Manichaikul et al., which exposes additional
capabilities not available through PLINK 2.

Advantages of standalone KING over PLINK 2 built-in:

| Feature | PLINK 2 built-in | Standalone KING |
|---------|-----------------|----------------|
| Pairwise kinship | Yes | Yes |
| Duplicate detection | Yes | Yes, dedicated --duplicate mode |
| Pedigree reconstruction | No | Yes, --pedigree |
| Family clustering | No | Yes, --cluster |
| Relationship inference | Basic | Full (PO, FS, HS, UN) |
| IBD segments | No | Yes |

The pedigree reconstruction module is particularly valuable for
cohorts with known family structure (e.g. founder populations,
biobanks with multi-generational recruitment) where cross-referencing
reconstructed pedigrees against recorded relationships is a powerful
QC check.

```text
KING_BINARY     = king           # Path to standalone KING binary
KING_CUTOFF     = 0.0884         # Second degree threshold
KING_PEDIGREE   = NO             # YES: reconstruct pedigree structure
KING_CLUSTER    = NO             # YES: group related individuals into families
```

Installation: `conda install -c bioconda king` or download from
https://www.kingrelatedness.com

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

## Version 1.3: AI-Powered Features

**Target:** Make CRISP the first genotype QC pipeline with genuine AI integration.

**Monomorphic SNP reference consultation (Step 4 update)**

Monomorphic SNPs in a cohort are not always truly monomorphic -- they
may be rare or absent in a specific sample but present in the broader
population. Blindly removing them can discard real variants,
particularly in non-European cohorts.

Step 4c is updated to cross-reference flagged monomorphic SNPs against
a reference population database before removal.

```text
MONO_ACTION              = REMOVE   # Default: REMOVE | CONSULT | FLAG_ONLY
MONO_REF_PANEL           = 1KG      # Options: 1KG, GNOMAD, HRC, CUSTOM
MONO_REF_POP             = EUR      # Population to check against
MONO_GLOBAL_MAF_THRESHOLD = 0.001   # Keep if globally above this MAF
```

| Action | Behaviour |
|--------|-----------|
| REMOVE | Current behaviour, remove all monomorphic SNPs |
| CONSULT | Cross-reference reference panel, keep variants with global MAF above threshold |
| FLAG_ONLY | Never remove, annotate in report only |

Supported reference panels:
- 1000 Genomes Project phase 3 (26 populations, 5 super-populations)
- gnomAD v4 (preferred for non-European cohorts)
- HRC r1.1
- CUSTOM (user-supplied allele frequency TSV)

CRISP analyses the missingness, heterozygosity and allele frequency
distributions in your dataset and recommends QC thresholds based on
where natural breaks occur in the data. A model trained on published
QC parameters from large biobanks (UK Biobank, EPIC, deCODE, FinnGen)
provides array-type-aware suggestions rather than generic defaults.

```text
AI_THRESHOLD_RECOMMEND = YES    # Default: NO
AI_ARRAY_DETECT        = YES    # Auto-detect array type from BIM file
```

Example output:
```
[CRISP AI] Array type detected      : UK Biobank Axiom
[CRISP AI] Recommended MIND         : 0.02 (your data: bimodal break at 0.019)
[CRISP AI] Recommended GENO         : 0.02 (array-type default)
[CRISP AI] Recommended MAF          : 0.01 (standard for this array)
[CRISP AI] Confidence               : HIGH
```

**Anomaly detection**

Rather than fixed Z-score cutoffs for heterozygosity and relatedness,
a model flags samples that look unusual relative to the rest of the
cohort. Catches outliers that fixed thresholds miss and reduces
over-exclusion of legitimate edge cases.

```text
AI_ANOMALY_DETECT = YES    # Default: NO
```

**Automated methods paragraph**

After all steps complete, CRISP generates a publication-ready
plain-English methods paragraph with correct citations for PLINK,
R, and CRISP itself. Ready to paste directly into a manuscript
methods section or supplementary materials.

```text
AI_METHODS_TEXT = YES    # Default: NO
```

Example output:
```
Genotype quality control was performed using CRISP v1.0 (Pupko, 2026).
Samples were excluded for excess missingness using a progressive cascade
threshold (mind=0.25, 0.20, 0.10, 0.05), removing 1,025 of 487,409
samples. Sex was verified against reported phenotype using X-chromosome
F-statistics (female: F < 0.2, male: F > 0.8), with 7 samples flagged
for sex mismatch or chromosomal aneuploidy...
```

**QC decision assistant**

An LLM-powered conversational mode where the user can ask questions
about their QC outputs and receive contextual answers grounded in
the actual pipeline results.

```text
crisp_assistant.sh --query "Should I be worried about this heterozygosity plot?"
```

**Relative-informed pseudo-imputation (Step 7b)**

When two related samples are identified in Step 7, the lower quality
sample is normally discarded. However it carries real genotype
information -- particularly in regions where the higher quality sample
has missing calls. Within IBD regions, the two samples share sequence
identical by descent and missing calls in the kept sample can be
filled from the discarded one with high confidence.

This is not population reference panel imputation. It is
haplotype-informed fill-in from a known relative within confirmed
IBD tracts. Confidence is much higher than standard imputation
because the relationship and IBD structure are explicitly known.

```text
RELATIVE_IMPUTE      = NO     # Default: NO | YES
RELATIVE_IMPUTE_IBD  = IBD2   # IBD2: both alleles shared (high confidence)
                               # IBD1: one allele shared (lower confidence)
                               # BOTH: fill IBD2 confidently, IBD1 with flag
RELATIVE_IMPUTE_LOG  = YES    # Log every filled position with source and IBD status
```

Fill confidence by IBD status:

| IBD status | Alleles shared | Fill confidence | Action |
|------------|---------------|----------------|--------|
| IBD2 | Both | Very high | Fill directly, log source |
| IBD1 | One | Moderate | Fill one allele, flag for review |
| IBD0 | None | None | Skip, standard imputation if needed |

Every filled position is recorded in `pseudo_imputation_log.txt`
with source sample IID, IBD status, chromosome, position, and
original vs filled genotype. Filled positions are clearly flagged
in all downstream output files.

Prerequisites: standalone KING (v1.2) for IBD segment detection,
SHAPEIT2 phasing (v1.1) for high-resolution IBD tract boundaries.

This is a novel methodological contribution with no equivalent in
any existing published pipeline. Particularly valuable for:
- Founder populations with long IBD tracts (Finnish, Pakistani, Amish)
- Rare variants poorly covered by population reference panels
- Clinical samples where maximising data completeness matters
- Arrays with poor coverage in specific regions

**AI-assisted sample retention scoring**

When first or second degree relatives are identified, the AI layer
computes a composite QC score across all preceding steps and
recommends which sample to keep with a full explanation.

```text
RELATEDNESS_KEEP = AI_SCORE   # Keep sample with higher composite QC score
```

Composite score weighted across:
- Call rate (Step 3, weight 0.30)
- Heterozygosity Z-score (Step 6, weight 0.25)
- ROH KB (Step 6, weight 0.20)
- Sex check status (Step 5, weight 0.15)
- Aneuploidy status (Step 5, weight 0.10)

```
[CRISP AI] First degree pair: IND234891 vs IND234892
           Composite QC score : 0.94 vs 0.71
           Recommendation     : KEEP IND234891
           Reason             : Higher call rate, normal heterozygosity,
                                normal ROH. IND234892 shows elevated
                                Z-score and ROH approaching thresholds.
```

`RELATEDNESS_KEEP = MANUAL` is also available from v1.0 for cohorts
where automated decisions are not appropriate. MANUAL mode flags
the pair with full QC scores and waits for explicit researcher
instruction before removing any samples.

---

## Version 1.4: Notification System

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

## Version 1.5: Methodology Preparation Modules

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
| `crisp_prepare_topmed.sh` | TOPMed pre-imputation: strand alignment, allele frequency checks, VCF formatting and upload preparation for the TOPMed Imputation Server |

**TOPMed pre-imputation module detail**

Prepares a CRISP-cleaned dataset for submission to the
[TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov).
Steps include:

- Strand alignment against the TOPMed reference panel (hg38)
- Liftover from hg19/GRCh37 if needed via UCSC liftOver or CrossMap
- Allele frequency comparison against TOPMed AF to flag discordant variants
- Removal of ambiguous (A/T, C/G) palindromic SNPs with MAF > 0.4
- Chromosome-split VCF output (one file per chromosome, bgzipped and tabix-indexed)
- Pre-flight check report confirming file readiness before upload

---

## Version 1.6: Multi-cohort Mode

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

## Version 1.7: Array Manifest Integration

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

**Native format conversion engine**

CRISP-py ships its own genotype format converter, removing the PLINK
dependency for Step 2. This is the centrepiece of v2.0 and a genuine
technical differentiator over PLINK-based pipelines.

Supported natively in v2.0:
- PED/MAP to BED/BIM/FAM (full native implementation)
- VCF/VCF.gz to BED/BIM/FAM (via cyvcf2 bindings)
- BGEN to BED/BIM/FAM (via bgen-reader, with dosage-to-hardcall)
- BCF to BED/BIM/FAM (via pysam bindings)

Advantages over PLINK conversion:

| Feature | PLINK | CRISP-py native |
|---------|-------|----------------|
| Transparent encoding report | No | Yes -- every encoding decision logged |
| Per-variant checksum | No | Yes -- JSON manifest with byte layout |
| Streaming large files | No | Yes -- chunk-based, low RAM footprint |
| Integrated QC during conversion | No | Yes -- single pass, no second scan |
| Round-trip verification | No | Yes -- silent back-conversion check |
| In-memory handoff to Step 3 | No | Yes -- zero disk I/O between steps |
| Full audit trail | No | Yes -- lossless JSON record |

```text
CONVERSION_ENGINE = PLINK     # Default in CRISP bash
                              # NATIVE: CRISP-py built-in converter
                              # AUTO: NATIVE where available, PLINK fallback
```

The native converter is not just a dependency removal -- it makes the
conversion step fully auditable and significantly faster on large
datasets by eliminating the intermediate disk write between Steps 2
and 3.

**Native QC operation reimplementations**

CRISP-py reimplements core QC operations natively using NumPy and
Zarr, eliminating file I/O overhead between steps. Genotypes are
stored as an in-memory array (samples x variants, uint8) and passed
directly between QC functions without writing intermediate files.

Operations reimplemented natively in v2.0:

| Step | Operation | Notes |
|------|-----------|-------|
| Step 3 | Missingness per sample (--mind) | NumPy sum over genotype array |
| Step 4 | Missingness per variant (--geno) | NumPy sum over genotype array |
| Step 4 | MAF calculation | Allele count via NumPy |
| Step 4 | Monomorphic detection | Single pass over allele counts |
| Step 5 | X chromosome F-statistic | Computed directly from X genotypes |
| Step 6 | Autosomal F-statistic (--het) | Observed vs expected het ratio |
| Step 6 | Heterozygosity Z-score | NumPy mean and SD |

Operations retained in PLINK 2 in v2.0:

| Step | Operation | Reason |
|------|-----------|--------|
| Step 7 | Relatedness / KING | Standalone KING retained (v1.2+) |
| Step 8 | HWE exact test | Exact p-value at scale, complex algorithm |
| Step 9 | PCA and LD pruning | flashPCA integration, hard to beat |

Estimated performance improvement at UK Biobank scale (500k samples,
800k variants) on a 32-core machine:

| Operation | PLINK time | CRISP-py native | Speedup |
|-----------|-----------|----------------|---------|
| --missing (samples) | ~45s | ~8s | ~5x |
| --missing (variants) | ~45s | ~8s | ~5x |
| --het | ~30s | ~5s | ~6x |
| --check-sex F-stat | ~20s | ~4s | ~5x |
| MAF calculation | ~25s | ~3s | ~8x |

Speedup comes primarily from eliminating disk I/O between steps, not
from faster algorithms. The algorithms are identical to PLINK -- the
results are numerically equivalent and fully validated against PLINK
output as part of the test suite.

```text
COMPUTE_ENGINE = PLINK      # Default in CRISP bash
                            # NATIVE: CRISP-py built-in implementations
                            # AUTO: NATIVE where available, PLINK fallback
```

The key narrative for the JOSS paper:

"CRISP-py reimplements core QC operations natively, eliminating file
I/O overhead between steps and reducing wall-clock time by 5-8x on
biobank-scale datasets. PLINK 2 is retained for computationally
complex operations where its highly optimised implementations remain
the gold standard."

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

**CODA: Post-analysis QC and Result Visualisation**

CODA is a post-analysis results processing toolkit and part of the
Compass Genomics suite. Planned after CRISP-py reaches v1.0.

Repository placeholder: `ipupko/CODA` (to be created after CRISP-py v1.0)

---

## The Compass Genomics Suite

CRISP is part of Compass Genomics, a collection of open-source tools
for end-to-end genomic data processing developed by Igor Pupko.

| Tool | Purpose | Repository |
|------|---------|------------|
| CRISP | Genotype quality control | ipupko/CRISP |
| CRISP-py | Python-native reimplementation of CRISP | ipupko/CRISP-py |
| CAIRN | Phenotype and metabolomics imputation | ipupko/CAIRN |
| CODA | Post-analysis QC and result visualisation | ipupko/CODA |

All tools share the same design philosophy: plain instruction files,
no GUI required, dual R and Python support, HPC-friendly, and
fully auditable outputs.

https://github.com/ipupko

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
*Developed by Igor Pupko | Compass Genomics | https://github.com/ipupko*

---

## Authors and Contributors

**Lead Developer**
Igor Pupko

**Co-authors (tentative)**
Dr Jared Maina
Dr Vincent Pascat

*Co-authorship subject to confirmation. Contributions and authorship
will be formalised prior to the JOSS paper submission at v1.1.*
