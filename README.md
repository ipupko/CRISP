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
| `crisp_snprate.sh` | Step 4: Variant call rate, monomorphic, MAF | Complete |
| `crisp_sexcheck.sh` | Step 5: Sex check and aneuploidy detection | Complete |
| `crisp_homozygosity.sh` | Step 6: Heterozygosity and homozygosity | Complete |
| `crisp_relatedness.sh` | Step 7: Relatedness and duplicates | Complete |
| `crisp_hwe.sh` | Step 8: HWE filtering | Under Review |
| `crisp_pca.sh` | Step 9: PCA and ancestry | In Progress |
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

**Tool citations (auto-generated in every run)**

`crisp_report.sh` produces a `CITATIONS.txt` file listing every tool
used in the run with the correct reference. Paste directly into your
methods section.

```
##########################################################################
  CRISP: STEP 11 CITATIONS
  Generated automatically from your run configuration
##########################################################################

CRISP
  Pupko I. CRISP: Comprehensive Robust Integrated SNP Processing.
  Journal of Open Source Software (in preparation). 2026.
  https://github.com/ipupko/CRISP

PLINK 1.9
  Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ.
  Second-generation PLINK: rising to the challenge of larger and richer
  datasets. GigaScience. 2015;4:7.

PLINK 2.0
  Chang CC et al. PLINK 2.0. https://www.cog-genomics.org/plink/2.0/

KING (relatedness estimation)
  Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM.
  Robust relationship inference in genome-wide association studies.
  Bioinformatics. 2010;26(22):2867-2873.

BCFtools / HTSlib
  Danecek P, Bonfield JK, Liddle J, et al.
  Twelve years of SAMtools and BCFtools.
  GigaScience. 2021;10(2):giab008.

R Statistical Computing Environment
  R Core Team. R: A Language and Environment for Statistical Computing.
  R Foundation for Statistical Computing, Vienna, Austria. 2024.

ggplot2
  Wickham H. ggplot2: Elegant Graphics for Data Analysis.
  Springer-Verlag New York, 2016.

Colour-blind safe palette (Okabe-Ito)
  Okabe M, Ito K. Color Universal Design (CUD): How to make figures
  and presentations that are friendly to colorblind people. 2008.
  Wong B. Points of view: Color blindness.
  Nature Methods. 2011;8(6):441.

[Additional citations added based on your run configuration:
 SHAPEIT2, liftOver, REGENIE, SAIGE, BOLT-LMM etc. as applicable]
##########################################################################
```

**Report output format**

The final report from `crisp_report.sh` can be produced in multiple
formats controlled by `REPORT_FORMAT` in the instruction file.

```text
REPORT_FORMAT = PDF        # Default: PDF summary report
                           # PPT: PowerPoint slide deck of key plots
                           # MP4: Animated slide show (ffmpeg required)
                           # ALL: All formats simultaneously
```

PowerPoint and MP4 modes assemble the most informative plots from each
step -- sex check scatter, homozygosity plot, relatedness IBD scatter,
HWE meta-analysis classification, PCA ancestry scatter -- into a
presentation-ready deck. Designed for group meetings and grant
reporting. Each slide includes the key QC numbers from that step.

```text
REPORT_HIGHLIGHT_ONLY = YES  # Only include flagged/anomalous findings
                             # in PPT/MP4 (default: YES)
                             # NO: include all plots from all steps
```

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

## Version 1.1: Distribution, Packaging and Data Security Framework

**Target:** Make CRISP installable, citable, and trusted by institutions.

- Bioconda recipe submission for `conda install -c bioconda crisp`
- JOSS (Journal of Open Source Software) paper submission
- GitHub Releases with versioned archives
- Automated test suite for CI
- Snakemake wrapper for HPC workflow integration
- Nextflow support as optional alternative
- Multi-allelic variant support in Step 2 via MULTIALLELIC_MODE parameter
- Chromosome splitting for large chromosomes prior to phasing
- SHAPEIT2 phasing integration with chunk stitching via --ligate

**Data Security Framework**

CRISP is designed for use with sensitive human genomic data. The
following security principles are built into the architecture from v1.0
and formally documented in v1.1.

Core principles:

- No raw genotype data ever leaves the user's machine under any setting
- Federated metadata sharing is opt-in only, default is NO, always
- Anonymisation occurs locally before any transmission
- One-way hashing of haplotype sequences (SHA-256) makes reverse
  engineering mathematically impossible
- Differential privacy applied to aggregate statistics -- individual
  contributions cannot be identified from the aggregate model
- Deletion on request with cryptographic confirmation of removal
- All anonymisation code is open source and independently auditable

Required documentation shipped with v1.1:

```
docs/data_governance.md     What is collected, what is not, who has
                            access, retention periods, deletion policy
docs/security.md            Technical security measures, hashing
                            approach, differential privacy parameters
docs/consent_template.md    Template for institutional ethics boards
                            covering CRISP_AI_SHARE opt-in
```

Independent security audit scheduled before v1.5 (federated learning)
ships. Audit report published openly on GitHub.

```text
CRISP_AI_SHARE = NO         # Default: NO. Never changes without
                            # explicit user action.
                            # AGGREGATE: share anonymised QC fingerprint
```

This framework is designed to satisfy the requirements of:
- UK Biobank data access agreements
- GDPR (EU/UK)
- HIPAA (US clinical contexts)
- Institutional ethics board review at major genomics centres

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

## CRISP-FILL: Population-Aware Pseudo-Imputation (Experimental)

**Status:** Optional, experimental. Not enabled by default.
**Prerequisites:** Step 7 relatedness, Step 9 PCA, SHAPEIT2 phasing (v1.1)

CRISP-FILL is a novel experimental module that fills missing genotype
regions using population-aware haplotype matching and IBD-confirmed
relative information. It is the first systematic implementation of
this approach in any published genotype QC pipeline.

**Workflow:**

```
Access CRISP (GitHub or Cloud)
  |
  +--> Step 9 PCA runs
  |      Population inferred from PCA clustering against reference panel
  |      Epidemiological context loaded:
  |        founder population, consanguinity rate,
  |        known rare variant burden for this ancestry
  |
  +--> Step 7 IBD segments identified from relatedness step
  |      Related pairs and their confirmed IBD regions known
  |
  +--> Missing chunks located in kept samples
  |
  +--> Fill sources identified (two independent lines of evidence):
  |      (a) Discarded relative within confirmed IBD region
  |          IBD2: both alleles shared, very high confidence
  |          IBD1: one allele shared, moderate confidence, flagged
  |      (b) Population-matched haplotypes from CRISP-REF (opt-in)
  |          Only used when no relative fill available
  |          Matched by PCA-inferred ancestry
  |
  +--> Gaps filled with per-position confidence score
  |      All fills logged in pseudo_imputation_log.txt
  |      Filled positions flagged in all downstream output
  |
  +--> Epidemiological metadata used to weight fill confidence
         Finnish cohort: long IBD tracts expected, high confidence fills
         Pakistani cohort: high consanguinity, rich relative fill network
         Admixed cohort: ancestry-matched CRISP-REF fills preferred
```

**Why two independent lines of evidence matter:**

Standard imputation uses population haplotypes only. CRISP-FILL adds
a second confirmation -- the filled region must also be in a confirmed
IBD segment with the source relative. Population match plus IBD
confirmation together make spurious fills extremely unlikely.

**Instruction file parameters:**

```text
CRISP_FILL           = NO      # Default: NO | YES (experimental)
CRISP_FILL_IBD       = IBD2    # IBD2: both alleles | IBD1: one allele | BOTH
CRISP_FILL_REF       = NO      # Use CRISP-REF population haplotypes as fallback
CRISP_FILL_LOG       = YES     # Log every filled position
CRISP_FILL_MIN_CONF  = 0.90    # Minimum confidence score to apply fill
```

**Output:**

```
pseudo_imputation_log.txt
  IID, CHR, POS, ORIGINAL_GT, FILLED_GT, SOURCE_IID, IBD_STATUS,
  CONFIDENCE, ANCESTRY_MATCH, FILL_METHOD
```

**Important caveats:**

- Filled positions are not directly measured genotypes
- All filled positions must be clearly flagged in downstream analysis
- Fills outside confirmed IBD regions are not performed
- This module is experimental -- validate against WGS truth data
  before using in published research
- A methods note is required when CRISP-FILL outputs are used in
  association analysis

**The methodological contribution:**

No existing published pipeline performs population-aware, IBD-confirmed
pseudo-imputation systematically. This is a genuine novel contribution
suitable for its own methods paper, separate from the main CRISP JOSS
paper.

**Advanced layers (planned extensions to CRISP-FILL):**

*Layer 1: Metadata-guided variant quality cross-referencing*

Every opt-in CRISP run contributes per-variant QC statistics indexed
by genomic coordinate -- MAF, HWE p-value, call rate, array type,
population. Over hundreds of runs, CRISP builds a per-variant quality
profile. When filling a missing region, the model asks: across all
similar cohorts that processed this locus, what was the typical QC
profile? A variant that consistently passes QC across 500 EUR runs is a
much stronger fill candidate than one with variable behaviour. Genomic
coordinates are the join key across all contributed metadata.

```text
CRISP_FILL_META_QC  = YES   # Cross-reference variant quality from
                            # accumulated run metadata before filling
```

*Layer 2: Suggestive population labelling from known markers*

Certain variants are population-specific fingerprints -- Finnish founder
mutations, South Asian-specific haplotypes, African-specific alleles.
CRISP cross-references detected variants against a curated database of
population-informative markers to assign a probabilistic ancestry label
per sample, refining the PCA clustering and improving fill matching.

A sample labelled as likely Finnish gets Finnish-specific haplotype
fills. A sample with ambiguous PCA position but strong founder marker
signal gets labelled accordingly with a confidence score.

```text
CRISP_FILL_MARKER_DB = YES  # Use population marker database for
                            # suggestive ancestry labelling
CRISP_FILL_MIN_LABEL_CONF = 0.80  # Minimum confidence for label use
```

*Layer 3: Sequence-level learning*

When WGS truth data is available for a subset of samples, CRISP-FILL
can compare its fills against ground truth and update its confidence
model. Fills confirmed by WGS increase the weight of that fill pattern.
Fills that disagree reduce it. Over time CRISP learns which genomic
regions in which populations filled by which IBD status are reliably
correct -- building an empirically grounded confidence model from
real outcomes rather than theoretical priors.

```text
CRISP_FILL_WGS_TRUTH = /path/to/wgs_vcf  # Optional WGS truth for
                                          # confidence model training
```

*Layer 4: Epigenetic integration*

Methylation state at a locus correlates with nearby genotype. When a
sample has RRBS or WGBS methylation data alongside array data,
methylation patterns in a missing region provide independent evidence
about the likely genotype. Combined with IBD confirmation and population
matching, the three-way evidence substantially improves fill reliability.

Telomere length as a biological age proxy (already planned in CRISP)
is a related epigenetic signal -- TL correlates with specific genomic
regions under selection and can provide additional context for fills
near known ageing loci.

```text
CRISP_FILL_EPIGENETIC = NO            # Default: NO
CRISP_FILL_METHYLATION = /path/to/methylation.bed
CRISP_FILL_TL_COVARIATE = YES         # Use telomere length as prior
```

**Evidence hierarchy for fill decisions:**

| Evidence | Weight | Notes |
|----------|--------|-------|
| IBD2 confirmation | Very high | Both alleles identical by descent |
| IBD1 confirmation | Moderate | One allele shared, other inferred |
| Population match (PCA) | High | Same ancestry cluster |
| Marker-based label | Moderate | Founder marker confirmation |
| Metadata QC profile | Moderate | Consistent QC across similar runs |
| Methylation pattern | Supplementary | Correlated, not causal |
| TL biological age | Supplementary | Regional selection signal |

**The long-term vision:**

CRISP-FILL becomes a self-improving imputation engine that learns from
every WGS-validated fill, incorporates population-specific knowledge
from the accumulated metadata, and integrates multi-modal evidence
including epigenetic signals. The confidence model is retrained
continuously as new evidence accumulates.

This is not a replacement for standard imputation panels. It is a
complementary approach that works where standard panels fail -- rare
variants, under-represented populations, regions with poor probe
coverage -- and provides a second independent line of evidence where
standard imputation succeeds.

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

**CRISP-REF: Community Haplotype Reference Resource**

Every time CRISP performs relative-informed pseudo-imputation it
generates high-confidence genotype fills within IBD tracts between
known relatives. These fills are not inferred from population
patterns -- they are grounded in confirmed biological relatedness.
Aggregated across hundreds of cohorts, stripped of all identifiers,
they represent a novel source of haplotype information that current
imputation panels have never had access to.

CRISP-REF is a community-built reference haplotype resource derived
from relatedness-confirmed IBD segments across globally diverse
cohorts. Entirely anonymous, entirely open, entirely novel.

What CRISP-REF captures that existing panels do not:
- Verified IBD segments rather than inferred population haplotypes
- Cross-cohort diversity from cohorts currently under-represented
  in TOPMed, HRC and 1000 Genomes
- Rare variant context confirmed by relatedness rather than frequency
- Haplotype phase information carried by IBD2 fills

The metadata contributed per run is fully anonymised:

```text
crisp_ibd_haplotype_v1 {
    chromosome         : 3
    start_bp           : 45000000
    end_bp             : 46200000
    ibd_status         : IBD2
    population         : EUR_FIN
    n_confirming_pairs : 847
    fill_confidence    : 0.991
    haplotype_hash     : sha256_of_sequence
    # No individual IDs, no sample identifiers,
    # no cohort names, no phenotype data
}
```

The feedback loop:

```
CRISP processes diverse cohorts worldwide
  --> identifies IBD2 pairs, generates high-confidence fills
  --> strips all identifiers, contributes to CRISP-REF
  --> CRISP-REF grows with global diversity
  --> imputation panels query CRISP-REF
  --> better imputation for rare variants and under-represented populations
  --> which improves CRISP's own pseudo-imputation in future runs
```

Why this matters for health equity: the single biggest problem in
genomics right now is that imputation panels work brilliantly for
Europeans and poorly for everyone else. CRISP-REF would naturally
improve with diversity because the cohorts most likely to have
consanguinity and long IBD tracts -- Pakistani, Finnish, Ashkenazi
Jewish, and African founder populations -- would contribute
disproportionately rich IBD metadata.

Prerequisites: relative-informed pseudo-imputation (v1.3), differential
privacy implementation for haplotype hashes, community governance
framework, collaboration with TOPMed or HRC panel maintainers.

```text
CRISP_REF_CONTRIBUTE = NO    # Default: NO
                              # YES: contribute anonymised haplotype
                              # metadata to CRISP-REF after each run
```

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
| CRISP-REF | Community IBD haplotype reference resource | ipupko/CRISP-REF |
| CRISP-RNA | Bulk and single-cell RNA-seq QC | ipupko/CRISP-RNA |
| COMPASS-GWAS | Multi-tool GWAS orchestration | ipupko/COMPASS-GWAS |
| CAIRN | Phenotype and metabolomics imputation | ipupko/CAIRN |
| CODA | Post-analysis QC and result visualisation | ipupko/CODA |

All tools share the same design philosophy: plain instruction files,
no GUI required, dual R and Python support, HPC-friendly, and
fully auditable outputs.

https://github.com/ipupko

**The end-to-end Compass Genomics pipeline:**

```
Raw genotype data
  |
  +--> CRISP (genotype QC)
  |      Steps 1-11: validation, conversion, call rate,
  |      sex check, homozygosity, relatedness, HWE, PCA,
  |      amendments, report
  |      Optional: CRISP-FILL pseudo-imputation
  |
  +--> CAIRN (phenotype and metabolomics imputation)
  |      Missing phenotype imputation
  |      NMR metabolomics QC and imputation
  |      CAIRN-imputed phenotypes fed directly to COMPASS-GWAS
  |
  +--> COMPASS-GWAS (multi-tool GWAS orchestration)
  |      Parallel REGENIE, SAIGE, PLINK2, BOLT-LMM, SNPTEST
  |      Harmonised summary statistics across all tools
  |      Cross-tool concordance check
  |      FUMA annotation (automatic API submission)
  |      PowerPoint / MP4 results presentation
  |
  +--> CODA (results visualisation and interpretation)
         Manhattan plots, QQ plots, locus zoom
         Cross-phenotype comparison
         LD structure visualisation

One instruction file per tool. Data flows automatically between steps.
```

---

**COMPASS-GWAS: Multi-tool GWAS Orchestration**

The single biggest friction point in GWAS analysis is that different
tools require different input formats, different covariate
specifications, and produce different output formats. Running REGENIE,
SAIGE and PLINK2 on the same dataset requires setting up three
separate analysis environments, submitting three sets of jobs, and
then manually harmonising the results.

COMPASS-GWAS removes that friction entirely.

```text
GWAS_TOOLS     = REGENIE,SAIGE,PLINK2   # Run all three in parallel
GWAS_PHENO     = /path/to/pheno.txt
GWAS_COVARIATES= /path/to/covariates.txt
GWAS_THREADS   = 16
FUMA_SUBMIT    = YES                    # Auto-submit to FUMA API
REPORT_FORMAT  = PPT                    # PowerPoint results deck
```

Steps COMPASS-GWAS performs:

1. Format conversion per tool from CRISP clean output
2. Parallel HPC job submission (SLURM array or local)
3. Job monitoring with email/Slack notification on completion
4. Collect and harmonise summary statistics to a standard format
5. Cross-tool concordance check -- hits appearing in only one
   tool are flagged for manual review
6. Meta-analysis across tools if requested
7. FUMA API submission and result retrieval
8. CODA visualisation of top hits
9. PowerPoint / MP4 assembly of key results

**Why cross-tool concordance matters:**

REGENIE, SAIGE and BOLT-LMM make different assumptions about
population structure and genetic architecture. A hit that appears
robustly across all three tools is far more credible than one that
appears in only one. Running all three and comparing is best practice
but almost nobody does it because the overhead is enormous.
COMPASS-GWAS makes cross-tool comparison the default.

**FUMA integration:**

FUMA (Functional Mapping and Annotation of Genome-Wide Association
Studies) maps GWAS hits to genes, pathways, eQTLs and tissue
expression. COMPASS-GWAS submits summary statistics to the FUMA API
automatically and pulls back annotated results for inclusion in
the CODA report and PowerPoint deck.

```text
FUMA_SUBMIT       = YES
FUMA_API_KEY      = your_fuma_api_key
FUMA_MAGMA        = YES   # Run MAGMA gene-set analysis
FUMA_EQTL         = YES   # Cross-reference GTEx eQTLs
FUMA_TISSUE       = Brain,Blood,Liver  # Tissue expression priority
```

**CAIRN connection:**

CAIRN-imputed phenotypes feed directly into COMPASS-GWAS via a
standard phenotype file format shared across both tools. Missing
phenotype data imputed by CAIRN is flagged in the COMPASS-GWAS
report so downstream interpretation accounts for imputed vs
directly measured phenotypes.

**Cross-tool concordance and meta-tool analysis**

When two or more tools are run on the same dataset, COMPASS-GWAS
computes a full concordance report per run:

```
[COMPASS-GWAS] Cross-tool concordance report:

  REGENIE vs SAIGE
  Genome-wide significant hits (p < 5e-8):
    REGENIE only    : 3
    SAIGE only      : 1
    Both tools      : 47    Concordance: 94.0%

  Direction of effect (suggestive hits p < 1e-5):
    Agreement       : 98.2%
    Disagreement    : 1.8%  (see concordance_discordant.txt)

  Lambda:
    REGENIE : 1.021   SAIGE : 1.019

  Novel hits per tool:
    REGENIE: rs1234567 chr3:45123456 p=2.1e-9 (not in SAIGE)
    SAIGE  : rs9876543 chr7:12345678 p=4.8e-9 (not in REGENIE)

  Recommendation (SAS/PAK n=12,847 CC 1:4):
    SAIGE preferred for this cohort type.
    Evidence from 847 similar runs in COMPASS-GWAS database.
```

Variants significant in only one tool are flagged for manual review.
Variants concordant across all tools are marked as high-confidence.
A meta-analysis across tools is available via Fisher's method on
the combined p-values -- the same approach used in CRISP HWE
meta-analysis, applied to association results.

**Population-tool performance matching (AI layer)**

Over hundreds of opt-in COMPASS-GWAS runs, concordance data
accumulates against a growing knowledge base. The AI layer learns
which tool performs best for which cohort profile -- population,
sample size, case-control ratio, phenotype distribution.

This does not exist anywhere in the literature as a systematic
resource. COMPASS-GWAS generates it as a byproduct of use.

Example output from the recommendation engine:

```
Cohort: SAS (Pakistani), n=4,200, case-control 1:10
  Recommended : SAIGE
  Second      : REGENIE
  Evidence    : 847 similar runs, concordance with WGS: 0.94
  Avoid       : BOLT-LMM (not designed for unbalanced CC)

Cohort: EUR (Finnish), n=450,000, quantitative trait
  Recommended : REGENIE (fastest, best calibrated at biobank scale)
  Second      : BOLT-LMM (good for quantitative traits)
  Evidence    : 2,341 similar runs

Cohort: AFR (Yoruba), n=8,000, admixed
  Recommended : SAIGE (handles structure best)
  Warning     : All tools show lambda inflation in this cohort type.
                Verify PCA ancestry components are adequate.
  Evidence    : 203 similar runs (limited -- contribute your run)
```

**The scientific output of the performance database:**

The accumulated performance data is itself a publishable finding.
A systematic cross-tool evaluation across thousands of real-world
cohorts -- stratified by population, sample size, phenotype type --
is a paper that does not yet exist. COMPASS-GWAS generates it
automatically as users opt in.

*"Systematic cross-tool evaluation of GWAS methods across diverse
cohorts reveals population-specific performance differences"*

That is a Nature Methods or AJHG paper generated as a byproduct
of building useful infrastructure.

**The feedback loop:**

```
Users run COMPASS-GWAS
  --> concordance data contributed (opt-in)
  --> population-tool performance database grows
  --> recommendation engine improves
  --> future users get better tool recommendations
  --> those runs contribute more data
  --> repeat
```

The same compounding logic as CRISP federated threshold learning
and CRISP-REF. The more it is used, the better it gets.
The better it gets, the more it is used.

Repository: `ipupko/COMPASS-GWAS` (to be created after CRISP v1.0)

---

## COMPASS-AI: Learned GWAS Engine (Long-term Vision)

**Status:** Long-term research goal. Requires COMPASS-GWAS performance
database at scale as prerequisite. No timeline assigned.

**The problem COMPASS-AI solves:**

Every existing GWAS tool applies a fixed statistical algorithm uniformly
across all variants, all populations, all cohort configurations. SAIGE
uses saddlepoint approximation. REGENIE uses whole-genome regression
with block-wise inversion. BOLT-LMM uses a Bayesian mixture model.
Each was designed with specific assumptions that hold well in some
contexts and poorly in others.

No tool knows when its own assumptions are being violated.
No tool adapts to the data in front of it.
No tool learns from previous runs on similar data.

COMPASS-AI does all three.

**The three levels of COMPASS-AI:**

*Level 1: Empirical performance mapping*

Which tool works best for which cohort profile. Generated directly
from the COMPASS-GWAS concordance database. Population, sample size,
case-control ratio, phenotype distribution, array type -- all mapped
to tool performance. This level is achievable as soon as the database
reaches sufficient scale (estimated 500+ diverse cohorts).

*Level 2: Mechanistic failure mode understanding*

When two tools disagree on a variant, why? Discordant hits plus
cohort metadata plus variant characteristics plus tool parameters
form a training dataset for learning which algorithmic choices cause
which failure patterns in which contexts.

Output is not just "tool A is better here" but "tool A is better
here because the saddlepoint approximation underestimates significance
at MAF < 0.02 in admixed cohorts with case-control imbalance > 1:8."
Mechanistic understanding derived from empirical observation at scale.

This is not available in any current tool or any current paper.

*Level 3: Learned combination and novel tool design*

A meta-learner that combines summary statistics from multiple tools
with learned, population-aware weights -- not fixed Fisher's method
but an adaptive weighting that knows, for this population and this
variant frequency class, how much to trust each tool's p-value.

The learned weights reflect everything the performance database has
accumulated. A variant in a Pakistani founder population at MAF 0.008
in a 1:12 case-control study gets weights derived from thousands of
similar runs across similar cohorts. The combination is demonstrably
better calibrated than any individual tool.

Eventually the meta-learner is replaced by a native learned mixed
model that has internalised what the performance database taught --
a GWAS engine that does not choose between existing approaches but
has learned a better one from observing all of them.

**Why this has not been done before:**

The training data does not exist. Cross-tool, cross-cohort,
population-stratified concordance data at the scale needed to train
something meaningful requires infrastructure that generates it
systematically as a byproduct of routine use. That infrastructure
is COMPASS-GWAS. COMPASS-GWAS requires clean data from CRISP.

The prerequisite chain:

```
CRISP (now)
  --> COMPASS-GWAS (after CRISP v1.0)
  --> Performance database (accumulates with adoption)
  --> COMPASS-AI Level 1 (empirical mapping, ~500 cohorts)
  --> COMPASS-AI Level 2 (failure mode analysis, ~2,000 cohorts)
  --> COMPASS-AI Level 3 (learned engine, ~5,000+ cohorts)
```

Each step is the foundation for the next.
None is possible without CRISP.

**The publication trajectory:**

```
CRISP JOSS paper          : Methods infrastructure
COMPASS-GWAS paper        : Multi-tool orchestration and concordance
Performance database paper : Population-stratified tool evaluation
                            (Nature Methods)
COMPASS-AI Level 1 paper  : Empirical tool selection framework
                            (AJHG or PLOS Genetics)
COMPASS-AI Level 2 paper  : Mechanistic failure mode analysis
                            (Nature Methods or Genome Research)
COMPASS-AI Level 3 paper  : Learned GWAS engine
                            (Nature or Nature Genetics)
```

Six papers from one coherent research programme. Each independently
strong. Together a body of work that reframes how the field does
association analysis.

**Novelty assessment:**

Tool comparison papers exist but are static, single-cohort, and do
not update. Ensemble methods exist in clinical prediction but not in
GWAS mixed models. Federated learning in genomics is nascent and
focused on privacy, not performance learning. AutoML frameworks
exist but have no concept of genomic architecture or population
structure.

The combination of live updating, population stratification,
mechanistic failure analysis, and learned combination trained on
infrastructure-generated data is genuinely unprecedented.

The reason it has not been done is that the data infrastructure
to do it has not existed. That infrastructure is being built now.

Repository: `ipupko/COMPASS-AI` (long-term, after COMPASS-GWAS)

---

*"It is not even cutting edge -- it is bleeding edge."*
*-- Tony Stark (and Igor Pupko, re: Compass Genomics, 2026)*

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

## Version 3.0: CRISP-RNA

**Target:** Standalone RNA-seq QC tool under the Compass Genomics suite.
Independent software with its own repository, JOSS paper, and release cycle.

CRISP-RNA applies the same design philosophy as CRISP -- plain instruction
files, dual R and Python support, fully auditable outputs -- to bulk RNA-seq
and single-cell RNA-seq data.

**Bulk RNA-seq QC modules:**

| Module | What it checks |
|--------|---------------|
| Mapping QC | Alignment rate, uniquely mapped reads, rRNA contamination |
| Sample QC | Read depth, RNA integrity, outlier detection via PCA |
| Feature QC | Low count genes, zero-inflation, mitochondrial proportion |
| Sex check | XIST expression vs Y chromosome gene expression |
| Batch effects | Detection and correction, surrogate variable analysis |
| GC and length bias | Library prep artefacts affecting quantification |
| Strand specificity | Library protocol verification |
| Duplication rate | PCR duplicate assessment |

**Single-cell RNA-seq QC modules (CRISP-SC):**

| Module | What it checks |
|--------|---------------|
| Doublet detection | Cell multiplets from droplet capture |
| Ambient RNA | Contamination from lysed cells |
| Cell-level QC | UMI count, gene count, mitochondrial fraction |
| Gene-level QC | Zero-inflation, highly variable genes |
| Batch integration | Cross-sample normalisation |

**Cross-modal validation with CRISP:**

When a cohort has both genotype array and RNA-seq data, CRISP and
CRISP-RNA cross-validate each other:

- RNA-seq sex check vs genotype sex check -- discordance flags sample swaps
- RNA-seq inferred genotype at common SNPs vs array genotype -- concordance
  check catches contamination and mislabelling
- Relatedness inferred from RNA-seq vs from genotype array -- catches
  incorrectly paired samples

This cross-modal QC is not available in any existing published pipeline
and represents a genuine methodological contribution for the CRISP-RNA
JOSS paper.

```text
CRISP_RNA_MODE    = BULK       # BULK | SINGLECELL
CROSS_MODAL_QC    = NO         # YES: cross-validate with CRISP genotype output
CRISP_OUTPUT_DIR  = ./results  # Path to CRISP genotype QC output for comparison
```

Repository: `ipupko/CRISP-RNA` (to be created at v3.0 development start)

---

## Authors and Contributors

**Lead Developer**
Igor Pupko

**Co-authors (tentative)**
Dr Jared Maina
Dr Vincent Pascat

*Co-authorship subject to confirmation. Contributions and authorship
will be formalised prior to the JOSS paper submission at v1.1.*
