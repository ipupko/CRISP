<img width="1408" height="768" alt="Firefly_Logo for CRISP (Comprehensive Robust Integrated SNP Processing), a genomic data quali 726701" src="https://github.com/user-attachments/assets/dda35904-238b-414b-b2b2-9579f60a6256" />

# CRISP Roadmap
**Comprehensive Robust Integrated SNP Processing**
*Part of the Compass Genomics suite*

## Who is CRISP for?

Honestly, anyone who works with genotype data.

If you are an experienced statistical geneticist, CRISP gives you a reproducible, well-documented pipeline with sensible defaults you can override at every step. You know what you are doing and CRISP stays out of your way while keeping a full audit trail of every decision.

If you are newer to genomic QC, the instruction file guides you through what parameters to set and why. The reports tell you what happened at each step. The PDF manual (coming with v1.0) explains the reasoning behind every design decision. You do not need to be an expert to run CRISP correctly.

If you are a novice or expert designing your study CRISP will aid you in your endeavour by showing you examples of studies that were conducted in ther past.

If you are a clinician or translational researcher who needs clean genotype data but does not want to become a bioinformatician, CRISP is designed for you too. One instruction file, one command, one report you can attach to a grant application or a publication. The technical complexity is inside the pipeline, not in front of it.

---



CRISP is a genotype quality control pipeline. It takes raw genomic data (PLINK, VCF, BCF, or BGEN format) and produces a clean, well-documented dataset ready for GWAS analysis.

That probably sounds straightforward. It isn't. The QC step is where most genomic studies quietly go wrong.

Every research group runs QC differently. Different call rate thresholds, different HWE cutoffs, different approaches to relatedness and ancestry. And almost nobody documents their decisions clearly enough for someone else to reproduce them. When results from different studies get combined in meta-analyses, nobody really knows whether the underlying data was processed consistently. Often it wasn't.

CRISP tries to fix that. One pipeline, one instruction file, every decision logged and justified, every result citable.

---

## Why population-aware?

Most QC pipelines were built for European ancestry cohorts and quietly assume that's what you're giving them. Applied to admixed or non-European cohorts, standard QC parameters introduce systematic biases: inflated false positive rates, lost true signals, and ancestry-confounded results that are hard to detect and harder to correct.

CRISP is built from the ground up to handle diverse cohorts properly. The two-pass stratified HWE framework, the Wahlund-aware filtering, and the gnomAD population-specific comparison are all there specifically because the field has been pretending diverse cohorts are fine with European defaults for too long.

---

## Where things stand

CRISP is actively developed and currently at v0.2.0. Eight of eleven steps are complete and tested. Steps 9, 10 and 11 (PCA-based ancestry inference, final QC pass, and report generation) are in progress.

| Step | Description | Status |
|------|-------------|--------|
| 1 | Upload and intake | Complete |
| 2 | Pre-processing and standardisation | Complete |
| 3 | Sample call rate filtering | Complete |
| 4 | Variant call rate filtering | Complete |
| 4.5 | gnomAD population-specific MAF comparison | Complete |
| 5 | Sex check and aneuploidy detection | Complete |
| 6 | PCA Pass 1 (ancestry assignment) | Complete |
| 7 | Heterozygosity / homozygosity QC | Complete |
| 8 | Wahlund-aware HWE filtering | Complete |
| 9 | PCA Pass 2 (final GWAS covariates) | In progress |
| 10 | Final QC pass and output formatting | In progress |
| 11 | Report generation | In progress |

A JOSS paper is in preparation. Once published, every CRISP-processed dataset will be independently citable by referencing the paper and the instruction file used.

---

## How it works

CRISP is driven by a plain text instruction file. You specify your input data, your project parameters, and any cohort-specific settings. CRISP does the rest.

```bash
# Run the full pipeline
bash crisp_run.sh my_cohort.instructions

# Or run individual steps
bash crisp_convert.sh my_cohort.instructions
bash crisp_callrate.sh my_cohort.instructions
bash crisp_hwe.sh my_cohort.instructions 2
```

Every step writes a log, produces a QC report section, and records its decisions to the metadata registry. Nothing happens silently.

---

## What makes it different

**Two-pass PCA architecture.** The classic chicken-and-egg problem in genomic QC is that you need ancestry labels to do proper stratified QC, but you need clean data to get reliable ancestry labels. CRISP solves this with two PCA passes: a lightweight pre-cleaning pass for ancestry assignment, and a post-cleaning pass for final GWAS covariates.

**Wahlund-aware HWE filtering.** Standard HWE testing on a mixed-ancestry cohort will flag real variants as failures because it mistakes population structure for genotyping error. CRISP tests HWE within ancestry strata and combines evidence across strata using Fisher's method. Variants that fail in cases but not controls are flagged and retained rather than dropped. They might be signals, not artefacts.

**Scalable filtering modes.** Fixed p-value HWE thresholds break down at large sample sizes. They remove an ever-growing fraction of variants as N increases, regardless of whether those variants are actually problematic. At cohort sizes above 100,000, CRISP automatically switches to a deviation-based filter that is stable regardless of sample size.

**Population-specific gnomAD comparison.** Step 4.5 compares cohort allele frequencies against gnomAD population-specific frequencies and identifies which reference population best matches the cohort, before PCA runs. Variants that diverge from all gnomAD populations simultaneously are flagged as strong genotyping error candidates.

**Metadata registry.** Every CRISP run writes a structured fingerprint to a local SQLite database: sample size, ancestry composition, relatedness density, QC metrics, and concordance scores. These fingerprints are what COMPASS-AI trains on. The more cohorts that run through CRISP, the smarter the recommendations get. The pipeline and the AI are the same system, not two separate things.

**Dual R/Python plotting.** Every QC plot has both an R/ggplot2 and a Python/matplotlib implementation. You choose which engine to use and CRISP produces equivalent outputs from either.

**LaTeX Report.** Coming as soon as possible. Probably after version 1.2 is running.

**Graphical Intreface and Real-Time data QC parameter alterantion.** Planned for Version 2.

---

## What CRISP is not

CRISP does not perform imputation. It does not run GWAS. It does not do post-GWAS analysis. It is a QC pipeline: it gets your data clean and documented, then hands it to whatever downstream tool you are using.

It also does not make decisions for you. CRISP flags, documents, and recommends. The analyst confirms. Every exclusion is auditable.

---

## Part of the Compass Genomics ecosystem (Upcoming)

CRISP is the foundation of the Aletheia programme, an independent research initiative building infrastructure to move genomic science beyond statistical association toward causal biological understanding.

A detailed instruction manual covering every parameter, every step, and every design decision is in preparation and will be released as a PDF alongside the JOSS paper. If you want to understand not just what CRISP does but why it does it that way, that document is what you want.

Here is what sits above CRISP in the planned ecosystem:

**AVAR (Ancestry Validation and Admixture Refinement)** uses ancestry-informative marker panels to compute continuous per-individual admixture proportions, enabling individual-level population stratification correction that is mathematically more precise than the cluster-based approaches currently used. It is the part of the pipeline that makes CRISP actually work for diverse cohorts rather than just claiming to. In active development. Open source, MIT licence.

**COMPASS-GWAS** runs REGENIE, SAIGE, PLINK2, and BOLT-LMM in parallel on CRISP-cleaned data and scores cross-tool concordance at every significant locus. If four independent tools with different mathematical assumptions all agree on a signal, that signal is more trustworthy than one any single tool finds alone. We are also keeping multi-phenotype models firmly in mind here. The many-to-many relationship between genetic variants and phenotypes carries biological information that one-at-a-time analysis discards, and COMPASS-GWAS is being designed with that in mind from the start. In development.

**COMPASS-AI** is the adaptive intelligence layer. Every cohort that passes through the pipeline contributes an anonymised fingerprint to a performance database. No individual-level data. No genotypes. Just cohort-level summary statistics describing what the data looked like and how different analytical approaches performed on it. Over time the system learns which approaches work best for which cohort types and makes population-specific recommendations. The data safety model is simple: CRISP processes your data on your infrastructure. What leaves is a fingerprint, not the data itself. COMPASS-AI Native (single institution) and COMPASS-AI Global (cross-institution aggregation across verified research groups only) are planned commercial products built on the open source infrastructure.

**COMPASS-NET** is the long-term research component. A genome-wide neural architecture investigation engine designed to learn joint genotype-phenotype relationships rather than testing them one variant at a time. Think of it as what comes after GWAS has told us where to look. Early stage, but it is what the rest of the ecosystem is ultimately building toward.

The goal is a unified framework that connects genetic variation through molecular intermediaries to phenotypic outcome. Not associations. Mechanisms. CRISP is where that starts.

---

## Installation

```bash
git clone https://github.com/ipupko/CRISP.git
cd CRISP

# Dependencies
# PLINK2: https://www.cog-genomics.org/plink/2.0/
# R (4.0 or later) with ggplot2, data.table
# Python (3.8 or later) with pandas, matplotlib, numpy
# bcftools (for VCF/BCF input)
# figlet (for section banners)
```

Full installation guide and dependency checker coming with v1.0.

---

## Contributing

CRISP is open source under the MIT licence and contributions are welcome. If you find a bug, open an issue. If you have a suggestion, open a discussion. If you want to contribute code, open a pull request.

The project is early stage and moving fast. The best way to follow development is to watch the repository and check the issues tab.

---

## Supporting the project

CRISP is built and maintained by a small independent team without institutional backing. If you find it useful and want to support continued development:

- **GitHub Sponsors** (coming soon)
- **Ko-fi** (coming soon)

Infrastructure costs are modest but real. Every contribution helps.

---

## Acknowledgements

The development of CRISP has been assisted by Claude (Anthropic) for architectural specification, code generation, and scientific communication. AI assistance was used throughout the development process.

---

## Licence

MIT. Free to use, modify, and distribute. Always.

---

## Contact

**Dr Igor Pupko**
Founder and Director, Compass Genomics Foundation
igor@compassgx.com
compassgx.com

---

*CRISP is part of the Aletheia programme. Revealing the mechanisms of biology.*
