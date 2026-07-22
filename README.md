# HSP and HSF Gene Family Analysis & RNA-seq Pipeline (`hsp_hsf`)

![R Version](https://img.shields.io/badge/R-4.5.1-blue?logo=r)
![Python Version](https://img.shields.io/badge/Python-3.12.12-blue?logo=python)
![Quarto Version](https://img.shields.io/badge/Quarto-1.9.37-blue)
![Mamba Version](https://img.shields.io/badge/Mamba-2.1.1-green)
![License](https://img.shields.io/badge/License-Internal-orange)

A reproducible bioinformatics pipeline designed for identifying, validating, and performing post-quantification quality control (QC) for Heat Shock Proteins (HSP) and Heat Shock Transcription Factors (HSF) across plant species under drought and elevated CO₂ conditions.

## 📑 Table of Contents
- [Target Species](#-target-species)
- [Prerequisites & Environment](#-prerequisites--environment)
- [Key Components & Workflow](#-key-components--workflow)
- [Repository Structure](#-repository-structure)
- [Testing & Verification](#-testing--verification)
- [Running the Exploratory QC Pipeline](#-running-the-exploratory-qc-pipeline)
- [Pipeline Outputs](#-pipeline-outputs)
- [License & Attribution](#-license--attribution)

---

## 🌿 Target Species

This repository covers comparative analysis across four key crop and model plant species:
- ***Setaria viridis*** (Green foxtail) — BioProject `PRJNA1496374`
- ***Sorghum bicolor*** (Sorghum) — BioProject `PRJNA742236`
- ***Triticum aestivum*** (Wheat) — BioProjects `PRJNA1249400` & `PRJNA793265`
- ***Glycine max*** (Soybean) — BioProjects `PRJNA295411` & `PRJNA1033144`

---

## ⚙️ Prerequisites & Environment

This pipeline requires specific tools which may be available via your HPC's `module` system or managed via `mamba`/`conda` environments.

### Core Dependencies:
- **R** (>= 4.5.1): Includes packages like `DESeq2`, `tximport`, `ggplot2`.
- **Python** (>= 3.12.12): For data validation and CLI wrappers.
- **Quarto** (>= 1.9.37): For generating interactive HTML reports.
- **Salmon**: For transcript-level quantification.
- **FastQC**: For quality control of raw reads (usually loaded via HPC module).

---

## 🛠️ Key Components & Workflow

### 1. Reference Identification & Domain Validation
- **HMMER & InterProScan Filtering**: Identifies candidate HSP and HSF proteins using Pfam domain profiles (`Pfam-A.hmm`).
- **Isoform Selection**: Selects representative protein isoforms using canonical transcript mappings and sequence length heuristics.

### 2. RNA-seq Quantification
- **Salmon Pipeline**: Automated quantification of transcript-level and gene-level abundances across 58 RNA-seq libraries.
- **Transcript-to-Gene Aggregation**: Uses `tximport` in R to generate unified gene count, TPM, and effective length matrices.

### 3. Exploratory Quality Control & Experimental Design
- **Experimental Design Normalization**: Classifies complex multi-factor conditions (`co2`, `water`, `dehydration_stage`).
- **DESeq2 Blind VST Normalization**: Unsupervised variance-stabilizing transformation for library diagnostics, sample Spearman correlation, Euclidean sample distance, and PCA plots.
- **Robust Advisory Flags**: Emits 3-MAD distributional flags without automated sample removal.
- **Reproducible Quarto Reporting**: Generates interactive HTML reports (`reports/rnaseq_exploratory_qc.html`) and paired PNG/SVG figures.

---

## 📂 Repository Structure

```text
hsp_hsf/
├── config/                  # Sample manifests and experimental designs
│   ├── rnaseq_samples.tsv
│   └── rnaseq_design.tsv
├── scripts/                 # Core analysis and CLI scripts (R & Python)
│   ├── build_rnaseq_design.R
│   ├── import_salmon_tximport.R
│   ├── rnaseq_exploratory_qc_lib.R
│   ├── run_rnaseq_exploratory_qc.R
│   └── validate_salmon_quant.py
├── reports/                 # Quarto presentation sources and compiled HTML reports
│   ├── rnaseq_exploratory_qc.qmd
│   └── rnaseq_exploratory_qc.html
├── results/                 # Versioned outputs and matrices
│   ├── rnaseq/matrices/final/   # Gene/transcript count and TPM matrices
│   └── rnaseq/exploratory_qc/   # QC tables, figures, contrasts, session info
├── tests/                   # Automated unit & integration tests (Python unittest & R)
└── docs/                    # Technical specifications and architectural design docs
```

---

## 🧪 Testing & Verification

The repository enforces strict test suites covering data alignment, schema validation, and staging rollback:

```bash
# Run Python unit tests
python3 -m unittest discover -s tests -v

# Run R test suites
Rscript tests/test_import_salmon_tximport.R
Rscript tests/test_build_rnaseq_design.R
Rscript tests/test_rnaseq_exploratory_qc.R
```

---

## 🚀 Running the Exploratory QC Pipeline

To regenerate the experimental design and run the exploratory QC workflow:

```bash
# 1. Build the normalized design table
Rscript scripts/build_rnaseq_design.R \
  config/rnaseq_samples.tsv \
  config/rnaseq_design.tsv

# 2. Execute exploratory QC and render the Quarto report
Rscript scripts/run_rnaseq_exploratory_qc.R \
  config/rnaseq_design.tsv \
  results/rnaseq/qc/quantification.tsv \
  results/rnaseq/matrices/final \
  results/rnaseq/exploratory_qc \
  reports/rnaseq_exploratory_qc.qmd \
  reports/rnaseq_exploratory_qc.html
```

---

## 📊 Pipeline Outputs

After running the pipeline, the following artifacts are generated:
- **`results/rnaseq/matrices/final/`**: Contains raw counts (`counts.tsv`), normalized TPM (`tpm.tsv`), and length matrices ready for downstream Differential Expression (DE) analysis.
- **`results/rnaseq/exploratory_qc/tables/`**: VST-normalized counts and summary statistics for the dataset.
- **`results/rnaseq/exploratory_qc/figures/`**: High-quality (SVG/PNG) figures covering read depth, PCA clustering, Euclidean distance heatmaps, and sample-to-sample correlations.
- **`reports/rnaseq_exploratory_qc.html`**: A fully interactive HTML report containing exploratory metadata, diagnostic warnings, and integrated plots for intuitive quality control.

---

## 📄 License & Attribution

Internal research pipeline for HSP/HSF bioinformatic analysis.
