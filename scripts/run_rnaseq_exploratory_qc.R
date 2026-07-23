#!/usr/bin/env Rscript
required_packages <- c("BiocManager", "DESeq2", "dplyr", "ggplot2", "pheatmap", "tidyverse")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
if (!nzchar(Sys.which("quarto"))) stop("Quarto executable is required")

source("scripts/rnaseq_design_lib.R")
source("scripts/rnaseq_exploratory_qc_lib.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: run_rnaseq_exploratory_qc.R DESIGN SALMON_QC MATRIX_ROOT OUTPUT_ROOT REPORT_QMD REPORT_HTML")
}
# read_character_tsv is provided by utils_io.R (via rnaseq_design_lib.R)
design <- read_character_tsv(args[1])
validate_design(design)
salmon_qc <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE,
                        na.strings = character())
if (!identical(design$sample_id, salmon_qc$sample_id)) stop("Design and Salmon QC order differ")
species_order <- unique(design$species)
counts_by_species <- setNames(lapply(species_order, function(species) {
  path <- file.path(args[3], species, "gene_counts.tsv")
  counts <- read_gene_counts(path)
  expected <- design$sample_id[design$species == species]
  if (!identical(colnames(counts), expected)) stop("Species count order differs for ", species)
  counts
}), species_order)
inputs <- list(
  design = design, salmon_qc = salmon_qc,
  counts_by_species = counts_by_species,
  contrasts = build_contrast_inventory(design)
)
dir.create(dirname(args[6]), recursive = TRUE, showWarnings = FALSE)
publish_qc_collection(inputs, args[4], args[6], args[5])
