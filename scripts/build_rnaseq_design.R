#!/usr/bin/env Rscript
source("scripts/rnaseq_design_lib.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: build_rnaseq_design.R MANIFEST OUTPUT")
manifest <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = "character", na.strings = character())
design <- normalize_design(manifest)
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.table(design, args[2], sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
