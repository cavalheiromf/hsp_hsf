#!/usr/bin/env Rscript

write_matrix <- function(matrix, path) {
  output <- data.frame(feature_id = rownames(matrix), matrix, check.names = FALSE)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

require_columns <- function(table, required, label) {
  missing <- setdiff(required, names(table))
  if (length(missing) > 0L) {
    stop(label, " lacks required columns: ", paste(missing, collapse = ", "))
  }
}

import_species <- function(sample_rows, quant_root, tx2gene_path, output_dir) {
  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop("Bioconductor package 'tximport' is required")
  }
  require_columns(sample_rows, c("sample_id", "species"), "sample metadata")
  if (nrow(sample_rows) == 0L) stop("sample metadata must contain at least one row")
  if (anyNA(sample_rows$sample_id) || any(!nzchar(sample_rows$sample_id))) {
    stop("sample_id values must be non-empty")
  }
  if (anyDuplicated(sample_rows$sample_id)) stop("sample_id values must be unique")
  species <- unique(sample_rows$species)
  if (length(species) != 1L || is.na(species) || !nzchar(species)) {
    stop("sample metadata must contain exactly one non-empty species")
  }

  files <- file.path(quant_root, species, sample_rows$sample_id, "quant.sf")
  names(files) <- sample_rows$sample_id
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("Missing quant.sf: ", paste(missing, collapse = ", "))
  }

  if (!file.exists(tx2gene_path)) stop("Missing tx2gene: ", tx2gene_path)
  tx2gene <- read.delim(
    tx2gene_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character"
  )
  if (!identical(names(tx2gene), c("transcript_id", "gene_id"))) {
    stop("tx2gene must contain exactly transcript_id and gene_id columns")
  }
  if (nrow(tx2gene) == 0L || anyNA(tx2gene) || any(!nzchar(tx2gene$transcript_id)) ||
      any(!nzchar(tx2gene$gene_id))) {
    stop("tx2gene identifiers must be non-empty")
  }
  if (anyDuplicated(tx2gene$transcript_id)) {
    stop("tx2gene transcript_id values must be unique")
  }

  expected_transcripts <- tx2gene$transcript_id
  for (file in files) {
    quantified <- read.delim(
      file, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = c(Name = "character")
    )
    require_columns(quantified, "Name", paste0("quant.sf ", file))
    if (anyNA(quantified$Name) || any(!nzchar(quantified$Name)) ||
        anyDuplicated(quantified$Name)) {
      stop("quant.sf transcript names must be non-empty and unique: ", file)
    }
    missing_mapping <- setdiff(quantified$Name, expected_transcripts)
    extra_mapping <- setdiff(expected_transcripts, quantified$Name)
    if (length(missing_mapping) > 0L || length(extra_mapping) > 0L) {
      stop(
        "Transcript/tx2gene mismatch for ", file,
        ": missing=", length(missing_mapping), ", extra=", length(extra_mapping)
      )
    }
  }

  transcript <- tximport::tximport(
    files, type = "salmon", txOut = TRUE, countsFromAbundance = "no"
  )
  gene <- tximport::tximport(
    files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no"
  )
  expected_samples <- sample_rows$sample_id
  matrices <- list(
    transcript$counts, transcript$abundance, gene$counts, gene$abundance, gene$length
  )
  if (!all(vapply(
    matrices, function(x) identical(colnames(x), expected_samples), logical(1)
  ))) {
    stop("tximport output columns are not aligned to sample metadata")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_matrix(transcript$counts, file.path(output_dir, "transcript_counts.tsv"))
  write_matrix(transcript$abundance, file.path(output_dir, "transcript_tpm.tsv"))
  write_matrix(gene$counts, file.path(output_dir, "gene_counts.tsv"))
  write_matrix(gene$abundance, file.path(output_dir, "gene_tpm.tsv"))
  write_matrix(gene$length, file.path(output_dir, "effective_lengths.tsv"))
  write.table(
    sample_rows, file.path(output_dir, "sample_metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(list(transcript = transcript, gene = gene))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5L) {
    stop(
      "Usage: import_salmon_tximport.R ",
      "MANIFEST QUANT_ROOT TX2GENE_ROOT OUTPUT_ROOT SCOPE"
    )
  }
  manifest <- read.delim(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
  require_columns(manifest, c("sample_id", "species", "canary"), "manifest")
  scope <- args[[5L]]
  if (!scope %in% c("canary", "all")) stop("Scope must be canary or all")
  if (scope == "canary") {
    manifest <- manifest[manifest$canary == "true", , drop = FALSE]
  }
  if (nrow(manifest) == 0L) stop("No samples selected for scope: ", scope)

  for (species in unique(manifest$species)) {
    rows <- manifest[manifest$species == species, , drop = FALSE]
    import_species(
      rows,
      args[[2L]],
      file.path(args[[3L]], paste0(species, ".tsv")),
      file.path(args[[4L]], species)
    )
  }
}

if (sys.nframe() == 0L) main()
