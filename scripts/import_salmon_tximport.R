#!/usr/bin/env Rscript

write_matrix <- function(matrix, path) {
  output <- data.frame(feature_id = rownames(matrix), matrix, check.names = FALSE)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_output_matrix <- function(path, expected_features, expected_samples, positive = FALSE) {
  table <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  expected_columns <- c("feature_id", expected_samples)
  if (!identical(names(table), expected_columns)) {
    stop("Matrix columns are not aligned: ", path)
  }
  if (!identical(table$feature_id, expected_features)) {
    stop("Matrix features are not aligned: ", path)
  }
  values <- as.matrix(table[, expected_samples, drop = FALSE])
  storage.mode(values) <- "numeric"
  if (any(!is.finite(values)) || any(values < 0) || (positive && any(values <= 0))) {
    stop("Matrix contains invalid numeric values: ", path)
  }
  invisible(table)
}

validate_species_output <- function(output_dir, sample_rows, tx2gene) {
  expected_files <- c(
    "transcript_counts.tsv", "transcript_tpm.tsv", "gene_counts.tsv",
    "gene_tpm.tsv", "effective_lengths.tsv", "sample_metadata.tsv"
  )
  paths <- file.path(output_dir, expected_files)
  missing <- paths[!file.exists(paths) | file.info(paths)$size == 0]
  if (length(missing) > 0L) {
    stop("Missing or empty staged output: ", paste(missing, collapse = ", "))
  }

  samples <- sample_rows$sample_id
  transcripts <- tx2gene$transcript_id
  genes <- unique(tx2gene$gene_id)
  read_output_matrix(paths[[1L]], transcripts, samples)
  read_output_matrix(paths[[2L]], transcripts, samples)
  read_output_matrix(paths[[3L]], genes, samples)
  read_output_matrix(paths[[4L]], genes, samples)
  read_output_matrix(paths[[5L]], genes, samples, positive = TRUE)

  metadata <- read.delim(paths[[6L]], stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(metadata), names(sample_rows)) || nrow(metadata) != nrow(sample_rows)) {
    stop("Staged sample metadata shape is not aligned: ", paths[[6L]])
  }
  aligned <- vapply(
    names(sample_rows),
    function(column) identical(as.character(metadata[[column]]), as.character(sample_rows[[column]])),
    logical(1)
  )
  if (!all(aligned)) stop("Staged sample metadata values are not aligned: ", paths[[6L]])
  invisible(TRUE)
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
  transcript_order <- match(tx2gene$transcript_id, rownames(transcript$counts))
  gene_ids <- unique(tx2gene$gene_id)
  gene_order <- match(gene_ids, rownames(gene$counts))
  if (anyNA(transcript_order) || anyNA(gene_order)) {
    stop("tximport output does not contain every expected feature")
  }
  for (name in c("abundance", "counts", "length")) {
    transcript[[name]] <- transcript[[name]][transcript_order, , drop = FALSE]
    gene[[name]] <- gene[[name]][gene_order, , drop = FALSE]
  }
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
  validate_species_output(output_dir, sample_rows, tx2gene)
  invisible(list(transcript = transcript, gene = gene))
}

publish_collection <- function(manifest, quant_root, tx2gene_root, output_root) {
  if (nrow(manifest) == 0L) stop("Cannot publish an empty sample collection")
  require_columns(manifest, c("sample_id", "species"), "selected manifest")
  if (anyNA(manifest$species) || any(!nzchar(manifest$species))) {
    stop("Selected manifest contains an empty species")
  }

  parent <- dirname(output_root)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent)) stop("Cannot create output parent: ", parent)
  prefix <- paste0(".", basename(output_root))
  stage <- tempfile(paste0(prefix, ".stage-"), tmpdir = parent)
  backup <- NULL
  dir.create(stage)
  if (!dir.exists(stage)) stop("Cannot create staging directory: ", stage)
  on.exit({
    if (!is.null(stage) && dir.exists(stage)) unlink(stage, recursive = TRUE)
  }, add = TRUE)

  species_order <- unique(manifest$species)
  mappings <- vector("list", length(species_order))
  names(mappings) <- species_order
  for (species in species_order) {
    rows <- manifest[manifest$species == species, , drop = FALSE]
    mapping_path <- file.path(tx2gene_root, paste0(species, ".tsv"))
    import_species(rows, quant_root, mapping_path, file.path(stage, species))
    mappings[[species]] <- read.delim(
      mapping_path, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = "character"
    )
  }
  for (species in species_order) {
    rows <- manifest[manifest$species == species, , drop = FALSE]
    validate_species_output(file.path(stage, species), rows, mappings[[species]])
  }

  if (file.exists(output_root) && !dir.exists(output_root)) {
    stop("Existing output root is not a directory: ", output_root)
  }
  if (dir.exists(output_root)) {
    backup <- tempfile(paste0(prefix, ".backup-"), tmpdir = parent)
    if (!file.rename(output_root, backup)) {
      stop("Cannot move existing collection to backup: ", output_root)
    }
  }
  if (!file.rename(stage, output_root)) {
    rollback_ok <- is.null(backup) || file.rename(backup, output_root)
    if (!rollback_ok) {
      stop(
        "Collection promotion and rollback failed; previous output remains at ", backup
      )
    }
    stop("Collection promotion failed; previous output was restored")
  }
  stage <- NULL
  if (!is.null(backup) && dir.exists(backup)) {
    cleanup_status <- unlink(backup, recursive = TRUE)
    if (cleanup_status != 0L || dir.exists(backup)) {
      stop("Collection published, but backup cleanup failed: ", backup)
    }
  }
  invisible(output_root)
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

  publish_collection(manifest, args[[2L]], args[[3L]], args[[4L]])
}

if (sys.nframe() == 0L) main()
