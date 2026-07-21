script_path <- normalizePath("scripts/import_salmon_tximport.R")
source(script_path)

root <- tempfile("tximport-fixture-")
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
quant_root <- file.path(root, "quant")
tx2gene_root <- file.path(root, "tx2gene")
dir.create(tx2gene_root)

species_names <- c("setaria_viridis", "glycine_max")
sample_rows <- data.frame(
  sample_id = c("setaria_c", "setaria_n", "soy_c", "soy_n"),
  species = rep(species_names, each = 2L),
  canary = rep(c("true", "false"), 2L),
  condition = c("control", "stress", "control", "stress"),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(root, "manifest.tsv")
write.table(sample_rows, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

write_quant <- function(species, sample, counts, tpm) {
  directory <- file.path(quant_root, species, sample)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  table <- data.frame(
    Name = c("transcript:t1", "transcript:t2", "transcript:t3"),
    Length = c(100, 200, 300),
    EffectiveLength = c(80, 180, 280),
    TPM = tpm,
    NumReads = counts,
    check.names = FALSE
  )
  write.table(
    table, file.path(directory, "quant.sf"), sep = "\t", quote = FALSE,
    row.names = FALSE
  )
}

for (species in species_names) {
  write.table(
    data.frame(
      transcript_id = c("transcript:t1", "transcript:t2", "transcript:t3"),
      gene_id = c("g1", "g1", "g2")
    ),
    file.path(tx2gene_root, paste0(species, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
write_quant("setaria_viridis", "setaria_c", c(10, 20, 30), c(100000, 200000, 700000))
write_quant("setaria_viridis", "setaria_n", c(5, 15, 40), c(50000, 150000, 800000))
write_quant("glycine_max", "soy_c", c(12, 18, 20), c(120000, 180000, 700000))
write_quant("glycine_max", "soy_n", c(7, 13, 30), c(70000, 130000, 800000))

run_cli <- function(output_root, scope, manifest = manifest_path) {
  output <- system2(
    "Rscript",
    c(script_path, manifest, quant_root, tx2gene_root, output_root, scope),
    stdout = TRUE, stderr = TRUE
  )
  list(status = if (is.null(attr(output, "status"))) 0L else attr(output, "status"), output = output)
}

read_matrix <- function(output_root, species, filename) {
  read.delim(file.path(output_root, species, filename), check.names = FALSE)
}

snapshot <- function(directory) {
  paths <- sort(list.files(directory, recursive = TRUE, full.names = TRUE))
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  setNames(unname(tools::md5sum(paths)), substring(paths, nchar(directory) + 2L))
}

# Direct import: exact values, multiple genes, deterministic feature/sample order.
direct_output <- file.path(root, "direct")
setaria_rows <- sample_rows[sample_rows$species == "setaria_viridis", , drop = FALSE]
result <- import_species(
  setaria_rows, quant_root,
  file.path(tx2gene_root, "setaria_viridis.tsv"), direct_output
)
expected_outputs <- c(
  "transcript_counts.tsv", "transcript_tpm.tsv", "gene_counts.tsv",
  "gene_tpm.tsv", "effective_lengths.tsv", "sample_metadata.tsv"
)
stopifnot(all(file.exists(file.path(direct_output, expected_outputs))))
transcript_counts <- read.delim(file.path(direct_output, "transcript_counts.tsv"), check.names = FALSE)
stopifnot(identical(transcript_counts$feature_id, paste0("transcript:t", 1:3)))
stopifnot(identical(names(transcript_counts), c("feature_id", "setaria_c", "setaria_n")))
stopifnot(identical(as.numeric(transcript_counts$setaria_c), c(10, 20, 30)))
transcript_tpm <- read.delim(file.path(direct_output, "transcript_tpm.tsv"), check.names = FALSE)
stopifnot(identical(as.numeric(transcript_tpm$setaria_n), c(50000, 150000, 800000)))
gene_counts <- read.delim(file.path(direct_output, "gene_counts.tsv"), check.names = FALSE)
stopifnot(identical(gene_counts$feature_id, c("g1", "g2")))
stopifnot(identical(as.numeric(gene_counts$setaria_c), c(30, 30)))
stopifnot(identical(as.numeric(gene_counts$setaria_n), c(20, 40)))
gene_tpm <- read.delim(file.path(direct_output, "gene_tpm.tsv"), check.names = FALSE)
stopifnot(identical(as.numeric(gene_tpm$setaria_c), c(300000, 700000)))
effective_lengths <- read.delim(file.path(direct_output, "effective_lengths.tsv"), check.names = FALSE)
stopifnot(identical(effective_lengths$feature_id, c("g1", "g2")))
stopifnot(all(is.finite(as.matrix(effective_lengths[-1]))))
stopifnot(abs(effective_lengths$setaria_c[2] - 280) < 1e-8)
metadata <- read.delim(file.path(direct_output, "sample_metadata.tsv"), check.names = FALSE)
stopifnot(identical(metadata$sample_id, c("setaria_c", "setaria_n")))
stopifnot(identical(colnames(result$transcript$counts), metadata$sample_id))

# CLI canary/all scope and species isolation.
canary_output <- file.path(root, "canary")
canary_run <- run_cli(canary_output, "canary")
stopifnot(canary_run$status == 0L)
stopifnot(identical(sort(list.dirs(canary_output, recursive = FALSE, full.names = FALSE)), sort(species_names)))
for (species in species_names) {
  species_metadata <- read.delim(
    file.path(canary_output, species, "sample_metadata.tsv"), check.names = FALSE
  )
  stopifnot(nrow(species_metadata) == 1L, identical(species_metadata$canary, "true"))
  counts <- read_matrix(canary_output, species, "transcript_counts.tsv")
  stopifnot(ncol(counts) == 2L)
}

all_output <- file.path(root, "all")
all_run <- run_cli(all_output, "all")
stopifnot(all_run$status == 0L)
for (species in species_names) {
  counts <- read_matrix(all_output, species, "transcript_counts.tsv")
  expected_samples <- sample_rows$sample_id[sample_rows$species == species]
  stopifnot(identical(names(counts), c("feature_id", expected_samples)))
}

# A successful rerun is content-idempotent and leaves no staging/backup siblings.
before_rerun <- snapshot(all_output)
rerun <- run_cli(all_output, "all")
stopifnot(rerun$status == 0L, identical(snapshot(all_output), before_rerun))
transaction_debris <- list.files(
  dirname(all_output),
  pattern = paste0("^\\.", basename(all_output), "\\.(stage|backup)-"),
  full.names = TRUE
)
stopifnot(length(transaction_debris) == 0L)

# RED regression: later-species failure must not publish changed earlier species.
write_quant("setaria_viridis", "setaria_c", c(999, 20, 30), c(100000, 200000, 700000))
unlink(file.path(quant_root, "glycine_max", "soy_n", "quant.sf"))
before_failure <- snapshot(all_output)
failed_run <- suppressWarnings(run_cli(all_output, "all"))
stopifnot(failed_run$status != 0L)
stopifnot(identical(snapshot(all_output), before_failure))

# Invalid scopes and empty selections fail without publishing output.
invalid_output <- file.path(root, "invalid")
invalid_run <- suppressWarnings(run_cli(invalid_output, "invalid"))
stopifnot(invalid_run$status != 0L, !file.exists(invalid_output))
empty_manifest <- file.path(root, "empty-selection.tsv")
empty_rows <- sample_rows
empty_rows$canary <- "false"
write.table(empty_rows, empty_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
empty_output <- file.path(root, "empty")
empty_run <- suppressWarnings(run_cli(empty_output, "canary", empty_manifest))
stopifnot(empty_run$status != 0L, !file.exists(empty_output))
