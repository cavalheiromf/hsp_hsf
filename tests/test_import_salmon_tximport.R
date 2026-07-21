script_path <- "scripts/import_salmon_tximport.R"
stopifnot(file.exists(script_path))
source(script_path)

root <- tempfile("tximport-fixture-")
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
quant_root <- file.path(root, "quant")
output_dir <- file.path(root, "matrices")
tx2gene_path <- file.path(root, "tx2gene.tsv")

write_quant <- function(sample, counts, tpm) {
  directory <- file.path(quant_root, "setaria_viridis", sample)
  dir.create(directory, recursive = TRUE)
  table <- data.frame(
    Name = c("transcript:t1", "transcript:t2"),
    Length = c(100, 200),
    EffectiveLength = c(80, 180),
    TPM = tpm,
    NumReads = counts,
    check.names = FALSE
  )
  write.table(
    table, file.path(directory, "quant.sf"), sep = "\t", quote = FALSE,
    row.names = FALSE
  )
}

write_quant("sample1", c(10, 20), c(400000, 600000))
write_quant("sample2", c(5, 15), c(300000, 700000))
write.table(
  data.frame(
    transcript_id = c("transcript:t1", "transcript:t2"),
    gene_id = c("g1", "g1")
  ),
  tx2gene_path, sep = "\t", quote = FALSE, row.names = FALSE
)
manifest <- data.frame(
  sample_id = c("sample1", "sample2"),
  species = c("setaria_viridis", "setaria_viridis"),
  canary = c("true", "false"),
  stringsAsFactors = FALSE
)

result <- import_species(manifest, quant_root, tx2gene_path, output_dir)

expected_outputs <- c(
  "transcript_counts.tsv", "transcript_tpm.tsv", "gene_counts.tsv",
  "gene_tpm.tsv", "effective_lengths.tsv", "sample_metadata.tsv"
)
stopifnot(all(file.exists(file.path(output_dir, expected_outputs))))

transcript_counts <- read.delim(
  file.path(output_dir, "transcript_counts.tsv"), check.names = FALSE
)
stopifnot(identical(transcript_counts$feature_id, c("transcript:t1", "transcript:t2")))
stopifnot(identical(names(transcript_counts), c("feature_id", "sample1", "sample2")))
stopifnot(all(abs(transcript_counts$sample1 - c(10, 20)) < 1e-8))

gene_counts <- read.delim(file.path(output_dir, "gene_counts.tsv"), check.names = FALSE)
stopifnot(identical(gene_counts$feature_id, "g1"))
stopifnot(abs(gene_counts$sample1 - 30) < 1e-8)
stopifnot(abs(gene_counts$sample2 - 20) < 1e-8)

gene_tpm <- read.delim(file.path(output_dir, "gene_tpm.tsv"), check.names = FALSE)
stopifnot(abs(gene_tpm$sample1 - 1e6) < 1e-8)
effective_lengths <- read.delim(
  file.path(output_dir, "effective_lengths.tsv"), check.names = FALSE
)
stopifnot(identical(names(effective_lengths), c("feature_id", "sample1", "sample2")))
stopifnot(all(is.finite(effective_lengths$sample1)))

metadata <- read.delim(file.path(output_dir, "sample_metadata.tsv"), check.names = FALSE)
stopifnot(identical(metadata$sample_id, c("sample1", "sample2")))
stopifnot(identical(colnames(result$transcript$counts), metadata$sample_id))

bad_mapping <- file.path(root, "bad-tx2gene.tsv")
write.table(
  data.frame(transcript_id = "transcript:t1", gene_id = "g1"),
  bad_mapping, sep = "\t", quote = FALSE, row.names = FALSE
)
mismatch <- tryCatch(
  import_species(manifest, quant_root, bad_mapping, file.path(root, "bad-output")),
  error = identity
)
stopifnot(inherits(mismatch, "error"))
stopifnot(grepl("Transcript/tx2gene mismatch", conditionMessage(mismatch), fixed = TRUE))

duplicate_samples <- rbind(manifest[1, ], manifest[1, ])
duplicate_error <- tryCatch(
  import_species(
    duplicate_samples, quant_root, tx2gene_path, file.path(root, "duplicate-output")
  ),
  error = identity
)
stopifnot(inherits(duplicate_error, "error"))
stopifnot(grepl("unique", conditionMessage(duplicate_error), fixed = TRUE))
