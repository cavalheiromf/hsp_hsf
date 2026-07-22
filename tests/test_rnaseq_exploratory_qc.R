source("scripts/rnaseq_design_lib.R")
source("scripts/rnaseq_exploratory_qc_lib.R")

counts <- read_gene_counts("tests/fixtures/rnaseq_qc/matrices/test_species/gene_counts.tsv")
metadata <- read.delim("tests/fixtures/rnaseq_qc/manifest.tsv", stringsAsFactors = FALSE,
                       check.names = FALSE, colClasses = "character")
salmon_qc <- read.delim("tests/fixtures/rnaseq_qc/quantification.tsv",
                        stringsAsFactors = FALSE, check.names = FALSE)

result <- calculate_group_qc(counts, metadata, salmon_qc)
stopifnot(identical(colnames(result$filtered_counts), metadata$sample_id))
stopifnot(identical(rownames(result$filtered_counts), paste0("gene_", letters[1:5])))
stopifnot(nrow(result$sample_metrics) == 6L)
stopifnot(nrow(result$pca_scores) == 6L)
stopifnot(all(c("PC1", "PC2") %in% names(result$pca_scores)))
stopifnot(nrow(result$correlation_long) == 36L)
stopifnot(nrow(result$distance_long) == 36L)
stopifnot(any(result$advisory_flags$sample_id == "s6" &
              result$advisory_flags$metric == "percent_mapped"))
stopifnot(all(metadata$sample_id %in% result$sample_metrics$sample_id))

misordered <- counts[, c(2, 1, 3, 4, 5, 6), drop = FALSE]
stopifnot(inherits(try(calculate_group_qc(misordered, metadata, salmon_qc), silent = TRUE), "try-error"))
negative <- counts; negative[1, 1] <- -1
stopifnot(inherits(try(calculate_group_qc(negative, metadata, salmon_qc), silent = TRUE), "try-error"))
nonfinite <- counts; nonfinite[1, 1] <- Inf
stopifnot(inherits(try(calculate_group_qc(nonfinite, metadata, salmon_qc), silent = TRUE), "try-error"))
