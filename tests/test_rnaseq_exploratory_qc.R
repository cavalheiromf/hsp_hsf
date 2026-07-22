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

output <- tempfile("rnaseq-qc-output-")
dir.create(output)
write_group_tables(result, "fixture_group", file.path(output, "tables"))
write_group_figures(result, "fixture_group", file.path(output, "figures"))

expected_figures <- as.vector(outer(
  paste0("fixture_group_", c(
    "library_size", "detected_genes", "vst_distribution", "pca",
    "correlation_heatmap", "distance_heatmap"
  )), c("png", "svg"), paste, sep = "."
))
stopifnot(identical(sort(list.files(file.path(output, "figures"))), sort(expected_figures)))
stopifnot(all(file.info(file.path(output, "figures", expected_figures))$size > 0))

metrics <- read.delim(file.path(output, "tables", "fixture_group_sample_metrics.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(metrics$sample_id, metadata$sample_id))
scores <- read.delim(file.path(output, "tables", "fixture_group_pca_scores.tsv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(scores$sample_id, metadata$sample_id))

# Rollback test: simulated failure during publication leaves existing root unchanged
old_root <- tempfile("published-qc-")
dir.create(old_root)
writeLines("old sentinel", file.path(old_root, "sentinel.txt"))
report_link <- tempfile("published-report-", fileext = ".html")
writeLines("old report", report_link)

fixture_inputs <- list(
  design = metadata,
  salmon_qc = salmon_qc,
  counts_by_species = list(test_species = counts),
  contrasts = data.frame(
    contrast_id = "c1", species = "test_species", bioproject = "TEST1",
    analysis_group = "fixture_group", model_formula = "~ 1",
    effect_type = "main_effect", numerator = "a", denominator = "b",
    reference_levels = "ref", estimable = "true", notes = "test",
    stringsAsFactors = FALSE
  )
)

failure_hook <- function(stage) stop("injected late failure")
message <- tryCatch({
  publish_qc_collection(fixture_inputs, old_root, report_link,
                        report_qmd = "reports/rnaseq_exploratory_qc.qmd",
                        failure_hook = failure_hook)
  NA_character_
}, error = conditionMessage)
stopifnot(grepl("injected late failure", message, fixed = TRUE))
stopifnot(identical(readLines(file.path(old_root, "sentinel.txt")), "old sentinel"))
stopifnot(identical(readLines(report_link), "old report"))
stopifnot(!length(Sys.glob(paste0(old_root, ".stage-*"))))
stopifnot(!length(Sys.glob(paste0(old_root, ".backup-*"))))

