suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

read_gene_counts <- function(path) {
  table <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                      colClasses = c("character", rep("numeric", length(read.delim(path, nrows = 1)) - 1L)))
  if (names(table)[1] != "feature_id") stop("First count column must be feature_id")
  if (anyDuplicated(table$feature_id) || any(!nzchar(table$feature_id))) stop("Invalid feature_id")
  matrix <- as.matrix(table[-1]); rownames(matrix) <- table$feature_id
  storage.mode(matrix) <- "double"
  matrix
}

validate_group_inputs <- function(counts, metadata, salmon_qc) {
  if (!identical(colnames(counts), metadata$sample_id)) stop("Count columns are not in metadata order")
  if (!identical(metadata$sample_id, salmon_qc$sample_id)) stop("Salmon QC rows are not in metadata order")
  if (any(!is.finite(counts)) || any(counts < 0)) stop("Counts must be finite and nonnegative")
  if (anyDuplicated(metadata$sample_id)) stop("Duplicate sample metadata")
  invisible(TRUE)
}

filter_gene_counts <- function(counts, metadata) {
  cell_sizes <- table(metadata$design_cell)
  minimum_replicates <- min(cell_sizes)
  keep <- rowSums(counts >= 10) >= minimum_replicates
  list(counts = counts[keep, , drop = FALSE], keep = keep,
       minimum_replicates = minimum_replicates)
}

mad_flags <- function(sample_ids, metric, values, direction = "both") {
  if (length(values) < 6L) return(data.frame())
  center <- median(values)
  spread <- mad(values, constant = 1.4826)
  if (!is.finite(spread) || spread == 0) return(data.frame())
  lower <- center - 3 * spread
  upper <- center + 3 * spread
  selected <- switch(direction,
    low = values < lower,
    high = values > upper,
    both = values < lower | values > upper
  )
  if (!any(selected)) return(data.frame())
  data.frame(
    sample_id = sample_ids[selected], metric = metric,
    observed_value = values[selected],
    rule = paste0(direction, " 3-MAD bounds [", signif(lower, 6), ", ", signif(upper, 6), "]"),
    severity = "review", explanation = paste("Advisory", metric, "deviation"),
    stringsAsFactors = FALSE
  )
}

build_advisory_flags <- function(sample_metrics, pca_matrix, correlation) {
  sample_ids <- sample_metrics$sample_id
  median_correlation <- vapply(seq_along(sample_ids), function(i) {
    median(correlation[i, -i])
  }, numeric(1))
  component_count <- min(5L, ncol(pca_matrix))
  centered <- sweep(pca_matrix[, seq_len(component_count), drop = FALSE], 2,
                    colMeans(pca_matrix[, seq_len(component_count), drop = FALSE]))
  pca_distance <- sqrt(rowSums(centered^2))
  bind_rows(
    mad_flags(sample_ids, "log10_library_size", log10(sample_metrics$raw_library_size + 1), "both"),
    mad_flags(sample_ids, "detected_ge_10", sample_metrics$detected_ge_10, "low"),
    mad_flags(sample_ids, "median_spearman_correlation", median_correlation, "low"),
    mad_flags(sample_ids, "pca_centroid_distance", pca_distance, "high")
  )
}

calculate_group_qc <- function(counts, metadata, salmon_qc) {
  validate_group_inputs(counts, metadata, salmon_qc)
  filtered <- filter_gene_counts(counts, metadata)
  dds <- DESeqDataSetFromMatrix(
    countData = round(filtered$counts), colData = metadata, design = ~ 1
  )
  dds <- estimateSizeFactors(dds)
  transformed <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  pca <- prcomp(t(transformed), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_scores <- data.frame(
    sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  ) %>% left_join(metadata, by = "sample_id")
  correlation <- cor(transformed, method = "spearman")
  distance <- as.matrix(dist(t(transformed)))
  sample_metrics <- metadata %>% transmute(
    sample_id, species, bioproject, analysis_group, design_cell,
    raw_library_size = colSums(counts),
    detected_ge_1 = colSums(counts >= 1),
    detected_ge_10 = colSums(counts >= 10),
    size_factor = sizeFactors(dds)
  ) %>% left_join(salmon_qc, by = c("sample_id", "species"))
  distributional_flags <- build_advisory_flags(sample_metrics, pca$x, correlation)
  mapping_flags <- data.frame(
    sample_id = salmon_qc$sample_id[salmon_qc$mapping_flag != "pass"],
    metric = "percent_mapped",
    observed_value = salmon_qc$percent_mapped[salmon_qc$mapping_flag != "pass"],
    rule = "mapping_flag != pass", severity = "warning",
    explanation = "Preserved Salmon mapping warning", stringsAsFactors = FALSE
  )
  flags <- bind_rows(distributional_flags, mapping_flags)
  list(
    filtered_counts = filtered$counts, keep = filtered$keep,
    minimum_replicates = filtered$minimum_replicates,
    transformed = transformed, sample_metrics = sample_metrics,
    pca_scores = pca_scores, pca_variance = variance,
    correlation = correlation, distance = distance,
    correlation_long = as.data.frame(as.table(correlation), stringsAsFactors = FALSE),
    distance_long = as.data.frame(as.table(distance), stringsAsFactors = FALSE),
    advisory_flags = flags
  )
}
