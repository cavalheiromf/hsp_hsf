source("scripts/utils_io.R")
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
  warn_idx <- which(salmon_qc$mapping_flag != "pass")
  mapping_flags <- if (length(warn_idx) > 0L) {
    data.frame(
      sample_id = salmon_qc$sample_id[warn_idx],
      metric = "percent_mapped",
      observed_value = salmon_qc$percent_mapped[warn_idx],
      rule = "mapping_flag != pass", severity = "warning",
      explanation = "Preserved Salmon mapping warning", stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      sample_id = character(), metric = character(), observed_value = numeric(),
      rule = character(), severity = character(), explanation = character(),
      stringsAsFactors = FALSE
    )
  }
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

# write_tsv_safe is provided by utils_io.R

write_group_tables <- function(result, group, table_dir) {
  correlation <- result$correlation_long
  names(correlation) <- c("sample_id_1", "sample_id_2", "spearman_correlation")
  distance <- result$distance_long
  names(distance) <- c("sample_id_1", "sample_id_2", "euclidean_distance")
  filter_summary <- data.frame(
    analysis_group = group,
    genes_input = length(result$keep),
    genes_retained = sum(result$keep),
    genes_removed = sum(!result$keep),
    minimum_replicates = result$minimum_replicates,
    stringsAsFactors = FALSE
  )
  write_tsv_safe(result$sample_metrics,
                 file.path(table_dir, paste0(group, "_sample_metrics.tsv")))
  write_tsv_safe(result$pca_scores,
                 file.path(table_dir, paste0(group, "_pca_scores.tsv")))
  write_tsv_safe(correlation,
                 file.path(table_dir, paste0(group, "_correlation.tsv")))
  write_tsv_safe(distance,
                 file.path(table_dir, paste0(group, "_distance.tsv")))
  write_tsv_safe(filter_summary,
                 file.path(table_dir, paste0(group, "_gene_filter_summary.tsv")))
  write_tsv_safe(result$advisory_flags,
                 file.path(table_dir, paste0(group, "_advisory_flags.tsv")))
}

write_combined_tables <- function(group_results, group_order, contrasts, output_root) {
  ordered <- group_results[group_order]
  combined <- function(field) bind_rows(lapply(ordered, `[[`, field))
  filter_summary <- bind_rows(lapply(group_order, function(group) {
    result <- ordered[[group]]
    data.frame(
      analysis_group = group, genes_input = length(result$keep),
      genes_retained = sum(result$keep), genes_removed = sum(!result$keep),
      minimum_replicates = result$minimum_replicates,
      stringsAsFactors = FALSE
    )
  }))
  table_dir <- file.path(output_root, "tables")
  write_tsv_safe(combined("sample_metrics"), file.path(table_dir, "sample_metrics.tsv"))
  write_tsv_safe(filter_summary, file.path(table_dir, "gene_filter_summary.tsv"))
  write_tsv_safe(combined("pca_scores"), file.path(table_dir, "pca_scores.tsv"))
  flags <- combined("advisory_flags")
  if (!ncol(flags)) flags <- data.frame(
    sample_id = character(), metric = character(), observed_value = numeric(),
    rule = character(), severity = character(), explanation = character()
  )
  write_tsv_safe(flags, file.path(table_dir, "advisory_flags.tsv"))
  write_tsv_safe(contrasts, file.path(output_root, "contrasts.tsv"))
}

save_plot_pair <- function(plot, stem, width = 8, height = 5) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         units = "in", dpi = 150, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height,
         units = "in", bg = "white")
}

heatmap_pair <- function(matrix, annotation, stem) {
  palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  for (extension in c("png", "svg")) {
    path <- paste0(stem, ".", extension)
    if (extension == "png") png(path, width = 1600, height = 1400, res = 180)
    else svg(path, width = 9, height = 8)
    pheatmap::pheatmap(matrix, annotation_col = annotation,
                       annotation_row = annotation, cluster_rows = FALSE,
                       cluster_cols = FALSE, color = palette, border_color = NA)
    dev.off()
  }
}

write_group_figures <- function(result, group, figure_dir) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  metrics <- result$sample_metrics
  base_theme <- theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  library_plot <- ggplot(metrics, aes(sample_id, raw_library_size, fill = design_cell)) +
    geom_col() + scale_y_log10() + labs(x = NULL, y = "Raw gene-count library size") + base_theme
  detected_plot <- ggplot(metrics, aes(sample_id, detected_ge_10, fill = design_cell)) +
    geom_col() + labs(x = NULL, y = "Genes with count >= 10") + base_theme
  vst_long <- as.data.frame(as.table(result$transformed), stringsAsFactors = FALSE)
  names(vst_long) <- c("feature_id", "sample_id", "vst")
  vst_plot <- ggplot(vst_long, aes(vst, color = sample_id)) + geom_density(linewidth = 0.5) +
    labs(x = "Blind VST expression", y = "Density") + theme_bw(base_size = 11)
  pca_plot <- ggplot(result$pca_scores, aes(PC1, PC2, color = design_cell, label = sample_id)) +
    geom_point(size = 3) + geom_text(vjust = -0.7, check_overlap = TRUE) +
    labs(x = sprintf("PC1 (%.1f%%)", result$pca_variance[1]),
         y = sprintf("PC2 (%.1f%%)", result$pca_variance[2])) + theme_bw(base_size = 11)
  save_plot_pair(library_plot, file.path(figure_dir, paste0(group, "_library_size")))
  save_plot_pair(detected_plot, file.path(figure_dir, paste0(group, "_detected_genes")))
  save_plot_pair(vst_plot, file.path(figure_dir, paste0(group, "_vst_distribution")))
  save_plot_pair(pca_plot, file.path(figure_dir, paste0(group, "_pca")))
  annotation <- data.frame(design_cell = metrics$design_cell, row.names = metrics$sample_id)
  heatmap_pair(result$correlation, annotation,
               file.path(figure_dir, paste0(group, "_correlation_heatmap")))
  heatmap_pair(result$distance, annotation,
               file.path(figure_dir, paste0(group, "_distance_heatmap")))
}

write_joint_species_pca <- function(counts, metadata, species_label, figure_dir) {
  if (!species_label %in% c("wheat", "soybean")) stop("Unsupported joint PCA label")
  if (!identical(colnames(counts), metadata$sample_id)) stop("Joint PCA sample order mismatch")
  minimum_project_size <- min(table(metadata$bioproject))
  keep <- rowSums(counts >= 10) >= minimum_project_size
  dds <- DESeqDataSetFromMatrix(round(counts[keep, , drop = FALSE]), metadata, ~ 1)
  dds <- estimateSizeFactors(dds)
  transformed <- assay(varianceStabilizingTransformation(dds, blind = TRUE))
  pca <- prcomp(t(transformed))
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2]) %>%
    left_join(metadata, by = "sample_id")
  plot <- ggplot(scores, aes(PC1, PC2, color = bioproject, label = sample_id)) +
    geom_point(size = 3) + geom_text(vjust = -0.7, check_overlap = TRUE) +
    labs(x = sprintf("PC1 (%.1f%%)", variance[1]),
         y = sprintf("PC2 (%.1f%%)", variance[2])) + theme_bw(base_size = 11)
  save_plot_pair(plot, file.path(figure_dir, paste0(species_label, "_joint_bioproject_pca")))
}

validate_qc_collection <- function(stage, expected_design) {
  required_tables <- c(
    "tables/sample_metrics.tsv", "tables/gene_filter_summary.tsv",
    "tables/pca_scores.tsv", "tables/advisory_flags.tsv", "contrasts.tsv"
  )
  groups <- unique(expected_design$analysis_group)
  stems <- c("library_size", "detected_genes", "vst_distribution", "pca",
             "correlation_heatmap", "distance_heatmap")
  group_figures <- unlist(lapply(groups, function(group) {
    as.vector(outer(paste0(group, "_", stems), c("png", "svg"), paste, sep = "."))
  }))
  joint_labels <- c(triticum_aestivum = "wheat", glycine_max = "soybean")
  joint_labels <- joint_labels[names(joint_labels) %in% expected_design$species]
  joint_labels <- joint_labels[vapply(names(joint_labels), function(species) {
    length(unique(expected_design$bioproject[expected_design$species == species])) > 1L
  }, logical(1))]
  joint <- if (length(joint_labels)) as.vector(outer(
    paste0(unname(joint_labels), "_joint_bioproject_pca"),
    c("png", "svg"), paste, sep = "."
  )) else character()
  required <- c(required_tables, "session_info.txt",
                file.path("figures", c(group_figures, joint)),
                "report/rnaseq_exploratory_qc.html")
  missing <- required[!file.exists(file.path(stage, required))]
  if (length(missing)) stop("Staged QC collection missing: ", paste(missing, collapse = ", "))
  if (any(file.info(file.path(stage, required))$size <= 0)) stop("Staged QC collection contains empty files")
  metrics <- read.delim(file.path(stage, "tables", "sample_metrics.tsv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  flags <- read.delim(file.path(stage, "tables", "advisory_flags.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(metrics$sample_id, expected_design$sample_id)) stop("QC metrics sample order mismatch")
  if (nrow(metrics) != nrow(expected_design) || anyDuplicated(metrics$sample_id)) {
    stop("QC metrics must contain every expected sample exactly once")
  }
  if (!identical(unique(metrics$analysis_group), groups)) stop("QC analysis-group order mismatch")
  numeric_columns <- vapply(metrics, is.numeric, logical(1))
  if (any(!is.finite(as.matrix(metrics[numeric_columns])))) stop("Nonfinite sample metric")
  accepted <- intersect(c("SRR39669466", "SRR39669467"), expected_design$sample_id)
  mapping_flags <- flags[flags$metric == "percent_mapped", "sample_id"]
  if (!all(accepted %in% mapping_flags)) stop("Accepted Setaria warnings are missing")
  invisible(TRUE)
}

relative_path <- function(target, from) {
  target_parts <- strsplit(normalizePath(target, mustWork = FALSE), "/", fixed = TRUE)[[1]]
  from_parts <- strsplit(normalizePath(from, mustWork = TRUE), "/", fixed = TRUE)[[1]]
  common <- 0L
  limit <- min(length(target_parts), length(from_parts))
  while (common < limit && target_parts[common + 1L] == from_parts[common + 1L]) {
    common <- common + 1L
  }
  up <- rep("..", length(from_parts) - common)
  down <- target_parts[(common + 1L):length(target_parts)]
  paste(c(up, down), collapse = "/")
}

publish_qc_collection <- function(inputs, output_root, report_html, report_qmd,
                                  failure_hook = function(stage) invisible(NULL),
                                  rename_fn = file.rename) {
  parent <- dirname(output_root)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(paste0(".", basename(output_root), ".stage-"), tmpdir = parent)
  dir.create(stage)
  output_backup <- NULL
  report_backup <- NULL
  on.exit(if (!is.null(stage) && dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)

  build_qc_collection(inputs, stage, report_qmd)
  validate_qc_collection(stage, inputs$design)
  failure_hook(stage)

  if (dir.exists(output_root)) {
    output_backup <- tempfile(paste0(".", basename(output_root), ".backup-"), tmpdir = parent)
    if (!rename_fn(output_root, output_backup)) stop("Cannot back up QC collection")
  }
  if (file.exists(report_html)) {
    dir.create(dirname(report_html), recursive = TRUE, showWarnings = FALSE)
    report_backup <- tempfile(".rnaseq-qc-report.backup-", tmpdir = dirname(report_html))
    if (!rename_fn(report_html, report_backup)) {
      if (!is.null(output_backup)) rename_fn(output_backup, output_root)
      stop("Cannot back up QC report link")
    }
  } else if (nzchar(Sys.readlink(report_html))) {
    unlink(report_html)
  }
  if (!rename_fn(stage, output_root)) {
    if (!is.null(output_backup)) rename_fn(output_backup, output_root)
    if (!is.null(report_backup)) rename_fn(report_backup, report_html)
    stop("Cannot promote QC collection")
  }
  stage <- NULL

  relative_target <- relative_path(
    file.path(output_root, "report", "rnaseq_exploratory_qc.html"),
    dirname(report_html)
  )
  temporary_link <- tempfile(".rnaseq-qc-report.link-", tmpdir = dirname(report_html))
  if (!file.symlink(relative_target, temporary_link) || !rename_fn(temporary_link, report_html)) {
    failed_new <- tempfile(".rnaseq-qc.failed-new-", tmpdir = parent)
    rename_fn(output_root, failed_new)
    if (!is.null(output_backup)) rename_fn(output_backup, output_root)
    if (!is.null(report_backup)) rename_fn(report_backup, report_html)
    unlink(failed_new, recursive = TRUE)
    stop("Report-link promotion failed; previous collection restored")
  }
  if (!is.null(output_backup)) unlink(output_backup, recursive = TRUE)
  if (!is.null(report_backup)) unlink(report_backup)
  invisible(output_root)
}

build_qc_collection <- function(inputs, stage, report_qmd) {
  group_order <- unique(inputs$design$analysis_group)
  results <- setNames(vector("list", length(group_order)), group_order)
  for (group in group_order) {
    metadata <- inputs$design[inputs$design$analysis_group == group, , drop = FALSE]
    species <- unique(metadata$species)
    if (length(species) != 1L) stop("Analysis group contains multiple species")
    counts <- inputs$counts_by_species[[species]][, metadata$sample_id, drop = FALSE]
    salmon <- inputs$salmon_qc[match(metadata$sample_id, inputs$salmon_qc$sample_id), , drop = FALSE]
    result <- calculate_group_qc(counts, metadata, salmon)
    results[[group]] <- result
    write_group_tables(result, group, file.path(stage, "tables", "groups"))
    write_group_figures(result, group, file.path(stage, "figures"))
  }
  write_combined_tables(results, group_order, inputs$contrasts, stage)
  for (entry in list(
    list(species = "triticum_aestivum", label = "wheat"),
    list(species = "glycine_max", label = "soybean")
  )) {
    metadata <- inputs$design[inputs$design$species == entry$species, , drop = FALSE]
    if (!nrow(metadata) || length(unique(metadata$bioproject)) < 2L) next
    counts <- inputs$counts_by_species[[entry$species]][, metadata$sample_id, drop = FALSE]
    write_joint_species_pca(counts, metadata, entry$label, file.path(stage, "figures"))
  }
  session <- c(
    capture.output(sessionInfo()),
    paste("Bioconductor", as.character(BiocManager::version())),
    paste("DESeq2", as.character(packageVersion("DESeq2"))),
    paste("tidyverse", as.character(packageVersion("tidyverse"))),
    paste("ggplot2", as.character(packageVersion("ggplot2"))),
    paste("pheatmap", as.character(packageVersion("pheatmap"))),
    paste("Quarto", system2("quarto", "--version", stdout = TRUE))
  )
  writeLines(session, file.path(stage, "session_info.txt"))
  report_dir <- file.path(stage, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  staged_qmd <- file.path(report_dir, basename(report_qmd))
  file.copy(report_qmd, staged_qmd, overwrite = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(report_dir)
  status <- system2("quarto", c(
    "render", basename(report_qmd), "--to", "html",
    "--output", "rnaseq_exploratory_qc.html",
    "-P", "data_root:.."
  ))
  unlink(staged_qmd)
  setwd(old_wd)
  if (!identical(status, 0L)) stop("Quarto rendering failed with status ", status)
  invisible(results)
}

