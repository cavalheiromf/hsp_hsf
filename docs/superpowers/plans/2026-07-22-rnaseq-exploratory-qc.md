# RNA-seq Exploratory QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible post-Salmon QC workflow that normalizes the experimental design, produces per-experiment exploratory diagnostics, documents estimable contrasts, and publishes validated HTML, TSV, PNG, and SVG outputs without excluding samples.

**Architecture:** Focused R libraries normalize metadata and calculate QC objects; thin CLI scripts expose deterministic production interfaces. A Quarto source consumes only published TSV and figure artifacts. The complete collection, including the rendered HTML, is staged and validated before same-filesystem promotion with rollback.

**Tech Stack:** R 4.5.1, Bioconductor DESeq2, tidyverse/readr/dplyr, ggplot2, pheatmap, Quarto, base R test scripts, Git.

## Global Constraints

- Differential-expression hypothesis testing is outside this plan.
- Analyze species independently; analyze wheat and soybean BioProjects independently for biological interpretation.
- Joint wheat and soybean PCA plots are batch visualizations only and cannot define biological contrasts.
- Preserve all 58 samples; flags are advisory and never trigger automatic removal.
- Preserve mapping warnings for `SRR39669466` and `SRR39669467`.
- Use estimated gene counts as input; round them only when constructing DESeq2 objects.
- Filter within each analysis group at count `>= 10` in at least the smallest biological replicate count.
- Use `DESeq2::varianceStabilizingTransformation(..., blind = TRUE)` for unsupervised diagnostics.
- Use explicit condition mappings for the six current BioProject analysis groups; reject unknown conditions.
- Write UTF-8 tab-separated tables with deterministic row and column order.
- Publish HTML, TSV, PNG, and SVG artifacts atomically; a failed rerun must preserve the previous valid collection.
- Use repository-relative defaults and record R, Bioconductor, package, and Quarto versions.

## File Structure

- Create `scripts/rnaseq_design_lib.R`: pure metadata parsing, design validation, reference levels, and contrast inventory.
- Create `scripts/build_rnaseq_design.R`: CLI that writes the normalized design table.
- Create `scripts/rnaseq_exploratory_qc_lib.R`: matrix validation, filtering, DESeq2 transformation, metrics, flags, figures, collection validation, and atomic promotion.
- Create `scripts/run_rnaseq_exploratory_qc.R`: production CLI and Quarto orchestration.
- Create `reports/rnaseq_exploratory_qc.qmd`: presentation-only report sourced from staged tables and figures.
- Create `tests/test_build_rnaseq_design.R`: metadata and contrast unit tests.
- Create `tests/test_rnaseq_exploratory_qc.R`: statistical, plotting, publication, and rollback tests.
- Create `tests/fixtures/rnaseq_qc/`: small deterministic manifest, QC, and count-matrix fixtures.
- Generate `config/rnaseq_design.tsv`.
- Generate `results/rnaseq/exploratory_qc/` and `reports/rnaseq_exploratory_qc.html`.

---

### Task 1: Normalize metadata and enumerate estimable contrasts

**Files:**
- Create: `scripts/rnaseq_design_lib.R`
- Create: `scripts/build_rnaseq_design.R`
- Create: `tests/test_build_rnaseq_design.R`
- Generate later: `config/rnaseq_design.tsv`

**Interfaces:**
- Consumes: `config/rnaseq_samples.tsv` with canonical sample metadata.
- Produces: `normalize_design(manifest) -> data.frame`, `validate_design(design) -> invisible(TRUE)`, `build_contrast_inventory(design) -> data.frame`, and CLI `build_rnaseq_design.R MANIFEST OUTPUT`.

- [ ] **Step 1: Write the failing metadata tests**

Create `tests/test_build_rnaseq_design.R` with the following test harness and assertions:

```r
source("scripts/rnaseq_design_lib.R")

assert_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}

manifest <- read.delim(
  "config/rnaseq_samples.tsv", stringsAsFactors = FALSE,
  check.names = FALSE, colClasses = "character"
)
design <- normalize_design(manifest)

stopifnot(nrow(design) == 58L)
stopifnot(identical(design$sample_id, manifest$sample_id))
stopifnot(identical(
  sort(unique(design$analysis_group)),
  sort(c(
    "setaria_prjna1496374", "sorghum_prjna742236",
    "wheat_prjna1249400", "wheat_prjna793265",
    "soybean_prjna295411", "soybean_prjna1033144"
  ))
))

expected_group_sizes <- c(
  setaria_prjna1496374 = 18L,
  sorghum_prjna742236 = 8L,
  wheat_prjna1249400 = 12L,
  wheat_prjna793265 = 6L,
  soybean_prjna295411 = 8L,
  soybean_prjna1033144 = 6L
)
observed_group_sizes <- table(design$analysis_group)
stopifnot(all(observed_group_sizes[names(expected_group_sizes)] == expected_group_sizes))

setaria <- design[design$analysis_group == "setaria_prjna1496374", ]
stopifnot(identical(sort(unique(setaria$co2)), c("ambient", "elevated")))
stopifnot(identical(
  sort(unique(setaria$dehydration_stage)),
  c("control", "cycle_1", "cycle_3")
))
stopifnot(all(table(setaria$co2, setaria$dehydration_stage) == 3L))

wheat_factorial <- design[design$analysis_group == "wheat_prjna1249400", ]
stopifnot(all(table(wheat_factorial$co2, wheat_factorial$water) == 3L))
soy_factorial <- design[design$analysis_group == "soybean_prjna295411", ]
stopifnot(all(table(soy_factorial$co2, soy_factorial$water) == 2L))

bad <- manifest
bad$condition[1] <- "unrecognized treatment"
assert_error(normalize_design(bad), "Unrecognized condition")

duplicate <- design
duplicate$replicate[2] <- duplicate$replicate[1]
duplicate$design_cell[2] <- duplicate$design_cell[1]
assert_error(validate_design(duplicate), "Duplicate replicate")

contrasts <- build_contrast_inventory(design)
stopifnot(all(c(
  "contrast_id", "species", "bioproject", "analysis_group",
  "model_formula", "effect_type", "numerator", "denominator",
  "reference_levels", "estimable", "notes"
) %in% names(contrasts)))
stopifnot(all(contrasts$estimable == "true"))
stopifnot(any(contrasts$contrast_id == "setaria_co2_x_cycle_1"))
stopifnot(any(contrasts$contrast_id == "wheat_prjna1249400_co2_x_drought"))
stopifnot(any(contrasts$contrast_id == "soybean_prjna295411_co2_x_air_drying"))
stopifnot(any(contrasts$contrast_id == "sorghum_elevated_vs_ambient"))
```

- [ ] **Step 2: Run the test and confirm RED**

Run:

```bash
Rscript tests/test_build_rnaseq_design.R
```

Expected: failure because `scripts/rnaseq_design_lib.R` does not exist.

- [ ] **Step 3: Implement explicit condition normalization**

Create `scripts/rnaseq_design_lib.R` with these public functions and exact factor vocabulary:

```r
suppressPackageStartupMessages(library(dplyr))

required_manifest_columns <- c(
  "sample_id", "species", "bioproject", "tissue", "condition",
  "replicate", "layout", "r1", "r2", "salmon_index", "canary"
)

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

parse_one_condition <- function(species, bioproject, condition) {
  result <- list(co2 = NA_character_, water = NA_character_,
                 dehydration_stage = NA_character_, replicate = NA_integer_)
  if (species == "setaria_viridis" && bioproject == "PRJNA1496374") {
    pattern <- "^(ambient|elevated) CO2, (control well-watered|dehydration cycle [13]), biological replicate ([123])$"
    match <- regexec(pattern, condition)
    pieces <- regmatches(condition, match)[[1]]
    if (length(pieces)) {
      result$co2 <- pieces[2]
      result$water <- if (pieces[3] == "control well-watered") "well_watered" else "dehydration"
      result$dehydration_stage <- c(
        "control well-watered" = "control",
        "dehydration cycle 1" = "cycle_1",
        "dehydration cycle 3" = "cycle_3"
      )[[pieces[3]]]
      result$replicate <- as.integer(pieces[4])
    }
  } else if (species == "sorghum_bicolor" && bioproject == "PRJNA742236") {
    pieces <- regmatches(condition, regexec("^(Elevated|Ambient) CO2, drought, rep([1-4])$", condition))[[1]]
    if (length(pieces)) {
      result$co2 <- tolower(pieces[2])
      result$water <- "drought"
      result$replicate <- as.integer(pieces[3])
    }
  } else if (species == "triticum_aestivum" && bioproject == "PRJNA1249400") {
    pieces <- regmatches(condition, regexec("^Alana_(elevated-CO2\\+drought|elevated-CO2|drought|control)_biol_rep([123])$", condition))[[1]]
    if (length(pieces)) {
      result$co2 <- if (grepl("elevated-CO2", pieces[2], fixed = TRUE)) "elevated" else "ambient"
      result$water <- if (grepl("drought", pieces[2], fixed = TRUE)) "drought" else "well_watered"
      result$replicate <- as.integer(pieces[3])
    }
  } else if (species == "triticum_aestivum" && bioproject == "PRJNA793265") {
    elevated <- regmatches(condition, regexec("^Elevated CO2\\. Rep([123])$", condition))[[1]]
    drought <- regmatches(condition, regexec("^Elevated CO2 and drought Rep([123])$", condition))[[1]]
    if (length(elevated)) {
      result$co2 <- "elevated"; result$water <- "well_watered"
      result$replicate <- as.integer(elevated[2])
    } else if (length(drought)) {
      result$co2 <- "elevated"; result$water <- "drought"
      result$replicate <- as.integer(drought[2])
    }
  } else if (species == "glycine_max" && bioproject == "PRJNA295411") {
    pieces <- regmatches(condition, regexec("^(Elevated|Ambient) CO2, (air-drying 50min|well-watered), rep([12])$", condition))[[1]]
    if (length(pieces)) {
      result$co2 <- tolower(pieces[2])
      result$water <- if (pieces[3] == "air-drying 50min") "air_drying" else "well_watered"
      result$replicate <- as.integer(pieces[4])
    }
  } else if (species == "glycine_max" && bioproject == "PRJNA1033144") {
    pieces <- regmatches(condition, regexec("^Elevated CO2 concentration \\+ (Drought stress|Controlled water treatment) ([123])$", condition))[[1]]
    if (length(pieces)) {
      result$co2 <- "elevated"
      result$water <- if (pieces[2] == "Drought stress") "drought" else "well_watered"
      result$replicate <- as.integer(pieces[3])
    }
  }
  if (is.na(result$replicate)) {
    stop("Unrecognized condition for ", species, "/", bioproject, ": ", condition)
  }
  result
}

analysis_group_for <- function(species, bioproject) {
  keys <- c(
    "setaria_viridis|PRJNA1496374" = "setaria_prjna1496374",
    "sorghum_bicolor|PRJNA742236" = "sorghum_prjna742236",
    "triticum_aestivum|PRJNA1249400" = "wheat_prjna1249400",
    "triticum_aestivum|PRJNA793265" = "wheat_prjna793265",
    "glycine_max|PRJNA295411" = "soybean_prjna295411",
    "glycine_max|PRJNA1033144" = "soybean_prjna1033144"
  )
  unname(keys[paste(species, bioproject, sep = "|")])
}

normalize_design <- function(manifest) {
  require_columns(manifest, required_manifest_columns, "manifest")
  if (anyDuplicated(manifest$sample_id)) stop("Duplicate sample_id in manifest")
  parsed <- lapply(seq_len(nrow(manifest)), function(i) {
    parse_one_condition(manifest$species[i], manifest$bioproject[i], manifest$condition[i])
  })
  design <- manifest %>%
    mutate(
      analysis_group = mapply(analysis_group_for, species, bioproject),
      co2 = vapply(parsed, `[[`, character(1), "co2"),
      water = vapply(parsed, `[[`, character(1), "water"),
      dehydration_stage = vapply(parsed, `[[`, character(1), "dehydration_stage"),
      replicate = as.character(vapply(parsed, `[[`, integer(1), "replicate")),
      design_cell = ifelse(
        analysis_group == "setaria_prjna1496374",
        paste(co2, dehydration_stage, sep = ":"),
        ifelse(grepl("prjna1249400|prjna295411", analysis_group),
               paste(co2, water, sep = ":"),
               ifelse(grepl("sorghum", analysis_group), co2, water))
      )
    )
  validate_design(design)
  design
}

validate_design <- function(design) {
  required <- c(required_manifest_columns, "analysis_group", "co2", "water",
                "dehydration_stage", "design_cell")
  require_columns(design, required, "design")
  if (anyNA(design$analysis_group) || any(!nzchar(design$analysis_group))) {
    stop("Unknown species/BioProject analysis group")
  }
  key <- paste(design$analysis_group, design$design_cell, design$replicate, sep = "|")
  if (anyDuplicated(key)) stop("Duplicate replicate within a design cell")
  expected <- c(
    setaria_prjna1496374 = 3L, sorghum_prjna742236 = 4L,
    wheat_prjna1249400 = 3L, wheat_prjna793265 = 3L,
    soybean_prjna295411 = 2L, soybean_prjna1033144 = 3L
  )
  observed <- table(design$analysis_group, design$design_cell)
  for (group in rownames(observed)) {
    cells <- observed[group, observed[group, ] > 0]
    if (any(cells != expected[[group]])) stop("Unexpected replicate count in ", group)
  }
  invisible(TRUE)
}
```

Append the complete contrast implementation to `scripts/rnaseq_design_lib.R`:

```r
contrast_row <- function(id, species, project, group, model, type,
                         numerator, denominator, references, notes) {
  data.frame(
    contrast_id = id, species = species, bioproject = project,
    analysis_group = group, model_formula = model, effect_type = type,
    numerator = numerator, denominator = denominator,
    reference_levels = references, estimable = "true", notes = notes,
    stringsAsFactors = FALSE
  )
}

build_contrast_inventory <- function(design) {
  validate_design(design)
  r <- contrast_row
  rows <- list(
    r("setaria_elevated_vs_ambient_control", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "elevated", "ambient", "co2=ambient;dehydration_stage=control", "CO2 effect at control"),
    r("setaria_elevated_vs_ambient_cycle_1", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "elevated", "ambient", "co2=ambient;dehydration_stage=control", "CO2 effect at cycle 1"),
    r("setaria_elevated_vs_ambient_cycle_3", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "elevated", "ambient", "co2=ambient;dehydration_stage=control", "CO2 effect at cycle 3"),
    r("setaria_cycle_1_vs_control_ambient", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "cycle_1", "control", "co2=ambient;dehydration_stage=control", "Cycle 1 effect at ambient CO2"),
    r("setaria_cycle_3_vs_control_ambient", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "cycle_3", "control", "co2=ambient;dehydration_stage=control", "Cycle 3 effect at ambient CO2"),
    r("setaria_cycle_1_vs_control_elevated", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "cycle_1", "control", "co2=ambient;dehydration_stage=control", "Cycle 1 effect at elevated CO2"),
    r("setaria_cycle_3_vs_control_elevated", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "simple_effect", "cycle_3", "control", "co2=ambient;dehydration_stage=control", "Cycle 3 effect at elevated CO2"),
    r("setaria_co2_x_cycle_1", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "interaction", "elevated:cycle_1", "ambient:cycle_1", "co2=ambient;dehydration_stage=control", "CO2 by cycle 1 interaction"),
    r("setaria_co2_x_cycle_3", "setaria_viridis", "PRJNA1496374", "setaria_prjna1496374", "~ co2 * dehydration_stage", "interaction", "elevated:cycle_3", "ambient:cycle_3", "co2=ambient;dehydration_stage=control", "CO2 by cycle 3 interaction"),
    r("sorghum_elevated_vs_ambient", "sorghum_bicolor", "PRJNA742236", "sorghum_prjna742236", "~ co2", "main_effect", "elevated", "ambient", "co2=ambient", "CO2 effect under drought"),
    r("wheat_prjna1249400_elevated_vs_ambient_well_watered", "triticum_aestivum", "PRJNA1249400", "wheat_prjna1249400", "~ co2 * water", "simple_effect", "elevated", "ambient", "co2=ambient;water=well_watered", "CO2 effect when well watered"),
    r("wheat_prjna1249400_elevated_vs_ambient_drought", "triticum_aestivum", "PRJNA1249400", "wheat_prjna1249400", "~ co2 * water", "simple_effect", "elevated", "ambient", "co2=ambient;water=well_watered", "CO2 effect under drought"),
    r("wheat_prjna1249400_drought_vs_well_watered_ambient", "triticum_aestivum", "PRJNA1249400", "wheat_prjna1249400", "~ co2 * water", "simple_effect", "drought", "well_watered", "co2=ambient;water=well_watered", "Drought effect at ambient CO2"),
    r("wheat_prjna1249400_drought_vs_well_watered_elevated", "triticum_aestivum", "PRJNA1249400", "wheat_prjna1249400", "~ co2 * water", "simple_effect", "drought", "well_watered", "co2=ambient;water=well_watered", "Drought effect at elevated CO2"),
    r("wheat_prjna1249400_co2_x_drought", "triticum_aestivum", "PRJNA1249400", "wheat_prjna1249400", "~ co2 * water", "interaction", "elevated:drought", "ambient:drought", "co2=ambient;water=well_watered", "CO2 by drought interaction"),
    r("wheat_prjna793265_drought_vs_well_watered", "triticum_aestivum", "PRJNA793265", "wheat_prjna793265", "~ water", "main_effect", "drought", "well_watered", "water=well_watered", "Drought effect under elevated CO2"),
    r("soybean_prjna295411_elevated_vs_ambient_well_watered", "glycine_max", "PRJNA295411", "soybean_prjna295411", "~ co2 * water", "simple_effect", "elevated", "ambient", "co2=ambient;water=well_watered", "CO2 effect when well watered"),
    r("soybean_prjna295411_elevated_vs_ambient_air_drying", "glycine_max", "PRJNA295411", "soybean_prjna295411", "~ co2 * water", "simple_effect", "elevated", "ambient", "co2=ambient;water=well_watered", "CO2 effect during air drying"),
    r("soybean_prjna295411_air_drying_vs_well_watered_ambient", "glycine_max", "PRJNA295411", "soybean_prjna295411", "~ co2 * water", "simple_effect", "air_drying", "well_watered", "co2=ambient;water=well_watered", "Air-drying effect at ambient CO2"),
    r("soybean_prjna295411_air_drying_vs_well_watered_elevated", "glycine_max", "PRJNA295411", "soybean_prjna295411", "~ co2 * water", "simple_effect", "air_drying", "well_watered", "co2=ambient;water=well_watered", "Air-drying effect at elevated CO2"),
    r("soybean_prjna295411_co2_x_air_drying", "glycine_max", "PRJNA295411", "soybean_prjna295411", "~ co2 * water", "interaction", "elevated:air_drying", "ambient:air_drying", "co2=ambient;water=well_watered", "CO2 by air-drying interaction"),
    r("soybean_prjna1033144_drought_vs_well_watered", "glycine_max", "PRJNA1033144", "soybean_prjna1033144", "~ water", "main_effect", "drought", "well_watered", "water=well_watered", "Drought effect under elevated CO2")
  )
  inventory <- bind_rows(rows)
  observed_groups <- unique(design$analysis_group)
  if (!all(inventory$analysis_group %in% observed_groups)) stop("Contrast references an absent analysis group")
  inventory
}
```

- [ ] **Step 4: Add the CLI**

Create executable `scripts/build_rnaseq_design.R`:

```r
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
```

- [ ] **Step 5: Run focused and regression tests**

Run:

```bash
Rscript tests/test_build_rnaseq_design.R
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
```

Expected: all tests pass, including 38 existing Python tests.

- [ ] **Step 6: Commit Task 1**

```bash
git add scripts/rnaseq_design_lib.R scripts/build_rnaseq_design.R tests/test_build_rnaseq_design.R
git commit -m "feat: normalize RNA-seq experimental design"
```

---

### Task 2: Calculate group-specific exploratory QC objects

**Files:**
- Create: `scripts/rnaseq_exploratory_qc_lib.R`
- Create: `tests/test_rnaseq_exploratory_qc.R`
- Create: `tests/fixtures/rnaseq_qc/manifest.tsv`
- Create: `tests/fixtures/rnaseq_qc/quantification.tsv`
- Create: `tests/fixtures/rnaseq_qc/matrices/test_species/gene_counts.tsv`

**Interfaces:**
- Consumes: normalized design, consolidated Salmon QC, and gene-count matrices.
- Produces: `read_gene_counts(path)`, `validate_qc_inputs(design, salmon_qc, counts_by_species)`, `calculate_group_qc(counts, metadata, salmon_qc) -> list`, and long-form metrics/PCA/correlation/distance tables.

- [ ] **Step 1: Create a deterministic synthetic fixture**

Create a six-sample fixture with two design cells and three replicates per cell. Use genes `gene_a` through `gene_f`; make `gene_f` low in every sample, and make sample `s6` lower in library size without making it invalid. The count matrix header must be `feature_id`, followed by `s1` through `s6`. Include valid Salmon QC rows and make `s6` a mapping warning at `65.0` percent.

Create `tests/fixtures/rnaseq_qc/manifest.tsv` with these tab-separated rows:

```text
sample_id	species	bioproject	condition	replicate	analysis_group	co2	water	dehydration_stage	design_cell
s1	test_species	TEST1	ambient control rep1	1	fixture_group	ambient	well_watered		ambient:well_watered
s2	test_species	TEST1	ambient control rep2	2	fixture_group	ambient	well_watered		ambient:well_watered
s3	test_species	TEST1	ambient control rep3	3	fixture_group	ambient	well_watered		ambient:well_watered
s4	test_species	TEST1	elevated drought rep1	1	fixture_group	elevated	drought		elevated:drought
s5	test_species	TEST1	elevated drought rep2	2	fixture_group	elevated	drought		elevated:drought
s6	test_species	TEST1	elevated drought rep3	3	fixture_group	elevated	drought		elevated:drought
```

Create `tests/fixtures/rnaseq_qc/quantification.tsv`:

```text
sample_id	species	status	mapping_flag	percent_mapped	library_types
s1	test_species	pass	pass	85.0	IU
s2	test_species	pass	pass	84.0	IU
s3	test_species	pass	pass	86.0	IU
s4	test_species	pass	pass	82.0	IU
s5	test_species	pass	pass	83.0	IU
s6	test_species	pass	warning	65.0	IU
```

The exact matrix values are:

```text
feature_id s1 s2 s3 s4 s5 s6
gene_a 100 110 120 300 320 150
gene_b 80 90 100 240 250 120
gene_c 300 290 310 100 110 50
gene_d 40 45 50 80 90 40
gene_e 20 25 30 35 40 20
gene_f 0 1 0 1 0 1
```

Write the actual fixture with tabs, not spaces.

- [ ] **Step 2: Write failing QC calculation tests**

Extend `tests/test_rnaseq_exploratory_qc.R`:

```r
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
```

- [ ] **Step 3: Run the test and confirm RED**

Run `Rscript tests/test_rnaseq_exploratory_qc.R`.

Expected: failure because `scripts/rnaseq_exploratory_qc_lib.R` does not exist.

- [ ] **Step 4: Implement matrix validation, filtering, and transformation**

Create `scripts/rnaseq_exploratory_qc_lib.R` with:

```r
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

robust_flag <- function(values, direction = c("both", "low", "high")) {
  direction <- match.arg(direction)
  if (length(values) < 6L) return(rep(FALSE, length(values)))
  center <- median(values)
  spread <- mad(values, constant = 1.4826)
  if (!is.finite(spread) || spread == 0) return(rep(FALSE, length(values)))
  low <- values < center - 3 * spread
  high <- values > center + 3 * spread
  switch(direction, both = low | high, low = low, high = high)
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
  flags <- data.frame(
    sample_id = salmon_qc$sample_id[salmon_qc$mapping_flag != "pass"],
    metric = "percent_mapped",
    observed_value = salmon_qc$percent_mapped[salmon_qc$mapping_flag != "pass"],
    rule = "mapping_flag != pass", severity = "warning",
    explanation = "Preserved Salmon mapping warning", stringsAsFactors = FALSE
  )
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
```

Add this exact advisory helper and replace the initial `flags` object with `build_advisory_flags(sample_metrics, pca$x, correlation)` followed by the mapping-warning rows:

```r
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
```

Append the Salmon rows with `bind_rows()` and assert that every metadata sample remains in `sample_metrics`, regardless of flags.

- [ ] **Step 5: Run focused and full tests**

```bash
Rscript tests/test_rnaseq_exploratory_qc.R
Rscript tests/test_build_rnaseq_design.R
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
```

Expected: all tests pass.

- [ ] **Step 6: Commit Task 2**

```bash
git add scripts/rnaseq_exploratory_qc_lib.R tests/test_rnaseq_exploratory_qc.R tests/fixtures/rnaseq_qc
git commit -m "feat: calculate exploratory RNA-seq QC"
```

---

### Task 3: Generate deterministic tables and figures

**Files:**
- Modify: `scripts/rnaseq_exploratory_qc_lib.R`
- Modify: `tests/test_rnaseq_exploratory_qc.R`

**Interfaces:**
- Consumes: the `calculate_group_qc()` result and a normalized design.
- Produces: `write_group_tables()`, `write_group_figures()`, `write_combined_tables()`, and exactly named PNG/SVG artifacts.

- [ ] **Step 1: Write failing output-inventory tests**

Add to `tests/test_rnaseq_exploratory_qc.R`:

```r
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
```

- [ ] **Step 2: Run the test and confirm RED**

Run `Rscript tests/test_rnaseq_exploratory_qc.R`.

Expected: failure because `write_group_tables()` is undefined.

- [ ] **Step 3: Implement stable table writers**

Add these table writers. They create parent directories, reject missing or unsafe character values, and preserve input order:

```r
write_tsv_safe <- function(x, path) {
  character_columns <- vapply(x, is.character, logical(1))
  if (anyNA(x)) stop("Cannot serialize NA values to ", path)
  if (any(vapply(x[character_columns], function(column) {
    any(grepl("[\\t\\r\\n]", column))
  }, logical(1)))) stop("Unsafe TSV character in ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

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
```

The resulting collection-level files are:

```text
tables/sample_metrics.tsv
tables/gene_filter_summary.tsv
tables/pca_scores.tsv
tables/advisory_flags.tsv
contrasts.tsv
```

Empty advisory output must retain its header.

- [ ] **Step 4: Implement paired PNG/SVG plotting**

Add:

```r
save_plot_pair <- function(plot, stem, width = 8, height = 5) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         units = "in", dpi = 150, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height,
         units = "in", bg = "white")
}
```

Add the complete group plot implementation:

```r
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
```

Use `sample_id` as the stable label and keep the metadata sample order in both heatmaps. The fixed palette and explicit devices ensure PNG/SVG parity.

- [ ] **Step 5: Run tests and inspect fixture graphics**

```bash
Rscript tests/test_rnaseq_exploratory_qc.R
find /tmp -path '*rnaseq-qc-output-*' -type f 2>/dev/null | head
python3 -m unittest discover -s tests -v
```

Expected: R and Python tests pass; the fixture test verifies exactly 12 nonempty figure files.

- [ ] **Step 6: Commit Task 3**

```bash
git add scripts/rnaseq_exploratory_qc_lib.R tests/test_rnaseq_exploratory_qc.R
git commit -m "feat: render exploratory RNA-seq QC figures"
```

---

### Task 4: Render and atomically publish the QC collection

**Files:**
- Create: `scripts/run_rnaseq_exploratory_qc.R`
- Create: `reports/rnaseq_exploratory_qc.qmd`
- Modify: `scripts/rnaseq_exploratory_qc_lib.R`
- Modify: `tests/test_rnaseq_exploratory_qc.R`

**Interfaces:**
- Consumes: repository root, normalized design, Salmon QC, per-species gene counts, and report source.
- Produces: CLI `run_rnaseq_exploratory_qc.R DESIGN SALMON_QC MATRIX_ROOT OUTPUT_ROOT REPORT_QMD REPORT_HTML`, validated collection, and stable HTML link.

- [ ] **Step 1: Write failing staging and rollback tests**

Build a temporary fixture repository and add:

```r
old_root <- tempfile("published-qc-")
dir.create(old_root)
writeLines("old sentinel", file.path(old_root, "sentinel.txt"))
report_link <- tempfile("published-report-", fileext = ".html")
writeLines("old report", report_link)

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
```

Add a success test that checks the required tables, 12 fixture group figures, nonempty `session_info.txt`, nonempty internal HTML, and a `reports/rnaseq_exploratory_qc.html` symlink resolving to the HTML inside the promoted collection.

- [ ] **Step 2: Run the test and confirm RED**

Run `Rscript tests/test_rnaseq_exploratory_qc.R`.

Expected: failure because `publish_qc_collection()` is undefined.

- [ ] **Step 3: Create the Quarto presentation source**

Create `reports/rnaseq_exploratory_qc.qmd` with YAML parameters and sections:

```yaml
---
title: "RNA-seq exploratory quality control"
format:
  html:
    toc: true
    embed-resources: true
execute:
  echo: false
  warning: false
params:
  data_root: "../results/rnaseq/exploratory_qc"
---
```

Append this report body after the YAML:

````markdown
```{r}
root <- params$data_root
read_qc_table <- function(name) read.delim(
  file.path(root, name), stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = character()
)
metrics <- read_qc_table(file.path("tables", "sample_metrics.tsv"))
filters <- read_qc_table(file.path("tables", "gene_filter_summary.tsv"))
scores <- read_qc_table(file.path("tables", "pca_scores.tsv"))
flags <- read_qc_table(file.path("tables", "advisory_flags.tsv"))
contrasts <- read_qc_table("contrasts.tsv")
groups <- unique(metrics$analysis_group)
```

## Methods

Gene-level Salmon estimated counts were filtered independently within each analysis group. Genes required a count of at least 10 in at least the smallest replicate-group size. Counts were rounded only for DESeq2 input; size factors and a blind variance-stabilizing transformation were used for unsupervised diagnostics. No sample was removed and no differential-expression test was fitted.

## Dataset summary

```{r}
knitr::kable(as.data.frame.matrix(table(metrics$species, metrics$bioproject)),
             caption = "Samples by species and BioProject")
knitr::kable(filters, caption = "Gene filtering by analysis group")
```

## Salmon QC

```{r}
salmon_columns <- intersect(c("sample_id", "species", "percent_mapped",
                              "mapping_flag", "library_types"), names(metrics))
knitr::kable(metrics[, salmon_columns, drop = FALSE],
             caption = "Salmon mapping and library diagnostics")
```

## Analysis groups

```{r, results='asis'}
for (group in groups) {
  cat("\n### ", group, "\n\n", sep = "")
  stems <- c("library_size", "detected_genes", "vst_distribution", "pca",
             "correlation_heatmap", "distance_heatmap")
  paths <- file.path(root, "figures", paste0(group, "_", stems, ".png"))
  knitr::include_graphics(paths)
  cat("\n")
}
```

## Advisory flags

```{r}
if (nrow(flags)) knitr::kable(flags) else cat("No advisory flags were emitted.")
```

## Joint BioProject diagnostics

These plots visualize BioProject separation only; they do not define cross-project biological contrasts.

```{r}
knitr::include_graphics(file.path(
  root, "figures", c("wheat_joint_bioproject_pca.png",
                     "soybean_joint_bioproject_pca.png")
))
```

## Estimable contrasts

The following contrasts are documented for a later differential-expression workflow and were not fitted here.

```{r}
knitr::kable(contrasts)
```

## Reproducibility

```{r, results='asis'}
cat("```text\n")
cat(readLines(file.path(root, "session_info.txt")), sep = "\n")
cat("\n```\n")
```
````

The report reads staged outputs only and never recalculates normalization or PCA.

- [ ] **Step 4: Implement collection validation and atomic promotion**

Add `validate_qc_collection(stage, expected_design)` to require:

- the five collection-level TSV files;
- `session_info.txt`;
- six group-specific figure stems in PNG and SVG for every analysis group;
- wheat and soybean joint PCA in both formats;
- an internal `report/rnaseq_exploratory_qc.html`;
- exactly 58 unique samples in `sample_metrics.tsv` in design order;
- all six analysis groups;
- both accepted Setaria warning rows;
- finite numeric metrics and no silent sample exclusion.

Use this implementation:

```r
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
```

Implement `publish_qc_collection()` with this transaction:

```r
relative_path <- function(target, from) {
  target_parts <- strsplit(normalizePath(target, mustWork = FALSE), "/", fixed = TRUE)[[1]]
  from_parts <- strsplit(normalizePath(from, mustWork = TRUE), "/", fixed = TRUE)[[1]]
  common <- 0L
  limit <- min(length(target_parts), length(from_parts))
  while (common < limit && target_parts[common + 1L] == from_parts[common + 1L]) {
    common <- common + 1L
  }
  file.path(c(rep("..", length(from_parts) - common),
              target_parts[(common + 1L):length(target_parts)]))
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
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)

  build_qc_collection(inputs, stage, report_qmd)
  validate_qc_collection(stage, inputs$design)
  failure_hook(stage)

  if (dir.exists(output_root)) {
    output_backup <- tempfile(paste0(".", basename(output_root), ".backup-"), tmpdir = parent)
    if (!rename_fn(output_root, output_backup)) stop("Cannot back up QC collection")
  }
  if (file.exists(report_html) || nzchar(Sys.readlink(report_html))) {
    report_backup <- tempfile(".rnaseq-qc-report.backup-", tmpdir = dirname(report_html))
    if (!rename_fn(report_html, report_backup)) {
      if (!is.null(output_backup)) rename_fn(output_backup, output_root)
      stop("Cannot back up QC report link")
    }
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
```

Harden the implementation by checking every rename/cleanup result, confirming that staging and output share a filesystem, and preserving a backup path in the error message if rollback itself fails. The `rename_fn` injection is required for deterministic promotion-failure tests.

- [ ] **Step 5: Implement the thin production CLI**

Add the collection builder to `scripts/rnaseq_exploratory_qc_lib.R`:

```r
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
  dir.create(report_dir)
  status <- system2("quarto", c(
    "render", report_qmd, "--to", "html",
    "--output", "rnaseq_exploratory_qc.html",
    "--output-dir", report_dir,
    "-P", paste0("data_root:", normalizePath(stage))
  ))
  if (!identical(status, 0L)) stop("Quarto rendering failed with status ", status)
  invisible(results)
}
```

Create executable `scripts/run_rnaseq_exploratory_qc.R`:

```r
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
read_character_tsv <- function(path) read.delim(
  path, stringsAsFactors = FALSE, check.names = FALSE,
  colClasses = "character", na.strings = character()
)
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
```

Run signature:

```text
run_rnaseq_exploratory_qc.R DESIGN SALMON_QC MATRIX_ROOT OUTPUT_ROOT REPORT_QMD REPORT_HTML
```

`build_qc_collection()` reads only its `inputs` object and writes only under `stage`; it never reads the production collection that it is about to replace.

- [ ] **Step 6: Run focused, rendering, and full tests**

```bash
command -v quarto
Rscript -e 'stopifnot(requireNamespace("DESeq2"), requireNamespace("tidyverse"), requireNamespace("ggplot2"), requireNamespace("pheatmap"))'
Rscript tests/test_build_rnaseq_design.R
Rscript tests/test_rnaseq_exploratory_qc.R
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
Rscript -e 'parse(file="scripts/run_rnaseq_exploratory_qc.R")'
git diff --check
```

Expected: all commands exit zero; the Quarto test renders a nonempty HTML file.

- [ ] **Step 7: Commit Task 4**

```bash
git add scripts/rnaseq_exploratory_qc_lib.R scripts/run_rnaseq_exploratory_qc.R reports/rnaseq_exploratory_qc.qmd tests/test_rnaseq_exploratory_qc.R
git commit -m "feat: publish reproducible RNA-seq QC report"
```

---

### Task 5: Run production QC and record validated artifacts

**Files:**
- Generate: `config/rnaseq_design.tsv`
- Generate: `results/rnaseq/exploratory_qc/`
- Generate: `reports/rnaseq_exploratory_qc.html`
- Verify: all Task 1–4 files and existing quantification artifacts.

**Interfaces:**
- Consumes: complete final gene matrices and consolidated Salmon QC.
- Produces: the approved production report collection and final provenance commit.

- [ ] **Step 1: Run the full preflight suite**

```bash
git status --short
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
Rscript tests/test_build_rnaseq_design.R
Rscript tests/test_rnaseq_exploratory_qc.R
quarto --version
```

Expected: clean worktree and all tests pass.

- [ ] **Step 2: Build and validate the production design**

```bash
Rscript scripts/build_rnaseq_design.R \
  config/rnaseq_samples.tsv \
  config/rnaseq_design.tsv
Rscript -e '
  x <- read.delim("config/rnaseq_design.tsv", stringsAsFactors=FALSE, check.names=FALSE)
  stopifnot(nrow(x) == 58L, length(unique(x$analysis_group)) == 6L)
  print(table(x$analysis_group))
'
```

Expected group sizes: 18 Setaria, 8 sorghum, 12 and 6 wheat, and 8 and 6 soybean.

- [ ] **Step 3: Publish the production collection**

```bash
Rscript scripts/run_rnaseq_exploratory_qc.R \
  config/rnaseq_design.tsv \
  results/rnaseq/qc/quantification.tsv \
  results/rnaseq/matrices/final \
  results/rnaseq/exploratory_qc \
  reports/rnaseq_exploratory_qc.qmd \
  reports/rnaseq_exploratory_qc.html
```

Expected: the collection and HTML link are promoted only after validation.

- [ ] **Step 4: Perform production acceptance checks**

```bash
Rscript -e '
  root <- "results/rnaseq/exploratory_qc"
  metrics <- read.delim(file.path(root,"tables","sample_metrics.tsv"), check.names=FALSE)
  flags <- read.delim(file.path(root,"tables","advisory_flags.tsv"), check.names=FALSE)
  contrasts <- read.delim(file.path(root,"contrasts.tsv"), check.names=FALSE)
  stopifnot(nrow(metrics) == 58L, !anyDuplicated(metrics$sample_id))
  stopifnot(all(c("SRR39669466","SRR39669467") %in% flags$sample_id))
  stopifnot(all(contrasts$estimable == "true"))
  stopifnot(file.info("reports/rnaseq_exploratory_qc.html")$size > 0)
  print(table(metrics$species, metrics$bioproject))
'
find results/rnaseq -maxdepth 1 -type d -name '.exploratory_qc.*' -print
find reports -maxdepth 1 -name '.rnaseq-qc-report.*' -print
```

Expected: 58 metrics rows, both accepted warnings, no staging/backup paths, and a nonempty HTML report.

- [ ] **Step 5: Re-run regression tests after production generation**

```bash
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
Rscript tests/test_build_rnaseq_design.R
Rscript tests/test_rnaseq_exploratory_qc.R
git diff --check
git status --short
```

Expected: all tests pass; only intended scripts, tests, design, report source/HTML, and QC artifacts appear.

- [ ] **Step 6: Commit production artifacts**

```bash
git add \
  config/rnaseq_design.tsv \
  scripts/rnaseq_design_lib.R \
  scripts/build_rnaseq_design.R \
  scripts/rnaseq_exploratory_qc_lib.R \
  scripts/run_rnaseq_exploratory_qc.R \
  tests/test_build_rnaseq_design.R \
  tests/test_rnaseq_exploratory_qc.R \
  tests/fixtures/rnaseq_qc \
  reports/rnaseq_exploratory_qc.qmd \
  reports/rnaseq_exploratory_qc.html \
  results/rnaseq/exploratory_qc
git diff --cached --check
git commit -m "feat: add exploratory RNA-seq quality control"
```

Before committing, verify that no FASTQ, SRA, transcriptome FASTA, Salmon index, or per-sample quantification directory is staged.

---

## Final Review Gate

Request an independent review across the complete implementation range. The reviewer must verify specification coverage, all six experimental designs, no sample exclusion, exact preservation of the two accepted mapping warnings, statistical method correctness, report reproducibility, publication rollback, production dimensions, test evidence, and Git scope. Resolve every Critical or Important finding before integration.
