source("scripts/utils_io.R")
suppressPackageStartupMessages(library(dplyr))

required_manifest_columns <- c(
  "sample_id", "species", "bioproject", "tissue", "condition",
  "replicate", "layout", "r1", "r2", "salmon_index", "canary"
)

# require_columns is provided by utils_io.R

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
