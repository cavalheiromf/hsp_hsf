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
