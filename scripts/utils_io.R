# Shared I/O utilities for the hsp_hsf pipeline.
# Source this file from any R script that reads or writes TSV data.

#' Read a TSV file with all columns as character (no automatic type coercion).
#'
#' This is the standard reader for manifests, design tables, and metadata files
#' where preserving exact string values (including empty strings) is essential.
#'
#' @param path File path to read.
#' @return A data.frame with character columns and no NA injection.
read_character_tsv <- function(path) {
  read.delim(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character", na.strings = character()
  )
}

#' Read a TSV file with all columns as character, suppressing quote parsing.
#'
#' Identical to \code{read_character_tsv} but also sets \code{quote = ""} so
#' that embedded quotation marks in metadata values are treated as literal text.
#'
#' @param path File path to read.
#' @return A data.frame with character columns and no NA injection.
read_metadata_tsv <- function(path) {
  read.delim(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character", na.strings = character(), quote = ""
  )
}

#' Validate that a data frame or named list contains the expected columns.
#'
#' Raises an error listing the missing columns if any are absent.
#'
#' @param x     A data.frame or named object.
#' @param required Character vector of required column names.
#' @param label   Human-readable label for error messages (e.g. "manifest").
require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

#' Validate that character columns in sample metadata contain no NA or unsafe
#' characters (tabs, carriage returns, newlines).
#'
#' @param sample_rows A data.frame of sample metadata.
validate_metadata_values <- function(sample_rows) {
  character_columns <- lapply(sample_rows, as.character)
  if (any(vapply(character_columns, anyNA, logical(1)))) {
    stop("sample metadata contains missing metadata values; use an empty string when intentional")
  }
  unsafe <- vapply(
    character_columns,
    function(values) any(grepl("[\t\r\n]", values, perl = TRUE)),
    logical(1)
  )
  if (any(unsafe)) stop("sample metadata values may not contain tabs or newlines")
  invisible(TRUE)
}

#' Write a data frame as a TSV file with safety checks.
#'
#' Raises an error if the data frame contains NA values or if character columns
#' contain tab, carriage-return, or newline characters. Parent directories are
#' created automatically.
#'
#' @param x    A data.frame to serialize.
#' @param path Output file path.
write_tsv_safe <- function(x, path) {
  if (anyNA(x)) stop("Cannot serialize NA values to ", path)
  char_cols <- names(x)[vapply(x, is.character, logical(1))]
  if (length(char_cols) > 0L) {
    if (any(vapply(x[char_cols], function(column) {
      any(grepl("[\t\r\n]", column))
    }, logical(1)))) stop("Unsafe TSV character in ", path)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
