# Category 2: Standards identification and enrichment with reference data
#
# Functions:
#   - identify_standards(): Label which rows are standards
#   - extract_standards(): Isolate standards, merge with reference DB
#
# Dependencies: tidyverse
# Used by: All pipelines
# Requires: Canonical format (from core_00_import)
# ============================================================================

library(tidyverse)

identify_standards <- function(df, standard_names, verbose = TRUE) {
  
  if (verbose) cat("Identifying standards\n")
  
  if (verbose) {
    cat("  Standards to find:", paste(standard_names, collapse = ", "), "\n")
  }
  
  # Validate input
  if (!is.character(standard_names) || length(standard_names) == 0) {
    stop("standard_names must be non-empty character vector")
  }
  
  # Mark standards
  df <- df %>%
    mutate(is_standard = sample_id %in% standard_names)
  
  n_stds <- sum(df$is_standard, na.rm = TRUE)
  
  if (verbose) {
    cat("  Found", n_stds, "standard measurements\n")
    cat("  ✓ Standards identified\n\n")
  }
  
  return(df)
}

# ============================================================================

#' Extract standards and merge with reference values
#'
#' Isolates standards from dataset and enriches them with certified
#' isotope values from reference database. This merged dataframe is
#' used to build calibration curves.
#'
#' @param df Canonical dataframe with is_standard column
#' @param stds_reference_df Reference database with certified values
#'   Must have at least: standard_id, delta_value, element_fraction, etc.
#'   Other columns (std.d15N.AIR, standard_type, etc.) preserved
#' @param config Configuration (optional, for element, etc.)
#' @param verbose Print progress (default: TRUE)
#'
#' @return Dataframe containing ONLY standards with reference values merged
#'   Columns include:
#'   - All canonical columns (run_id, sample_id, delta_value, etc.)
#'   - Reference columns (delta_value, element_fraction, etc.)
#'   - Warnings if standards not found in reference DB
#'
#' @details
#' Merge is by:
#' - canonical: sample_id
#' - reference: standard_id
#'
#' A warning is issued for unmatched standards.
#' These unmatched rows will have NA in reference columns.
#' Consider this an error to investigate.
#'
extract_standards <- function(df, stds_reference_df, config = NULL, verbose = TRUE) {
  
  if (verbose) cat("Extracting standards\n")
  
  # Filter to standards only
  stds_df <- df %>%
    filter(is_standard == TRUE)
  
  if (verbose) {
    cat("  Standards to process:", nrow(stds_df), "measurements\n")
  }
  
  # Validate reference database
  if (!is.data.frame(stds_reference_df)) {
    stop("stds_reference_df must be a dataframe")
  }
  
  if (!"standard_id" %in% colnames(stds_reference_df)) {
    stop("stds_reference_df must have 'standard_id' column")
  }
  
  # Rename reference value before joining so measured and reference values are explicit
  stds_reference_df <- stds_reference_df %>%
    rename(delta_value_reference = delta_value)
  
  # Merge with reference database
  stds_df <- left_join(
    stds_df,
    stds_reference_df,
    by = c("sample_id" = "standard_id"),
    relationship = "many-to-one"
  )
  
  # Check for unmatched standards
  unmatched <- sum(is.na(stds_df$delta_value_reference))
  
  if (unmatched > 0) {
    unmatched_names <- stds_df %>%
      filter(is.na(delta_value_reference)) %>%
      pull(sample_id) %>%
      unique()
    
    warning(paste(
      unmatched, "standard measurements not found in reference database:",
      paste(unmatched_names, collapse = ", ")
    ))
  }
  
  if (verbose) {
    cat("  Merged with reference database\n")
    if (unmatched > 0) {
      cat("  ⚠ Warning:", unmatched, "unmatched standards\n")
    }
    cat("  ✓ Standards extracted\n\n")
  }
  
  return(stds_df)
}
