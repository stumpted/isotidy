# Category 1: Removing noise/artifacts and filtering to relevant data
#
# Functions:
#   - clean_data(): Remove test runs, select peak/compound if needed
#   - identify_blanks(): Label blank samples
#
# Dependencies:
# Used by: All pipelines
# Requires: Canonical format (from core_00_import)

#' Clean canonical data: remove test runs, select peak/compound
#'
#' This function removes noise and artifacts from the dataset:
#' - Removes test/quality control runs (linearity, on-off, etc.)
#' - Optionally filters to specific peak number (EA)
#' - Optionally filters to specific compound (GC/LC)
#'
clean_data <- function(df, config, verbose = TRUE) {

  if (verbose) cat("Cleaning data\n")

  original_rows <- nrow(df)

  # Step 1: Remove test runs by pattern matching

  exclude_patterns <- config$exclude_patterns %||% c("on_off", "linearity", "test")

  if (!is.null(exclude_patterns) && "sample_id" %in% colnames(df) && length(exclude_patterns) > 0) {
    # Build pattern string: "pattern1|pattern2|pattern3"
    pattern_string <- paste(exclude_patterns, collapse = "|")

    # Count rows to be removed
    rows_to_remove <- df %>%
      dplyr::filter(stringr::str_detect(sample_id, pattern_string)) %>%
      nrow()

    # Remove test runs
    df <- df %>%
      dplyr::filter(!stringr::str_detect(sample_id, pattern_string))

    if (verbose && rows_to_remove > 0) {
      cat("  Removed test runs:", rows_to_remove, "rows\n")
    }
  }

  # Step 2: Filter to specific peak if specified (EA-like instruments)

  if (!is.null(config$peak_number) && "peak_number" %in% colnames(df)) {
    rows_before <- nrow(df)

    df <- df %>%
      dplyr::filter(peak_number == config$peak_number)

    rows_removed <- rows_before - nrow(df)

    if (verbose) {
      cat("  Selected peak", config$peak_number, ":",
          nrow(df), "rows remaining\n")
    }
  }

  # Step 3: Filter to specific compound if specified (GC/LC-like instruments)

  if (!is.null(config$select_compound) && "compound_name" %in% colnames(df)) {
    rows_before <- nrow(df)

    df <- df %>%
      dplyr::filter(compound_name == config$select_compound)

    rows_removed <- rows_before - nrow(df)

    if (verbose) {
      cat("  Selected compound", config$select_compound, ":",
          nrow(df), "rows remaining\n")
    }
  }

  if (verbose) {
    total_removed <- original_rows - nrow(df)
    cat("  Total removed:", total_removed, "rows\n")
    cat("  ✓ Data cleaned\n\n")
  }

  return(df)
}

identify_blanks <- function(df, blank_identifier = "Blank", verbose = TRUE) {

  if (verbose) cat("Identifying blanks\n")

  # Allow multiple patterns
  if (length(blank_identifier) > 1) {
    pattern_string <- paste(blank_identifier, collapse = "|")
  } else {
    pattern_string <- blank_identifier
  }

  # Mark blanks
  df <- df %>%
    dplyr::mutate(
      is_blank = stringr::str_detect(sample_id, pattern_string, negate = FALSE)
    )

  n_blanks <- sum(df$is_blank, na.rm = TRUE)

  if (verbose) {
    cat("  Found", n_blanks, "blank measurements\n")
    cat("  ✓ Blanks identified\n\n")
  }

  return(df)
}
