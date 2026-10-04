# Cleaning canonical data: remove test runs, select peak/compound.

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

  # Remove test runs by sample_id pattern

  exclude_patterns <- config$exclude_patterns %||% c("on_off", "linearity", "test")

  if (!is.null(exclude_patterns) && "sample_id" %in% colnames(df) && length(exclude_patterns) > 0) {
    pattern_string <- paste(exclude_patterns, collapse = "|")

    rows_to_remove <- df %>%
      dplyr::filter(stringr::str_detect(.data$sample_id, pattern_string)) %>%
      nrow()

    df <- df %>%
      dplyr::filter(!stringr::str_detect(.data$sample_id, pattern_string))

    if (verbose && rows_to_remove > 0) {
      cat("  Removed test runs:", rows_to_remove, "rows\n")
    }
  }

  # Filter to one peak (EA)

  if (!is.null(config$peak_number) && "peak_number" %in% colnames(df)) {
    rows_before <- nrow(df)

    df <- df %>%
      dplyr::filter(.data$peak_number == config$peak_number)

    rows_removed <- rows_before - nrow(df)

    if (verbose) {
      cat("  Selected peak", config$peak_number, ":",
          nrow(df), "rows remaining\n")
    }
  }

  # Filter to one compound (GC/LC)

  if (!is.null(config$select_compound) && "compound_name" %in% colnames(df)) {
    rows_before <- nrow(df)

    df <- df %>%
      dplyr::filter(.data$compound_name == config$select_compound)

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
