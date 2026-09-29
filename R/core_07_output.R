# ============================================================================
# CORE IRMS FUNCTIONS - OUTPUT & UTILITIES
# ============================================================================
# Category 7: Final output preparation and diagnostic utilities
#
# Functions:
#   - prepare_output(): Select columns, format for export
#   - summarize_calibrations(): Quality metrics summary
#
# Dependencies:
# Used by: All pipelines (final step)
# Purpose: Format results for saving and publication
# ============================================================================

# ============================================================================

#' Prepare output dataframe with final corrected values
#'
#' Selects specific columns and formats dataframe for export.
#' This is the last step before saving to Excel/CSV.
#'
#' @param df Fully corrected data dataframe
#'   Should have all correction columns: delta_value_final, etc.
#' @param output_cols Character vector of columns to include in output
#'   If NULL, uses defaults: run_id, sample_id, element, delta_value_final, amount
#'   You can specify any columns that exist in df
#' @param verbose Print progress (default: TRUE)
#'
#' @return Dataframe with only specified columns
#'   Ready to export to Excel/CSV/database
#'   Rows and data unchanged, only columns selected
#'
#' @details
#' Common output column sets:
#'
#' Minimal (just the results):
#' c("sample_id", "delta_value_final", "amount")
#'
#' Standard (with metadata):
#' c("run_id", "sample_id", "element", "delta_value_final", "amount")
#'
#' Detailed (with intermediate corrections):
#' c("run_id", "sample_id", "delta_value",
#'   "delta_value_blank_corrected", "delta_value_linearity_corrected",
#'   "delta_value_final", "amount")
#'
#' With calibration info:
#' c("run_id", "sample_id", "C.Percent", "delta_value_final",
#'   "area_or_voltage", "instrument")
#'
prepare_output <- function(df, output_cols = NULL, verbose = TRUE) {

  if (verbose) cat("Preparing output\n")

  # Set defaults if not specified
  if (is.null(output_cols)) {
    output_cols <- c("run_id", "sample_id", "element",
                     "delta_value_final", "amount")
  }

  if (verbose) {
    cat("  Output columns:", paste(output_cols, collapse = ", "), "\n")
  }

  # Validate that all requested columns exist
  missing_cols <- setdiff(output_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing output columns:", paste(missing_cols, collapse = ", "),
               "\nAvailable columns:", paste(colnames(df), collapse = ", ")))
  }

  # Select output columns
  output_df <- df[output_cols]

  if (verbose) {
    cat("  Rows:", nrow(output_df), "\n")
    cat("  Columns:", ncol(output_df), "\n")
    cat("  ✓ Output prepared\n\n")
  }

  return(output_df)
}

# ============================================================================

#' Summarize calibration quality metrics
#'
#' Creates a summary table of all calibration models built during processing.
#' Useful for assessing data quality and detecting problems.
#'
#' @param models List of calibration models
#'   Each model should have: slope, intercept, r_squared
#'   Typically: models$mass, models$linearity, models$scale
#'
#' @return Dataframe with one row per calibration type
#'   Columns: calibration, slope, intercept, r_squared, n_standards
#'   Ready to display or export
#'
#' @details
#' Useful for quality control:
#' - Check R² > 0.99 for all models
#' - Check slope ≈ 1 and intercept ≈ 0 for scale model
#' - Compare R² between runs
#'
#' Example output:
#' # A tibble: 3 × 5
#'   calibration  slope intercept r_squared n_stds
#'   <chr>        <dbl>     <dbl>     <dbl>  <int>
#' 1 mass        0.0234    -0.123     0.999      5
#' 2 linearity   0.0001    -0.052     0.876      5
#' 3 scale       1.0125     0.456     0.998      5
#'
summarize_calibrations <- function(models) {

  # Handle empty list
  if (length(models) == 0) {
    return(tibble(
      calibration = character(),
      slope = numeric(),
      intercept = numeric(),
      r_squared = numeric(),
      n_stds = integer()
    ))
  }

  # Build summary row by row
  summary_rows <- lapply(names(models), function(name) {
    model_info <- models[[name]]

    tibble(
      calibration = name,
      slope = round(model_info$slope %||% NA, 8),
      intercept = round(model_info$intercept %||% NA, 4),
      r_squared = round(model_info$r_squared %||% NA, 4),
      n_stds = model_info$n_stds %||% NA_integer_
    )
  })

  # Combine all rows
  cal_summary <- bind_rows(summary_rows)

  return(cal_summary)
}

# ============================================================================

#' Print a formatted summary of processing results
#'
#' Displays key information from the pipeline run in a readable format.
#'
#' @param results Results list from pipeline
#' @param show_calibrations Show calibration summary (default: TRUE)
#' @param show_blanks Show blank correction info (default: TRUE)
#'
#' @return Invisible - prints to console
#'
print_results_summary <- function(results, show_calibrations = TRUE,
                                  show_blanks = TRUE) {

  cat("\n")
  cat("════════════════════════════════════════════════════════\n")
  cat("PROCESSING RESULTS SUMMARY\n")
  cat("════════════════════════════════════════════════════════\n\n")

  # Experiment info
  if (!is.null(results$config)) {
    cat("EXPERIMENT:\n")
    cat("  Name:", results$config$experiment_name, "\n")
    cat("  File:", results$config$file_name, "\n")
    cat("  Instrument:", results$data$instrument[1], "\n\n")
  }

  # Output summary
  cat("OUTPUT:\n")
  cat("  Samples processed:", nrow(results$data), "\n")
  cat("  Columns:", ncol(results$data), "\n")
  cat("  Saved to:", results$config$output_file_path, "\n\n")

  # Blank correction info
  if (show_blanks && !is.null(results$blank_stats)) {
    cat("BLANK CORRECTION:\n")
    cat("  Mean blank area:", round(results$blank_stats$mean_area, 2), "\n")
    cat("  Mean blank δ:", round(results$blank_stats$mean_delta, 4), "‰\n\n")
  }

  # Calibration summary
  if (show_calibrations && !is.null(results$models)) {
    cat("CALIBRATION MODELS:\n")
    cal_summary <- summarize_calibrations(results$models)
    print(cal_summary, n = Inf)
    cat("\n")
  }

  # Quality checks
  cat("QUALITY CHECKS:\n")

  if (!is.null(results$models)) {
    # Check mass model
    if (!is.null(results$models$mass)) {
      r2 <- results$models$mass$r_squared
      status <- if (r2 > 0.99) "✓ Good" else if (r2 > 0.95) "⚠ OK" else "✗ Poor"
      cat("  Mass calibration R²:", round(r2, 4), status, "\n")
    }

    # Check linearity model
    if (!is.null(results$models$linearity)) {
      r2 <- results$models$linearity$r_squared
      status <- if (r2 > 0.95) "✓ Good" else if (r2 > 0.80) "⚠ OK" else "✗ Poor"
      cat("  Linearity calibration R²:", round(r2, 4), status, "\n")
    }

    # Check scale model
    if (!is.null(results$models$scale)) {
      r2 <- results$models$scale$r_squared
      slope <- results$models$scale$slope
      intercept <- results$models$scale$intercept

      slope_ok <- abs(slope - 1) < 0.05
      intercept_ok <- abs(intercept) < 1

      r2_status <- if (r2 > 0.99) "✓" else if (r2 > 0.95) "⚠" else "✗"
      slope_status <- if (slope_ok) "✓" else "⚠"
      intercept_status <- if (intercept_ok) "✓" else "⚠"

      cat("  Scale calibration R²:", round(r2, 4), r2_status, "\n")
      cat("    Slope:", round(slope, 6), slope_status, "\n")
      cat("    Intercept:", round(intercept, 4), intercept_status, "\n")
    }
  }

  cat("\n")
  cat("════════════════════════════════════════════════════════\n\n")

  return(invisible(NULL))
}

# ============================================================================
# End of core_07_output.R
# ============================================================================
