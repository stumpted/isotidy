# When I get a chance - modify the EA method of blank correction to accommodate for multiple caps worth of blanks in some blank runs, and divide by the number of blanks put in.

# Category 4: Blank subtraction and contamination removal
#
# Functions:
#   - calculate_blank_stats(): Summarize blank measurements
#   - apply_blank_correction(): Subtract blank effects (multiple methods)
#
# Dependencies:
# Used by: All pipelines
# Supports multiple methods: "ea_area_subtraction", "mean_subtraction", "normalized"
# ============================================================================

# ============================================================================

#' Calculate blank statistics
#'
#' Summarizes the blank measurements to determine correction factors.
#' Different correction methods may use different statistics.
#'
#' @param blank_df Dataframe containing ONLY blank measurements
#' @param n_blanks Number of blank replicates to average
#'   If NULL, uses all rows in blank_df
#' @param verbose Print statistics (default: TRUE)
#'
#' @return List containing:
#'   - mean_area: mean area_or_voltage of blanks
#'   - sd_area: standard deviation of blank area
#'   - mean_delta: mean delta_value of blanks
#'   - sd_delta: standard deviation of blank delta
#'   - n_blanks: number of blanks used
#'
#' @details
#' These statistics are used by apply_blank_correction()
#' Different instruments may average blanks differently:
#' - EA: average within a run, apply to that run
#' - TC/EA: average all blanks, apply globally
#' - GC/IRMS: average per compound, apply to that compound
#'
calculate_blank_stats <- function(blank_df, n_blanks = NULL, verbose = TRUE) {

  if (verbose) cat("Calculating blank statistics\n")

  if (nrow(blank_df) == 0) {
    stop("blank_df has no rows. Did you filter correctly?")
  }

  # Use number provided, or use all rows
  if (is.null(n_blanks)) {
    n_blanks <- nrow(blank_df)
  }

  # Validate blank data
  if (any(is.na(blank_df$area_or_voltage))) {
    warning("Some blank area measurements are NA")
  }

  if (any(is.na(blank_df$delta_value))) {
    warning("Some blank delta measurements are NA")
  }

  # Calculate statistics
  blank_stats <- list(
    mean_area = sum(blank_df$area_or_voltage, na.rm = TRUE) / n_blanks,
    sd_area = sd(blank_df$area_or_voltage, na.rm = TRUE),
    mean_delta = mean(blank_df$delta_value, na.rm = TRUE),
    sd_delta = sd(blank_df$delta_value, na.rm = TRUE),
    n_blanks = n_blanks
  )

  if (verbose) {
    cat("  Number of blanks:", blank_stats$n_blanks, "\n")
    cat("  Mean area:", round(blank_stats$mean_area, 2), "\n")
    cat("  Mean δ13C:", round(blank_stats$mean_delta, 4), "per mil\n")
    cat("  SD area:", round(blank_stats$sd_area, 2), "\n")
    cat("  SD δ13C:", round(blank_stats$sd_delta, 4), "per mil\n")
    cat("  ✓ Blank statistics calculated\n\n")
  }

  return(blank_stats)
}

# ============================================================================

#' Apply blank correction using specified method
#'
#' Removes the blank signal from all sample measurements.
#' Different instruments use different methods.
#'
#' @param df Canonical dataframe (all data)
#' @param blank_stats Blank statistics (from calculate_blank_stats)
#' @param method Correction method to use:
#'   - "ea_area_subtraction" (default): Account for varying blank per run
#'     Formula: δ_corrected = (area*δ - blank_area*blank_δ) / (area - blank_area)
#'   - "mean_subtraction": Simple subtraction
#'     Formula: δ_corrected = δ - blank_δ
#'   - "normalized": Normalize by area ratio
#'     Formula: δ_corrected = δ - (blank_δ * blank_area/area)
#' @param verbose Print progress (default: TRUE)
#'
#' @return Same dataframe with new column: delta_value_blank_corrected
#'   This column replaces delta_value for subsequent corrections
#'   Original delta_value is preserved
#'
#' @details
#' EA METHOD ("ea_area_subtraction"):
#' The EA typically has multiple blanks run throughout the sample sequence.
#' Each blank produces a peak with area and δ value.
#' This method accounts for the fact that blank peak area varies.
#' Used when there are multiple blanks per run.
#'
#' MEAN SUBTRACTION METHOD:
#' Simplest method: just subtract mean blank δ from all samples.
#' Used when you have a single global blank value.
#' Assumes area doesn't matter.
#'
#' NORMALIZED METHOD:
#' Corrects for blank area relative to sample area.
#' Intermediate between EA and mean subtraction.
#'
apply_blank_correction <- function(df, blank_stats, method = "ea_area_subtraction",
                                   verbose = TRUE) {

  if (verbose) cat("Applying blank correction\n")
  if (verbose) cat("  Method:", method, "\n")

  # Validate method
  valid_methods <- c("ea_area_subtraction", "mean_subtraction", "normalized")
  if (!method %in% valid_methods) {
    stop(paste("Unknown method. Valid options:",
               paste(valid_methods, collapse = ", ")))
  }

  # Apply correction based on method
  if (method == "ea_area_subtraction") {
    # EA method: accounts for different blank areas
    # Formula: (area*delta - blank_area*blank_delta) / (area - blank_area)

    df$delta_value_blank_corrected <- (
      (df$area_or_voltage * df$delta_value -
       blank_stats$mean_area * blank_stats$mean_delta) /
      (df$area_or_voltage - blank_stats$mean_area)
    )

  } else if (method == "mean_subtraction") {
    # Simple: subtract mean blank delta from all samples
    df$delta_value_blank_corrected <- df$delta_value - blank_stats$mean_delta

  } else if (method == "normalized") {
    # Normalize by area ratio: subtract (blank_delta * area_ratio)
    area_ratio <- blank_stats$mean_area / df$area_or_voltage
    df$delta_value_blank_corrected <- df$delta_value - (blank_stats$mean_delta * area_ratio)
  }

  if (verbose) {
    # Compare before and after
    cat("  Before correction: δ range",
        round(min(df$delta_value, na.rm = TRUE), 2), "to",
        round(max(df$delta_value, na.rm = TRUE), 2), "per mil\n")

    cat("  After correction: δ range",
        round(min(df$delta_value_blank_corrected, na.rm = TRUE), 2), "to",
        round(max(df$delta_value_blank_corrected, na.rm = TRUE), 2), "per mil\n")

    cat("  ✓ Blank correction applied\n\n")
  }

  return(df)
}

