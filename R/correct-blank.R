# Blank handling: identify blanks, summarize them, subtract their effect.

# TODO: EA blank correction should handle multiple caps' worth of blanks per run, dividing by the number of blanks.

identify_blanks <- function(df, blank_identifier = "Blank", verbose = TRUE) {

  if (verbose) cat("Identifying blanks\n")

  # Combine multiple identifiers into one regex
  if (length(blank_identifier) > 1) {
    pattern_string <- paste(blank_identifier, collapse = "|")
  } else {
    pattern_string <- blank_identifier
  }

  df <- df %>%
    dplyr::mutate(
      is_blank = stringr::str_detect(.data$sample_id, pattern_string, negate = FALSE)
    )

  n_blanks <- sum(df$is_blank, na.rm = TRUE)

  if (verbose) {
    cat("  Found", n_blanks, "blank measurements\n")
    cat("  ✓ Blanks identified\n\n")
  }

  return(df)
}


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
calculate_blank_stats <- function(
    blank_df,
    verbose = TRUE
) {

  if (verbose) cat("Calculating blank statistics\n")

  if (!is.data.frame(blank_df) || nrow(blank_df) == 0) {
    stop("blank_df has no rows. Did you filter correctly?")
  }

  required_columns <- c(
    "sample_id",
    "area_or_voltage",
    "delta_value"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(blank_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "blank_df is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Capsule count is the number in the blank ID (default 1)
  blank_df <- blank_df %>%
    dplyr::mutate(
      capsule_count = as.numeric(
        stringr::str_extract(
          .data$sample_id,
          "\\d+"
        )
      ),
      capsule_count = dplyr::coalesce(
        .data$capsule_count,
        1
      ),
      normalized_area = .data$area_or_voltage / .data$capsule_count
    )

  if (any(is.na(blank_df$area_or_voltage))) {
    warning("Some blank area measurements are NA")
  }

  if (any(is.na(blank_df$delta_value))) {
    warning("Some blank delta measurements are NA")
  }

  # Stats use per-capsule normalized areas
  blank_stats <- list(
    mean_area = mean(
      blank_df$normalized_area,
      na.rm = TRUE
    ),
    sd_area = stats::sd(
      blank_df$normalized_area,
      na.rm = TRUE
    ),
    mean_delta = mean(
      blank_df$delta_value,
      na.rm = TRUE
    ),
    sd_delta = stats::sd(
      blank_df$delta_value,
      na.rm = TRUE
    ),
    n_blanks = nrow(blank_df)
  )

  if (verbose) {
    cat(
      "  Number of blank measurements:",
      blank_stats$n_blanks,
      "\n"
    )

    cat(
      "  Mean normalized area:",
      round(blank_stats$mean_area, 2),
      "\n"
    )

    cat(
      "  Mean δ13C:",
      round(blank_stats$mean_delta, 4),
      "per mil\n"
    )

    cat(
      "  SD normalized area:",
      round(blank_stats$sd_area, 2),
      "\n"
    )

    cat(
      "  SD δ13C:",
      round(blank_stats$sd_delta, 4),
      "per mil\n"
    )

    cat("  ✓ Blank statistics calculated\n\n")
  }

  return(blank_stats)
}


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

  valid_methods <- c("ea_area_subtraction", "mean_subtraction", "normalized")
  if (!method %in% valid_methods) {
    stop(paste("Unknown method. Valid options:",
               paste(valid_methods, collapse = ", ")))
  }

  if (method == "ea_area_subtraction") {

    df$delta_value_blank_corrected <- (
      (df$area_or_voltage * df$delta_value -
       blank_stats$mean_area * blank_stats$mean_delta) /
      (df$area_or_voltage - blank_stats$mean_area)
    )

  } else if (method == "mean_subtraction") {
    df$delta_value_blank_corrected <- df$delta_value - blank_stats$mean_delta

  } else if (method == "normalized") {
    area_ratio <- blank_stats$mean_area / df$area_or_voltage
    df$delta_value_blank_corrected <- df$delta_value - (blank_stats$mean_delta * area_ratio)
  }

  if (verbose) {
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
