# ============================================================================
# CORE IRMS FUNCTIONS - SCALE CALIBRATION
# ============================================================================
# Category 6: Conversion to reference scales (VPDB, NIST, PDB, etc.)
#
# Functions:
#   - build_scale_model(): Linear regression: reference_delta ~ measured_delta
#   - apply_scale_correction(): Convert to reference scale
#
# Dependencies:
# Used by: All pipelines
# Requires: Standards with all previous corrections applied
# Purpose: Convert machine-specific measurements to universal reference scale
# ============================================================================

# ============================================================================

#' Build scale calibration model
#'
#' Creates linear regression: certified_delta ~ measured_delta
#'
#' This converts your machine's measurements to a universal reference scale
#' (VPDB, NIST, PDB, etc.). Without this, your δ values are only useful
#' relative to themselves; with this, they're comparable to published data.
#'
#' @param corrected_stds_df Standards dataframe with all corrections applied
#'   Should have: delta_value_linearity_corrected (or similar)
#' @param stds_reference_df Reference standards database with certified values
#'   Must have: std.d13C.VPDB (or appropriate reference column)
#' @param delta_column Which corrected delta column to use
#'   Options: "delta_value", "delta_value_blank_corrected",
#'   "delta_value_linearity_corrected" (default: "delta_value_linearity_corrected")
#' @param reference_column Which reference column to use
#'   Options: "std.d13C.VPDB" (default), "std.d15N.AIR", "std.d34S.VCDT", etc.
#' @param config Configuration
#' @param verbose Print progress and model summary (default: TRUE)
#'
#' @return List containing:
#'   - model: lm object (linear regression)
#'   - slope: scale calibration slope
#'   - intercept: scale calibration intercept
#'   - r_squared: R² (model fit quality)
#'   - delta_column: which column was used
#'   - reference_column: which reference column was used
#'   - stds_data: the dataframe used to build model (for diagnostics)
#'
#' @details
#' Formula: reference_delta = slope * measured_delta + intercept
#'
#' Example interpretation:
#' If slope = 1.05, intercept = 0.5:
#' - Measured δ of 0 → reported δ of 0.5 ‰
#' - Measured δ of 10 → reported δ of 10.5 + 0.5 = 11 ‰
#' - Machine is biased (too low) and scaled differently than reference
#'
#' Quality check:
#' - Slope should be close to 1 (ideally 0.95-1.05)
#' - Intercept should be close to 0 (ideally < 1 ‰)
#' - R² should be high (> 0.99)
#'
build_scale_model <- function(
    stds_df,
    delta_column = "delta_value_linearity_corrected",
    reference_column = "delta_value_reference",
    config = NULL,
    verbose = TRUE
) {

  if (verbose) {
    cat("Building scale calibration model\n")
    cat("  Measured column:", delta_column, "\n")
    cat("  Reference column:", reference_column, "\n")
  }

  required_columns <- c(
    "sample_id",
    delta_column,
    reference_column
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(stds_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Standards data are missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  model_df <- stds_df %>%
    dplyr::select(
      sample_id,
      dplyr::all_of(c(reference_column, delta_column))
    ) %>%
    dplyr::filter(
      !is.na(.data[[reference_column]]),
      !is.na(.data[[delta_column]])
    )

  if (nrow(model_df) < 3) {
    stop(
      paste(
        "Need at least 3 matched standards, got",
        nrow(model_df)
      )
    )
  }

  if (verbose) {
    cat("  Standards used:", nrow(model_df), "\n")
  }

  # Build linear model: reference δ ~ measured δ
  formula <- stats::as.formula(
    paste(
      reference_column,
      "~",
      delta_column
    )
  )

  model <- stats::lm(
    formula,
    data = model_df
  )

  summary_coefs <- summary(model)$coefficients

  slope <- summary_coefs[delta_column, "Estimate"]
  intercept <- summary_coefs["(Intercept)", "Estimate"]

  model_summary <- summary(model)

  r_squared <- model_summary$r.squared
  rse <- model_summary$sigma

  if (verbose) {
    cat("  Slope:", round(slope, 6), "\n")
    cat("  Intercept:", round(intercept, 4), "per mil\n")
    cat("  R²:", round(r_squared, 4), "\n")
    cat("  RSE:", round(rse, 4), "per mil\n")

    if (abs(slope - 1) > 0.05) {
      warning(
        paste(
          "⚠ Slope deviates from 1:",
          round(slope, 4)
        )
      )
    }

    if (abs(intercept) > 1) {
      warning(
        paste(
          "⚠ Large intercept:",
          round(intercept, 4),
          "‰"
        )
      )
    }

    if (r_squared < 0.95) {
      warning(
        paste(
          "⚠ R² is low:",
          round(r_squared, 4)
        )
      )
    }

    cat("  ✓ Scale model built\n\n")
  }

  return(
    list(
      model = model,
      slope = slope,
      intercept = intercept,
      r_squared = r_squared,
      residual_std_error = rse,
      delta_column = delta_column,
      reference_column = reference_column,
      stds_data = model_df,
      n_stds = nrow(model_df)
    )
  )
}

# ============================================================================

#' Apply scale calibration
#'
#' Converts measured delta values to reference scale (VPDB, etc).
#' This is the final correction step, producing publishable results.
#'
#' @param df Canonical dataframe
#' @param scale_model Scale model object (from build_scale_model)
#' @param verbose Print progress (default: TRUE)
#'
#' @return Same dataframe with new column: delta_value_final
#'   This is the final, publishable δ value on reference scale.
#'
#' @details
#' Formula applied: δ_final = δ_measured * slope + intercept
#'
#' The result is what you report in publications.
#' This column is what goes into prepare_output().
#'
apply_scale_correction <- function(df, scale_model, verbose = TRUE) {

  if (verbose) cat("Applying scale correction\n")

  # Use the appropriate delta column from the model
  delta_col <- scale_model$delta_column
  ref_col <- scale_model$reference_column

  if (!delta_col %in% colnames(df)) {
    stop(paste("Delta column not found:", delta_col))
  }

  if (verbose) {
    cat("  Reference scale:", scale_model$reference_column, "\n")
  }

  # Apply scale correction
  df$delta_value_final <- (
    df[[delta_col]] * scale_model$slope + scale_model$intercept
  )

  if (verbose) {
    cat("  Before:",
        round(min(df[[delta_col]], na.rm = TRUE), 2), "to",
        round(max(df[[delta_col]], na.rm = TRUE), 2), "per mil\n")

    cat("  After:",
        round(min(df$delta_value_final, na.rm = TRUE), 2), "to",
        round(max(df$delta_value_final, na.rm = TRUE), 2), "per mil\n")

    cat("  ✓ Scale correction applied\n\n")
  }

  return(df)
}
