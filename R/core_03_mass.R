# Category 3: Elemental mass calibration (optional, instrument-dependent)
#
# Functions:
#   - prepare_mass_standards(): Calculate elemental amounts
#   - build_mass_model(): Calibrate amount against area using one standard
#   - apply_mass_correction(): Apply model to all samples
#
# Dependencies:
# Used by: EA, TC/EA (maybe)
# NOT used by: GC/IRMS, LC/IRMS

#' Calculate elemental amounts for standards
#'
#' Calculates the amount of the selected element contained in each standard.
#'
#' @param stds_df Standards dataframe from extract_standards().
#' @param element Element to calibrate (e.g. "C", "N", "S", "H").
#' @param config Configuration.
#' @param verbose Print progress.
#'
#' @return Standards dataframe with a new {element}.Amount column.
prepare_mass_standards <- function(
    stds_df,
    element,
    config = NULL,
    verbose = TRUE
) {

  if (verbose) {
    cat(
      "Preparing mass standards for element",
      element,
      "\n"
    )
  }

  required_columns <- c(
    "amount",
    "element_fraction"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(stds_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Mass standards are missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (any(
    !is.na(stds_df$element_fraction) &
    (
      stds_df$element_fraction < 0 |
      stds_df$element_fraction > 1
    )
  )) {
    stop(
      "'element_fraction' must contain values between 0 and 1."
    )
  }

  amount_col <- paste0(element, ".Amount")

  stds_df[[amount_col]] <-
    stds_df$amount * stds_df$element_fraction

  if (verbose) {

    valid_amounts <- stds_df[[amount_col]][
      !is.na(stds_df[[amount_col]])
    ]

    if (length(valid_amounts) > 0) {
      cat("  Elemental amounts calculated\n")
      cat(
        "  Range:",
        round(min(valid_amounts, na.rm = TRUE), 4),
        "-",
        round(max(valid_amounts, na.rm = TRUE), 4),
        "\n"
      )
    }

    cat("  ✓ Mass standards prepared\n\n")
  }

  return(stds_df)
}


#' Build mass calibration model
#'
#' Builds a mass calibration using the configured calibration standard
#' measured at multiple sample amounts.
#'
#' @param stds_df Standards dataframe from prepare_mass_standards().
#' @param element Element to calibrate.
#' @param config Configuration containing calibration_standards.
#' @param verbose Print progress and model summary.
#' @param create_plot Create diagnostic plot.
#'
#' @return List containing the fitted model, coefficients, statistics,
#'   calibration standard, and optional diagnostic plot.

build_mass_model <- function(
    stds_df,
    element,
    config = NULL,
    verbose = TRUE,
    create_plot = TRUE
) {

  if (verbose) {
    cat(
      "Building mass calibration model for",
      element,
      "\n"
    )
  }

  # Check configuration
  if (is.null(config)) {
    stop(
      "An experiment configuration is required."
    )
  }

  if (is.null(config$standards)) {
    stop(
      "Configuration must contain a 'standards' section."
    )
  }

  if (is.null(config$standards$calibration_standards)) {
    stop(
      "Configuration must contain ",
      "'standards$calibration_standards'."
    )
  }

  # Get the configured calibration standard
  calibration_standard <-
    config$standards$calibration_standards[[element]]

  if (
    is.null(calibration_standard) ||
    is.na(calibration_standard) ||
    calibration_standard == ""
  ) {
    stop(
      "No calibration standard configured for element ",
      element,
      ". Add it under ",
      "config$standards$calibration_standards."
    )
  }

  # Check standard data
  if (!is.data.frame(stds_df)) {
    stop(
      "'stds_df' must be a data frame."
    )
  }

  if (nrow(stds_df) == 0) {
    stop(
      "'stds_df' contains no rows."
    )
  }

  # Check sample identifier
  if (!"sample_id" %in% colnames(stds_df)) {
    stop(
      "Mass calibration data must contain 'sample_id'."
    )
  }

  # Determine amount column
  amount_col <- paste0(
    element,
    ".Amount"
  )

  required_columns <- c(
    amount_col,
    "area_or_voltage",
    "sample_id"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(stds_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Mass calibration data are missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  # Select the configured calibration standard
  model_df <- stds_df %>%
    dplyr::filter(
      sample_id == calibration_standard
    ) %>%
    dplyr::select(
      dplyr::all_of(amount_col),
      area_or_voltage,
      sample_id
    ) %>%
    dplyr::filter(
      !is.na(.data[[amount_col]]),
      !is.na(area_or_voltage)
    )

  # Check number of valid measurements
  if (nrow(model_df) < 3) {
    stop(
      paste0(
        "Need at least 3 valid measurements of calibration ",
        "standard '",
        calibration_standard,
        "' for ",
        element,
        " mass calibration, got ",
        nrow(model_df),
        "."
      )
    )
  }

  # Check variation in peak area
  if (
    dplyr::n_distinct(
      model_df$area_or_voltage
    ) < 2
  ) {
    stop(
      paste0(
        "Calibration standard '",
        calibration_standard,
        "' does not contain enough variation ",
        "in area_or_voltage."
      )
    )
  }

  # Build the mass calibration model
  formula <- stats::as.formula(
    paste(
      amount_col,
      "~ area_or_voltage"
    )
  )

  model <- stats::lm(
    formula,
    data = model_df
  )

  # Extract model coefficients
  summary_coefs <- summary(model)$coefficients

  slope <- summary_coefs[
    "area_or_voltage",
    "Estimate"
  ]

  intercept <- summary_coefs[
    "(Intercept)",
    "Estimate"
  ]

  # Extract model statistics
  model_summary <- summary(model)

  r_squared <- model_summary$r.squared
  rse <- model_summary$sigma

  # Print model information
  if (verbose) {

    cat(
      "  Calibration standard:",
      calibration_standard,
      "\n"
    )

    cat(
      "  Slope:",
      round(slope, 8),
      "\n"
    )

    cat(
      "  Intercept:",
      round(intercept, 6),
      "\n"
    )

    cat(
      "  R²:",
      round(r_squared, 4),
      "(",
      nrow(model_df),
      "measurements)\n"
    )

    cat(
      "  RSE:",
      round(rse, 4),
      "\n"
    )

    if (r_squared < 0.95) {
      warning(
        paste(
          "R² is low:",
          round(r_squared, 4),
          "< 0.95"
        )
      )
    }

    cat(
      "  ✓ Mass model built\n\n"
    )
  }

  # Create diagnostic plot
  plot <- NULL

  if (create_plot) {

    plot_df <- model_df %>%
      dplyr::mutate(
        predicted = stats::predict(
          model,
          newdata = .
        ),
        residuals =
          .data[[amount_col]] - predicted
      )

    plot <- ggplot(
      plot_df,
      aes(
        x = area_or_voltage,
        y = .data[[amount_col]]
      )
    ) +

      geom_point(
        size = 3,
        alpha = 0.7
      ) +

      geom_smooth(
        method = "lm",
        se = TRUE
      ) +

      labs(
        x = "Peak Area",
        y = paste0(
          element,
          " Amount"
        ),
        title = paste0(
          "Mass Calibration - ",
          element
        ),
        subtitle = paste0(
          "Standard: ",
          calibration_standard,
          "   R² = ",
          round(r_squared, 4)
        )
      ) +

      theme_minimal() +

      theme(
        plot.title = element_text(
          face = "bold",
          size = 12
        ),
        plot.subtitle = element_text(
          size = 10
        )
      )
  }

  # Return model and calibration information
  return(
    list(
      model = model,
      slope = slope,
      intercept = intercept,
      element = element,
      calibration_standard = calibration_standard,
      r_squared = r_squared,
      residual_std_error = rse,
      n_measurements = nrow(model_df),
      plot = plot
    )
  )
}


#' Apply mass calibration to all samples
#'
#' Uses the mass calibration model to predict elemental amounts for
#' all samples and calculates elemental composition.
#'
#' @param df Canonical dataframe.
#' @param mass_model Mass model from build_mass_model().
#' @param config Configuration.
#' @param verbose Print progress.
#'
#' @return Canonical dataframe with {element}.Amount and {element}.Percent.
apply_mass_correction <- function(
    df,
    mass_model,
    config = NULL,
    verbose = TRUE
) {

  element <- mass_model$element

  if (verbose) {
    cat(
      "Applying mass correction for",
      element,
      "\n"
    )

    cat(
      "  Calibration standard:",
      mass_model$calibration_standard,
      "\n"
    )
  }

  amount_col <- paste0(
    element,
    ".Amount"
  )

  percent_col <- paste0(
    element,
    ".Percent"
  )

  required_columns <- c(
    "area_or_voltage",
    "amount"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Data are missing required canonical columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  pred_df <- data.frame(
    area_or_voltage = df$area_or_voltage
  )

  df[[amount_col]] <- stats::predict(
    mass_model$model,
    newdata = pred_df
  )

  df[[percent_col]] <-
    (df[[amount_col]] / df$amount) * 100

  if (verbose) {

    cat(
      "  Applied to",
      nrow(df),
      "samples\n"
    )

    sample_amounts <- df[[amount_col]][
      !is.na(df[[amount_col]])
    ]

    if (length(sample_amounts) > 0) {
      cat(
        "  Range:",
        round(min(sample_amounts, na.rm = TRUE), 4),
        "-",
        round(max(sample_amounts, na.rm = TRUE), 4),
        "\n"
      )
    }

    cat("  ✓ Mass correction applied\n\n")
  }

  return(df)
}
