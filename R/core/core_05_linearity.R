# Category 5: Area/voltage-dependent isotope-ratio bias correction
#
# Functions:
#   - build_linearity_model(): Calibrate δ error against area
#   - apply_linearity_correction(): Remove area-dependent δ bias
#
# Dependencies: tidyverse
# Used by: EA, maybe TC/EA
# Requires: One calibration standard measured at multiple amounts

library(tidyverse)

#' Build linearity correction model
#'
#' Creates a linearity model relating measured delta values to peak area.
#'
#' @param stds_df Standards dataframe.
#' @param config Configuration.
#' @param verbose Print progress and model summary (default: TRUE).
#' @param create_plot Create a diagnostic plot (default: TRUE).
#'
#' @return List containing the model, slope, element, R², delta column,
#'   number of standards, and diagnostic plot.
#'
#' @details Uses delta_value_blank_corrected if available, otherwise delta_value.
build_linearity_model <- function(
    stds_df,
    config = NULL,
    verbose = TRUE,
    create_plot = TRUE
) {
  
  if (verbose) cat("Building linearity correction model\n")
  
  # Select the measured delta column.
  if ("delta_value_blank_corrected" %in% colnames(stds_df)) {
    delta_col <- "delta_value_blank_corrected"
    using_col <- "blank-corrected δ"
  } else {
    delta_col <- "delta_value"
    using_col <- "raw δ"
  }
  
  if (verbose) cat("  Using", using_col, "values\n")
  
  # Prepare data for the model using the selected calibration standards.
  model_df <- stds_df %>%
    dplyr::select(
      dplyr::all_of(delta_col),
      area_or_voltage,
      sample_id
    ) %>%
    dplyr::filter(
      !is.na(.data[[delta_col]]),
      !is.na(area_or_voltage)
    )
  
  if (nrow(model_df) < 3) {
    stop(
      paste(
        "Need at least 3 valid calibration standards, got",
        nrow(model_df)
      )
    )
  }
  
  # Fit the linearity model.
  formula <- stats::as.formula(
    paste(delta_col, "~ area_or_voltage")
  )
  
  model <- stats::lm(
    formula,
    data = model_df
  )
  
  # Extract model statistics.
  summary_coefs <- summary(model)$coefficients
  
  slope <- summary_coefs[
    "area_or_voltage",
    "Estimate"
  ]
  
  intercept <- summary_coefs[
    "(Intercept)",
    "Estimate"
  ]
  
  model_summary <- summary(model)
  
  r_squared <- model_summary$r.squared
  rse <- model_summary$sigma
  
  if (verbose) {
    
    cat(
      "  Linearity slope:",
      round(slope, 8),
      "‰ per area unit\n"
    )
    
    cat(
      "  Intercept:",
      round(intercept, 4),
      "‰\n"
    )
    
    cat(
      "  R²:",
      round(r_squared, 4),
      "(",
      nrow(model_df),
      "calibration standards)\n"
    )
    
    cat(
      "  RSE:",
      round(rse, 4),
      "‰\n"
    )
    
    # Interpret the slope.
    if (abs(slope) < 0.0001) {
      
      cat(
        "  Interpretation: Negligible linearity effect\n"
      )
      
    } else {
      
      effect_per_1000 <- slope * 1000
      
      direction <- if (slope > 0) {
        "more positive"
      } else {
        "more negative"
      }
      
      cat(
        "  Interpretation: 1000 area units →",
        round(effect_per_1000, 2),
        "‰",
        direction,
        "\n"
      )
    }
    
    if (r_squared < 0.50) {
      
      warning(
        paste(
          "  ⚠ R² is low:",
          round(r_squared, 4),
          "Linearity correction may not be reliable"
        )
      )
    }
  }
  
  # Create the calibration plot.
  plot <- NULL
  
  if (create_plot) {
    
    # Use exactly the data used to fit the linearity model.
    plot_df <- model_df
    
    plot <- ggplot(
      plot_df,
      aes(
        x = area_or_voltage,
        y = .data[[delta_col]]
      )
    ) +
      
      geom_point(
        size = 3,
        alpha = 0.75
      ) +
      
      geom_smooth(
        method = "lm",
        se = FALSE
      ) +
      
      labs(
        x = "Peak Area",
        y = "Measured δ (‰)",
        title = paste0(
          "Linearity Calibration - ",
          config$element %||% "C"
        )
      ) +
      
      theme_minimal() +
      
      theme(
        plot.title = element_text(
          face = "bold",
          size = 12
        )
      )
  }
  
  if (verbose) {
    cat("  ✓ Linearity model built\n\n")
  }
  
  return(
    list(
      model = model,
      slope = slope,
      intercept = intercept,
      element = config$element %||% "C",
      r_squared = r_squared,
      residual_std_error = rse,
      delta_column = delta_col,
      n_stds = nrow(model_df),
      plot = plot
    )
  )
}

#' Apply linearity correction
#'
#' Removes the area-dependent δ bias predicted by the linearity model.
#'
#' @param df Canonical dataframe.
#' @param linearity_model Linearity model from build_linearity_model().
#' @param verbose Print progress.
#'
#' @return Same dataframe with delta_value_linearity_corrected added.
apply_linearity_correction <- function(
    df,
    linearity_model,
    verbose = TRUE
) {
  
  if (verbose) {
    cat("Applying linearity correction\n")
    cat(
      "  Calibration standard:",
      linearity_model$calibration_standard,
      "\n"
    )
  }
  
  delta_col <- linearity_model$delta_column
  
  if (!delta_col %in% colnames(df)) {
    stop(
      paste(
        "Delta column not found:",
        delta_col,
        "Make sure blank correction was applied first"
      )
    )
  }
  
  required_columns <- c(
    delta_col,
    "area_or_voltage"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Data are missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Predict the area-dependent δ error.
  predicted_error <- stats::predict(
    linearity_model$model,
    newdata = data.frame(
      area_or_voltage = df$area_or_voltage
    )
  )
  
  # Remove the predicted error.
  df$delta_value_linearity_corrected <-
    df[[delta_col]] - predicted_error
  
  if (verbose) {
    
    cat(
      "  Before:",
      round(
        min(
          df[[delta_col]],
          na.rm = TRUE
        ),
        2
      ),
      "to",
      round(
        max(
          df[[delta_col]],
          na.rm = TRUE
        ),
        2
      ),
      "per mil\n"
    )
    
    cat(
      "  After:",
      round(
        min(
          df$delta_value_linearity_corrected,
          na.rm = TRUE
        ),
        2
      ),
      "to",
      round(
        max(
          df$delta_value_linearity_corrected,
          na.rm = TRUE
        ),
        2
      ),
      "per mil\n"
    )
    
    cat("  ✓ Linearity correction applied\n\n")
  }
  
  return(df)
}