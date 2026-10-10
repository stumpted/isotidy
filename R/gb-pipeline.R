# GB-IRMS processing: adapt, filter, and correct run drift by analysis number.

#' Process GB-IRMS data with analysis-number drift correction
#'
#' Adapts a GB export and applies name exclusions, then selects sample peaks
#' using the inclusive configured low/high range. Other peaks remain available
#' in all_peak_data and excluded_peaks for inspection.
#'
#' @param raw_df Raw GB-IRMS data frame.
#' @param config Experiment config returned by load_IRMS_config().
#' @param schema Optional GB peripheral schema.
#' @param isotope Isotope system to process, e.g. "C13".
#' @param run_id Optional run identifier. Defaults to the experiment name.
#' @param exclude_patterns Name patterns to exclude. If NULL, uses
#'   config$exclude_patterns or the GB defaults.
#' @param exclude_analysis_numbers Analysis numbers to exclude. If NULL, uses
#'   config$processing$exclude_analysis_numbers or excludes none.
#' @param drift_standard_id Standard used to estimate run drift. Defaults to
#'   config$processing$drift_standard_id or "ETH4".
#' @param apply_drift_correction Apply a linear analysis-number correction
#'   when at least three drift-standard analyses are available.
#' @param sample_peak_numbers Optional explicit peak-number vector overriding
#'   the configured low/high range.
#' @param verbose Print progress information.
#' @return Selected sample peak data and output, all name-filtered peak rows,
#'   excluded rows and peaks, and the settings used.
#' @export
process_gb <- function(
    raw_df,
    config,
    schema = NULL,
    isotope = "C13",
    run_id = NULL,
    exclude_patterns = NULL,
    exclude_analysis_numbers = NULL,
    drift_standard_id = NULL,
    apply_drift_correction = TRUE,
    sample_peak_numbers = NULL,
    verbose = TRUE
) {
  if (is.null(config)) stop("An experiment config is required.")

  if (is.null(drift_standard_id)) {
    drift_standard_id <- config$processing$drift_standard_id %||% "ETH4"
  }
  if (length(drift_standard_id) != 1L || is.na(drift_standard_id) ||
      !nzchar(as.character(drift_standard_id))) {
    stop("'drift_standard_id' must be one non-empty standard ID.")
  }
  drift_standard_id <- as.character(drift_standard_id)
  if (!is.logical(apply_drift_correction) || length(apply_drift_correction) != 1L ||
      is.na(apply_drift_correction)) {
    stop("'apply_drift_correction' must be TRUE or FALSE.")
  }

  if (is.null(exclude_analysis_numbers)) {
    exclude_analysis_numbers <- config$processing$exclude_analysis_numbers %||%
      character()
  }
  exclude_analysis_numbers <- as.character(unlist(
    exclude_analysis_numbers,
    use.names = FALSE
  ))
  if (anyNA(exclude_analysis_numbers) || any(!nzchar(exclude_analysis_numbers))) {
    stop("'exclude_analysis_numbers' cannot contain missing or empty values.")
  }

  if (is.null(sample_peak_numbers)) {
    low <- config$processing$sample_peak_number_low
    high <- config$processing$sample_peak_number_high
    if (is.null(low) || is.null(high)) {
      stop(
        "Configure processing.sample_peak_number_low and ",
        "processing.sample_peak_number_high, or supply sample_peak_numbers."
      )
    }
    low <- suppressWarnings(as.numeric(low))
    high <- suppressWarnings(as.numeric(high))
    if (length(low) != 1L || length(high) != 1L ||
        !is.finite(low) || !is.finite(high) ||
        low != floor(low) || high != floor(high) || low < 1 || high < low) {
      stop("GB sample peak range must be two positive, ordered integers.")
    }
    sample_peak_numbers <- seq.int(low, high)
  } else {
    sample_peak_numbers <- suppressWarnings(
      as.numeric(unlist(sample_peak_numbers, use.names = FALSE))
    )
    if (!length(sample_peak_numbers) || any(!is.finite(sample_peak_numbers)) ||
        any(sample_peak_numbers < 1) ||
        any(sample_peak_numbers != floor(sample_peak_numbers))) {
      stop("'sample_peak_numbers' must contain positive integers.")
    }
    sample_peak_numbers <- unique(sample_peak_numbers)
  }

  if (is.null(exclude_patterns)) {
    exclude_patterns <- config$exclude_patterns %||%
      c("on_off", "linearity", "blank", "CIT", "f2_si_26087")
  }
  exclude_patterns <- as.character(unlist(exclude_patterns, use.names = FALSE))
  if (anyNA(exclude_patterns) || any(!nzchar(exclude_patterns))) {
    stop("'exclude_patterns' cannot contain missing or empty values.")
  }

  data <- adapt_gb_data(
    raw_df = raw_df,
    config = config,
    schema = schema,
    isotope = isotope,
    run_id = run_id,
    verbose = verbose
  )

  exclude_analysis_rows <- as.character(data$analysis_no) %in%
    exclude_analysis_numbers
  excluded_analyses <- data[exclude_analysis_rows, , drop = FALSE]
  excluded_analyses$analysis_exclusion_reason <- "analysis number listed for exclusion"
  candidate_data <- data[!exclude_analysis_rows, , drop = FALSE]

  matched_patterns <- lapply(exclude_patterns, function(pattern) {
    stringr::str_detect(candidate_data$sample_id, pattern) %in% TRUE
  })
  excluded_rows <- if (length(matched_patterns)) {
    Reduce(function(left, right) left | right, matched_patterns)
  } else {
    rep(FALSE, nrow(candidate_data))
  }

  excluded_data <- candidate_data[excluded_rows, , drop = FALSE]
  if (nrow(excluded_data)) {
    excluded_data$exclusion_reason <- vapply(
      which(excluded_rows),
      function(i) {
        hits <- exclude_patterns[vapply(
          exclude_patterns,
          function(pattern) isTRUE(stringr::str_detect(candidate_data$sample_id[[i]], pattern)),
          logical(1)
        )]
        paste(hits, collapse = "; ")
      },
      character(1)
    )
  } else {
    excluded_data$exclusion_reason <- character()
  }

  # clean_data() is shared with EA/GC. Disable its EA-only peak selector so
  # that every GB row for a retained injection remains available downstream.
  cleaning_config <- config
  cleaning_config$exclude_patterns <- exclude_patterns
  cleaning_config$peak_number <- NULL
  cleaned_data <- clean_data(
    df = candidate_data,
    config = cleaning_config,
    verbose = verbose
  )

  output_columns <- c(
    "run_id",
    "injection_id",
    "analysis_no",
    "sample_id",
    "raw_sample_id",
    "standard_id",
    "standard_variant",
    "element",
    "isotope",
    "peak_number",
    "delta_value",
    "delta_value_raw",
    "drift_offset",
    "drift_correction_applied",
    "drift_extrapolated",
    "area_or_voltage",
    "amount",
    "amount_unit",
    "is_standard",
    "is_blank",
    "is_quality_control"
  )
  all_peak_data <- cleaned_data
  in_sample_range <- !is.na(all_peak_data$peak_number) &
    all_peak_data$peak_number %in% sample_peak_numbers
  excluded_peaks <- all_peak_data[!in_sample_range, , drop = FALSE]
  excluded_peaks$peak_exclusion_reason <- "outside configured sample peak numbers"
  sample_data <- all_peak_data[in_sample_range, , drop = FALSE]

  # Fit against one mean per standard analysis so repeated peaks from a single
  # vial are not treated as independent drift observations.
  sample_data$delta_value_raw <- sample_data$delta_value
  sample_data$drift_offset <- 0
  sample_data$drift_correction_applied <- FALSE
  sample_data$drift_extrapolated <- FALSE
  drift_analysis <- data.frame(
    analysis_no = numeric(),
    mean_delta = numeric(),
    n_peaks = integer()
  )
  drift_correction <- list(
    applied = FALSE,
    reason = if (apply_drift_correction) "not enough reference analyses" else "disabled",
    standard_id = drift_standard_id,
    n_standard_analyses = 0L,
    extrapolated_rows = 0L,
    analysis_range = c(NA_real_, NA_real_),
    reference_delta = NA_real_,
    slope_per_analysis = NA_real_,
    intercept = NA_real_,
    r_squared = NA_real_,
    slope_p_value = NA_real_,
    residual_sd = NA_real_,
    standard_analyses = drift_analysis
  )

  if (apply_drift_correction) {
    drift_rows <- sample_data[
      !is.na(sample_data$standard_id) &
        sample_data$standard_id == drift_standard_id &
        is.finite(sample_data$delta_value),
      ,
      drop = FALSE
    ]
    drift_rows$analysis_no_numeric <- suppressWarnings(
      as.numeric(as.character(drift_rows$analysis_no))
    )
    drift_rows <- drift_rows[is.finite(drift_rows$analysis_no_numeric), , drop = FALSE]
    if (nrow(drift_rows)) {
      drift_analysis <- dplyr::summarise(
        dplyr::group_by(drift_rows, analysis_no_numeric),
        mean_delta = mean(delta_value),
        n_peaks = dplyr::n(),
        .groups = "drop"
      )
      names(drift_analysis)[names(drift_analysis) == "analysis_no_numeric"] <- "analysis_no"
      drift_analysis <- drift_analysis[order(drift_analysis$analysis_no), , drop = FALSE]
    }

    drift_correction$n_standard_analyses <- nrow(drift_analysis)
    drift_correction$analysis_range <- if (nrow(drift_analysis)) {
      range(drift_analysis$analysis_no)
    } else {
      c(NA_real_, NA_real_)
    }
    drift_correction$standard_analyses <- drift_analysis

    if (nrow(drift_analysis) >= 3L && length(unique(drift_analysis$analysis_no)) >= 3L) {
      drift_fit <- stats::lm(mean_delta ~ analysis_no, data = drift_analysis)
      reference_delta <- mean(stats::fitted(drift_fit))
      drift_analysis$fitted_delta <- stats::fitted(drift_fit)
      drift_analysis$drift_offset <- drift_analysis$fitted_delta - reference_delta
      drift_analysis$corrected_delta <-
        drift_analysis$mean_delta - drift_analysis$drift_offset
      drift_analysis$residual <- stats::residuals(drift_fit)
      drift_correction$standard_analyses <- drift_analysis
      sample_analysis_no <- suppressWarnings(
        as.numeric(as.character(sample_data$analysis_no))
      )
      valid_analysis <- is.finite(sample_analysis_no)
      fitted_delta <- rep(NA_real_, nrow(sample_data))
      fitted_delta[valid_analysis] <- stats::predict(
        drift_fit,
        newdata = data.frame(analysis_no = sample_analysis_no[valid_analysis])
      )
      drift_offset <- fitted_delta - reference_delta
      correctable <- is.finite(drift_offset) & is.finite(sample_data$delta_value)
      sample_data$drift_offset[correctable] <- drift_offset[correctable]
      sample_data$delta_value[correctable] <-
        sample_data$delta_value_raw[correctable] - drift_offset[correctable]
      sample_data$drift_correction_applied <- correctable
      sample_data$drift_extrapolated[correctable] <-
        sample_analysis_no[correctable] < min(drift_analysis$analysis_no) |
        sample_analysis_no[correctable] > max(drift_analysis$analysis_no)

      drift_correction$applied <- any(correctable)
      drift_correction$reason <- if (any(correctable)) "applied" else "no numeric sample analysis numbers"
      drift_correction$reference_delta <- reference_delta
      drift_correction$slope_per_analysis <- unname(stats::coef(drift_fit)[["analysis_no"]])
      drift_correction$intercept <- unname(stats::coef(drift_fit)[["(Intercept)"]])
      drift_fit_summary <- summary(drift_fit)
      drift_correction$r_squared <- drift_fit_summary$r.squared
      drift_correction$slope_p_value <- drift_fit_summary$coefficients[
        "analysis_no", "Pr(>|t|)"
      ]
      drift_correction$residual_sd <- drift_fit_summary$sigma
      drift_correction$extrapolated_rows <- sum(sample_data$drift_extrapolated)
    }
  }

  output <- prepare_output(
    df = sample_data,
    output_cols = intersect(output_columns, names(sample_data)),
    verbose = verbose
  )

  if (verbose) {
    message(
      "GB exclusions retained ", nrow(all_peak_data),
      " peak rows; sample peak selection retained ", nrow(sample_data),
      " rows (peak numbers ", paste(sample_peak_numbers, collapse = ", "),
      ") and set aside ", nrow(excluded_peaks), " other peaks."
    )
    if (apply_drift_correction) {
      if (isTRUE(drift_correction$applied)) {
        message(
          "Applied linear drift correction from ", drift_standard_id,
          " using ", drift_correction$n_standard_analyses,
          " analyses; slope = ", signif(drift_correction$slope_per_analysis, 4),
          " delta units per analysis; ", drift_correction$extrapolated_rows,
          " sample peak rows are outside the standard-analysis range."
        )
      } else {
        message("GB drift correction not applied: ", drift_correction$reason, ".")
      }
    }
  }

  list(
    data = sample_data,
    output = output,
    all_peak_data = all_peak_data,
    excluded_data = excluded_data,
    excluded_analyses = excluded_analyses,
    excluded_peaks = excluded_peaks,
    exclude_patterns = exclude_patterns,
    exclude_analysis_numbers = exclude_analysis_numbers,
    sample_peak_numbers = sample_peak_numbers,
    drift_correction = drift_correction
  )
}
