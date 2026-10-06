assign_gc_compounds <- function(
    canonical_df,
    method_config,
    verbose = TRUE
) {
  if (is.null(method_config$compounds) ||
      is.null(names(method_config$compounds))) {
    stop("GC method config must contain named 'compounds'.")
  }

  compounds <- method_config$compounds
  compound_ids <- names(compounds)

  get_setting <- function(compound, setting) {
    value <- compound[[setting]]
    if (is.null(value) || length(value) != 1) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(value))
  }

  retention_times <- vapply(
    compounds,
    get_setting,
    numeric(1),
    setting = "retention_time"
  )

  tolerances <- vapply(
    compounds,
    get_setting,
    numeric(1),
    setting = "tolerance"
  )

  valid <- is.finite(retention_times) &
    is.finite(tolerances) &
    tolerances >= 0

  if (any(!valid) && verbose) {
    warning(
      "Skipping GC compounds with missing or invalid retention_time ",
      "or tolerance: ",
      paste(compound_ids[!valid], collapse = ", ")
    )
  }

  if (!any(valid)) {
    stop("GC method config has no compounds with usable retention windows.")
  }

  compound_ids <- compound_ids[valid]
  retention_times <- retention_times[valid]
  tolerances <- tolerances[valid]

  canonical_df$compound_id <- NA_character_
  canonical_df$compound_match_status <- "unmatched"

  missing_rt <- is.na(canonical_df$retention_time) |
    !is.finite(canonical_df$retention_time)
  canonical_df$compound_match_status[missing_rt] <- "missing_retention_time"

  for (i in which(!missing_rt)) {
    matches <- which(
      abs(canonical_df$retention_time[[i]] - retention_times) <= tolerances
    )

    if (length(matches) == 1) {
      canonical_df$compound_id[[i]] <- compound_ids[[matches]]
      canonical_df$compound_match_status[[i]] <- "matched"
    } else if (length(matches) > 1) {
      canonical_df$compound_match_status[[i]] <- "ambiguous"
    }
  }

  if (verbose) {
    n_matched <- sum(canonical_df$compound_match_status == "matched")
    n_ambiguous <- sum(canonical_df$compound_match_status == "ambiguous")
    n_unmatched <- sum(
      canonical_df$compound_match_status %in%
        c("unmatched", "missing_retention_time")
    )

    message(
      "GC compound assignment: ",
      n_matched, " matched, ",
      n_ambiguous, " ambiguous, ",
      n_unmatched, " unmatched."
    )
  }

  canonical_df
}

#' @param gc_methods GC methods list read from the methods YAML, or its file path.
#' @param method_name Method key in the methods YAML. Defaults to
#'   config$processing$gc_method.

adapt_gc_data <- function(
    raw_df,
    config = NULL,
    isotope = "C13",
    verbose = TRUE,
    gc_methods = NULL,
    method_name = NULL
) {

  excel_datetime <- function(x) {
    x <- as.numeric(x)

    as.POSIXct(
      x * 86400,
      origin = "1899-12-30",
      tz = "UTC"
    )
  }

  isotope_system <- get_isotope_system(
    isotope = isotope,
    config = config
  )

  gc_schema <- get_peripheral_schema(
    peripheral = "GC",
    config = config
  )

  element <- isotope_system$element
  isotope_mass <- isotope_system$isotope_mass
  reference_mass <- isotope_system$reference_mass

  resolve_column <- function(template) {
    template |>
      stringr::str_replace_all(
        fixed("{element}"),
        element
      ) |>
      stringr::str_replace_all(
        fixed("{isotope_mass}"),
        as.character(isotope_mass)
      ) |>
      stringr::str_replace_all(
        fixed("{reference_mass}"),
        as.character(reference_mass)
      )
  }

  column_mapping <- gc_schema$column_mapping

  isotope_ratio_column <- resolve_column(
    column_mapping$isotope_ratio
  )

  retention_time_column <- resolve_column(
    column_mapping$retention_time
  )

  peak_amplitude_column <- resolve_column(
    column_mapping$peak_amplitude
  )

  peak_area_column <- resolve_column(
    column_mapping$peak_area
  )

  required_columns <- c(
    column_mapping$injection_id,
    column_mapping$sample_id,
    column_mapping$analysis_no,
    column_mapping$comment,
    column_mapping$sample_type,
    column_mapping$status,
    column_mapping$start_time,
    column_mapping$stop_time,
    column_mapping$peak_number,
    isotope_ratio_column,
    retention_time_column,
    peak_amplitude_column,
    peak_area_column
  )

  missing_columns <- setdiff(
    required_columns,
    names(raw_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing required GC columns for isotope ",
      isotope,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data <- raw_df |>
    dplyr::transmute(
      injection_id = .data[[column_mapping$injection_id]],
      sample_id = .data[[column_mapping$sample_id]],
      analysis_no = .data[[column_mapping$analysis_no]],
      comment = .data[[column_mapping$comment]],
      sample_type = .data[[column_mapping$sample_type]],
      status = .data[[column_mapping$status]],
      start_time = excel_datetime(
        .data[[column_mapping$start_time]]
      ),
      stop_time = excel_datetime(
        .data[[column_mapping$stop_time]]
      ),
      peak_number = .data[[column_mapping$peak_number]],
      retention_time = as.numeric(
        .data[[retention_time_column]]
      ),
      delta_value = as.numeric(
        .data[[isotope_ratio_column]]
      ),
      peak_amplitude = as.numeric(
        .data[[peak_amplitude_column]]
      ),
      area_or_voltage = as.numeric(
        .data[[peak_area_column]]
      ),
      element = element,
      isotope = isotope,
      amount = NA_real_,
      instrument = "GC-IRMS",
      compound_id = NA_character_,
      compound_match_status = "not_configured"
    )

  if (!is.null(gc_methods)) {
    if (is.character(gc_methods) &&
        length(gc_methods) == 1 &&
        file.exists(gc_methods)) {
      gc_methods <- yaml::read_yaml(gc_methods)
    }

    if (is.null(method_name)) {
      method_name <- config$processing$gc_method
    }

    # Accept either a direct method definition or a catalog keyed by method name.
    if (!is.null(gc_methods$compounds)) {
      method_config <- gc_methods
    } else {
      method_catalog <- gc_methods$methods %||% gc_methods

      if (is.null(method_name) ||
          length(method_name) != 1 ||
          is.na(method_name) ||
          !nzchar(method_name)) {
        stop(
          "Supply 'method_name' or set config$processing$gc_method."
        )
      }

      method_config <- method_catalog[[method_name]]

      if (is.null(method_config)) {
        stop("GC method not found in methods config: ", method_name)
      }
    }

    data <- assign_gc_compounds(
      canonical_df = data,
      method_config = method_config,
      verbose = verbose
    )
  }

  if (verbose) {
    cat(
      "Adapted ",
      nrow(data),
      " GC rows for element ",
      element,
      " using isotope ",
      isotope,
      " and reference mass ",
      reference_mass,
      "\n",
      sep = ""
    )
  }

  data
}

correct_gc_external_standard <- function(
    canonical_df,
    stds_reference_df,
    standard_id,
    date_column = "start_time",
    compound_column = "compound_id",
    delta_column = "delta_value",
    verbose = TRUE
) {
  required_data <- c(
    "injection_id",
    "sample_id",
    date_column,
    compound_column,
    delta_column
  )

  missing_data <- setdiff(required_data, names(canonical_df))
  if (length(missing_data) > 0) {
    stop(
      "GC data are missing required columns: ",
      paste(missing_data, collapse = ", ")
    )
  }

  required_reference <- c(
    "standard_id",
    "compound_id",
    "delta_value"
  )

  missing_reference <- setdiff(required_reference, names(stds_reference_df))
  if (length(missing_reference) > 0) {
    stop(
      "Standards reference data are missing required columns: ",
      paste(missing_reference, collapse = ", ")
    )
  }

  if (
    length(standard_id) != 1 ||
    is.na(standard_id) ||
    !nzchar(standard_id)
  ) {
    stop("'standard_id' must be one non-empty value.")
  }

  data <- canonical_df
  data$correction_date <- as.Date(data[[date_column]])

  is_external <- !is.na(data$sample_id) &
    trimws(as.character(data$sample_id)) == standard_id

  if (!any(is_external)) {
    stop(
      "No GC injections found for external standard '",
      standard_id,
      "'."
    )
  }

  references <- stds_reference_df |>
    dplyr::filter(
      as.character(.data$standard_id) == standard_id,
      !is.na(.data$compound_id)
    ) |>
    dplyr::transmute(
      compound_id = as.character(.data$compound_id),
      reference_delta = as.numeric(.data$delta_value)
    )

  if (nrow(references) == 0) {
    stop(
      "No compound-specific reference values found for standard '",
      standard_id,
      "'."
    )
  }

  if (anyDuplicated(references$compound_id)) {
    stop(
      "Reference values contain duplicate compound IDs for standard '",
      standard_id,
      "'."
    )
  }

  standard_measurements <- data[is_external, , drop = FALSE] |>
    dplyr::transmute(
      injection_id = .data$injection_id,
      correction_date = .data$correction_date,
      compound_id = as.character(.data[[compound_column]]),
      measured_delta = as.numeric(.data[[delta_column]])
    ) |>
    dplyr::left_join(
      references,
      by = "compound_id",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      offset = .data$reference_delta - .data$measured_delta
    )

  unmatched <- sum(
    !is.na(standard_measurements$compound_id) &
      is.na(standard_measurements$reference_delta)
  )

  if (unmatched > 0) {
    warning(
      unmatched,
      " external-standard peak measurements had no matching ",
      "compound reference value."
    )
  }

  injection_offsets <- standard_measurements |>
    dplyr::filter(
      !is.na(.data$correction_date),
      !is.na(.data$injection_id),
      is.finite(.data$offset)
    ) |>
    dplyr::group_by(
      .data$injection_id,
      .data$correction_date
    ) |>
    dplyr::summarise(
      delta_offset = mean(.data$offset),
      n_compounds = dplyr::n(),
      .groups = "drop"
    )

  daily_offsets <- injection_offsets |>
    dplyr::group_by(.data$correction_date) |>
    dplyr::summarise(
      delta_offset = mean(.data$delta_offset),
      n_injections = dplyr::n(),
      .groups = "drop"
    )

  if (nrow(daily_offsets) == 0) {
    stop(
      "No usable compound matches were found for external standard '",
      standard_id,
      "'."
    )
  }

  data$delta_value_external_corrected <- as.numeric(
    data[[delta_column]]
  )

  # Leave the external-standard measurements unchanged. Correct all other
  # injections using the mean offset for their calendar date.
  daily_index <- match(
    data$correction_date,
    daily_offsets$correction_date
  )

  has_daily_offset <- !is_external & !is.na(daily_index)

  data$delta_value_external_corrected[has_daily_offset] <-
    as.numeric(data[[delta_column]][has_daily_offset]) +
    daily_offsets$delta_offset[daily_index[has_daily_offset]]

  # Match the Python workflow: non-standard injections without a same-day
  # external standard receive missing corrected values.
  no_daily_offset <- !is_external & is.na(daily_index)
  data$delta_value_external_corrected[no_daily_offset] <- NA_real_

  if (verbose) {
    message(
      "Applied external-standard correction from ",
      standard_id,
      " on ",
      nrow(daily_offsets),
      " date(s)."
    )
  }

  list(
    data = data,
    injection_offsets = injection_offsets,
    daily_offsets = daily_offsets
  )
}
