# GC-IRMS adapter: convert raw GC data to the canonical format.


#' Adapt raw GC-IRMS data to the canonical format
#'
#' Uses the GC peripheral schema and isotope-system configuration to map
#' raw instrument columns into the package-wide canonical IRMS format.
#'
#' @param raw_df Raw GC-IRMS data frame.
#' @param config Experiment configuration.
#' @param isotope Isotope system used for processing, e.g. "C13".
#' @param verbose Print progress messages.
#'
#' @return A canonical GC-IRMS data frame.
#' @export

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

#' Adapt raw GC-IRMS data to the canonical format
#'
#' Maps raw instrument columns into canonical fields and, when method
#' retention-time data are supplied, assigns compound IDs to uniquely
#' matching peaks.
#'
#' @param raw_df Raw GC-IRMS data frame.
#' @param config Experiment configuration.
#' @param isotope Isotope system used for processing, e.g. "C13".
#' @param verbose Print progress messages.
#' @param gc_methods GC methods list read from YAML, or its file path.
#' @param method_name Method key in the methods YAML. Defaults to
#'   config$processing$gc_method.
#'
#' @return A canonical GC-IRMS data frame.
#' @export
adapt_gc_data <- function(
    raw_df,
    config = NULL,
    isotope = "C13",
    verbose = TRUE,
    gc_methods = NULL,
    method_name = NULL
) {
  if (!is.data.frame(raw_df)) {
    stop("'raw_df' must be a data frame.")
  }

  if (nrow(raw_df) == 0) {
    stop("'raw_df' contains no rows.")
  }

  if (is.null(config)) {
    stop("An experiment config is required.")
  }

  if (
    length(isotope) != 1 ||
    is.na(isotope) ||
    !nzchar(isotope)
  ) {
    stop("'isotope' must contain exactly one non-empty value.")
  }

  isotope_system <- get_isotope_system(
    isotope = isotope,
    config = config
  )

  gc_schemas <- load_peripheral_schemas()

  gc_schema <- get_peripheral_schema(
    schemas = gc_schemas,
    peripheral = "GC"
  )

  if (is.null(gc_schema$column_mapping)) {
    stop("The GC schema does not contain 'column_mapping'.")
  }

  mapping <- gc_schema$column_mapping
  element <- isotope_system$element
  isotope_mass <- isotope_system$isotope_mass
  reference_mass <- isotope_system$reference_mass

  if (is.null(element)) {
    stop("Isotope system '", isotope, "' does not define 'element'.")
  }

  if (is.null(isotope_mass)) {
    stop("Isotope system '", isotope, "' does not define 'isotope_mass'.")
  }

  if (is.null(reference_mass)) {
    stop("Isotope system '", isotope, "' does not define 'reference_mass'.")
  }

  resolve_template <- function(template) {
    if (
      length(template) != 1 ||
      is.na(template) ||
      !nzchar(template)
    ) {
      stop(
        "GC column mapping templates must contain ",
        "one non-empty value."
      )
    }

    template <- stringr::str_replace_all(
      template,
      stringr::fixed("{element}"),
      element
    )

    template <- stringr::str_replace_all(
      template,
      stringr::fixed("{isotope_mass}"),
      as.character(isotope_mass)
    )

    stringr::str_replace_all(
      template,
      stringr::fixed("{reference_mass}"),
      as.character(reference_mass)
    )
  }

  resolved_mapping <- list(
    injection_id = mapping$injection_id,
    sample_id = mapping$sample_id,
    analysis_no = mapping$analysis_no,
    comment = mapping$comment,
    sample_type = mapping$sample_type,
    status = mapping$status,
    start_time = mapping$start_time,
    stop_time = mapping$stop_time,
    peak_number = mapping$peak_number,
    isotope_ratio = resolve_template(mapping$isotope_ratio),
    retention_time = resolve_template(mapping$retention_time),
    peak_amplitude = resolve_template(mapping$peak_amplitude),
    peak_area = resolve_template(mapping$peak_area)
  )

  required_mapping <- c(
    "injection_id",
    "sample_id",
    "analysis_no",
    "comment",
    "sample_type",
    "status",
    "start_time",
    "stop_time",
    "peak_number",
    "isotope_ratio",
    "retention_time",
    "peak_amplitude",
    "peak_area"
  )

  missing_mapping <- required_mapping[
    vapply(
      resolved_mapping[required_mapping],
      is.null,
      logical(1)
    )
  ]

  if (length(missing_mapping) > 0) {
    stop(
      "The GC schema is missing required mappings: ",
      paste(missing_mapping, collapse = ", ")
    )
  }

  missing_columns <- setdiff(
    unlist(resolved_mapping, use.names = FALSE),
    names(raw_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "The raw GC-IRMS data are missing mapped columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  excel_datetime <- function(x) {
    if (inherits(x, "POSIXt")) {
      return(x)
    }

    as.POSIXct(
      as.numeric(x) * 86400,
      origin = "1899-12-30",
      tz = "UTC"
    )
  }

  data <- raw_df |>
    dplyr::transmute(
      injection_id = .data[[resolved_mapping$injection_id]],
      sample_id = .data[[resolved_mapping$sample_id]],
      analysis_no = .data[[resolved_mapping$analysis_no]],
      comment = .data[[resolved_mapping$comment]],
      sample_type = .data[[resolved_mapping$sample_type]],
      status = .data[[resolved_mapping$status]],
      start_time = excel_datetime(
        .data[[resolved_mapping$start_time]]
      ),
      stop_time = excel_datetime(
        .data[[resolved_mapping$stop_time]]
      ),
      peak_number = .data[[resolved_mapping$peak_number]],
      retention_time = as.numeric(
        .data[[resolved_mapping$retention_time]]
      ),
      delta_value = as.numeric(
        .data[[resolved_mapping$isotope_ratio]]
      ),
      peak_amplitude = as.numeric(
        .data[[resolved_mapping$peak_amplitude]]
      ),
      area_or_voltage = as.numeric(
        .data[[resolved_mapping$peak_area]]
      ),
      element = element,
      isotope = isotope,
      amount = NA_real_,
      instrument = "GC-IRMS",
      compound_id = NA_character_,
      compound_match_status = "not_configured"
    )

  if (!is.null(gc_methods)) {
    if (is.character(gc_methods) && length(gc_methods) == 1) {
      if (!file.exists(gc_methods)) {
        stop("GC methods file not found: ", gc_methods)
      }

      gc_methods <- yaml::read_yaml(gc_methods)
    }

    if (is.null(method_name)) {
      method_name <- config$processing$gc_method
    }

    if (!is.null(gc_methods$compounds)) {
      method_config <- gc_methods
    } else {
      method_catalog <- gc_methods$methods %||% gc_methods

      if (
        is.null(method_name) ||
        length(method_name) != 1 ||
        is.na(method_name) ||
        !nzchar(method_name)
      ) {
        stop(
          "Supply 'method_name' or set ",
          "config$processing$gc_method."
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
      " GC-IRMS rows for element ",
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

#' Validate canonical GC-IRMS data
#'
#' Checks that the GC adapter produced the required canonical fields
#' and that their basic data types and identifiers are valid.
#'
#' @param canonical_df Canonical GC-IRMS data frame.
#'
#' @return Invisibly returns TRUE if validation passes.
#' @export
validate_gc_canonical_data <- function(
    canonical_df
) {

  if (!is.data.frame(canonical_df)) {
    stop(
      "'canonical_df' must be a data frame."
    )
  }

  required_columns <- c(
    "injection_id",
    "sample_id",
    "analysis_no",
    "comment",
    "sample_type",
    "status",
    "start_time",
    "stop_time",
    "peak_number",
    "retention_time",
    "delta_value",
    "peak_amplitude",
    "area_or_voltage",
    "element",
    "isotope",
    "amount",
    "instrument"
  )

  missing_columns <- setdiff(
    required_columns,
    names(canonical_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Canonical GC data are missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }


  identifier_columns <- c(
    "injection_id",
    "sample_id",
    "element",
    "isotope",
    "instrument"
  )

  for (column in identifier_columns) {

    values <- canonical_df[[column]]

    if (all(is.na(values))) {
      warning(
        "Canonical GC column '",
        column,
        "' contains only missing values."
      )
    }
  }


  if (
    !is.numeric(canonical_df$injection_id) &&
    !is.character(canonical_df$injection_id)
  ) {
    stop(
      "'injection_id' must be numeric or character."
    )
  }


  if (
    !is.numeric(canonical_df$peak_number) &&
    !is.character(canonical_df$peak_number)
  ) {
    stop(
      "'peak_number' must be numeric or character."
    )
  }


  numeric_columns <- c(
    "retention_time",
    "delta_value",
    "peak_amplitude",
    "area_or_voltage",
    "amount"
  )

  non_numeric <- numeric_columns[
    !vapply(
      canonical_df[numeric_columns],
      is.numeric,
      logical(1)
    )
  ]

  if (length(non_numeric) > 0) {
    stop(
      "The following canonical GC fields must be numeric: ",
      paste(
        non_numeric,
        collapse = ", "
      )
    )
  }


  datetime_columns <- c(
    "start_time",
    "stop_time"
  )

  invalid_datetime <- datetime_columns[
    !vapply(
      canonical_df[datetime_columns],
      function(x) {
        inherits(x, "POSIXct") ||
          inherits(x, "POSIXlt") ||
          inherits(x, "POSIXt")
      },
      logical(1)
    )
  ]

  if (length(invalid_datetime) > 0) {
    stop(
      "The following canonical GC fields must contain ",
      "date-time values: ",
      paste(
        invalid_datetime,
        collapse = ", "
      )
    )
  }


  if (any(
    !is.na(canonical_df$peak_number) &
    suppressWarnings(
      as.numeric(
        as.character(canonical_df$peak_number)
      ) < 1
    ),
    na.rm = TRUE
  )) {
    stop(
      "Canonical GC data contain peak numbers less than 1."
    )
  }


  peak_keys <- paste(
    canonical_df$injection_id,
    canonical_df$peak_number,
    sep = "::"
  )

  duplicate_peaks <- duplicated(
    peak_keys
  ) &
    !is.na(canonical_df$injection_id) &
    !is.na(canonical_df$peak_number)

  if (any(duplicate_peaks)) {
    stop(
      "Canonical GC data contain duplicate ",
      "injection_id + peak_number combinations."
    )
  }


  if (any(
    !is.na(canonical_df$element) &
    canonical_df$element == ""
  )) {
    stop(
      "Canonical GC data contain empty element values."
    )
  }


  if (any(
    !is.na(canonical_df$isotope) &
    canonical_df$isotope == ""
  )) {
    stop(
      "Canonical GC data contain empty isotope values."
    )
  }


  if (any(
    !is.na(canonical_df$instrument) &
    canonical_df$instrument == ""
  )) {
    stop(
      "Canonical GC data contain empty instrument values."
    )
  }


  invisible(TRUE)

}

correct_gc_external_standard <- function(
    canonical_df,
    stds_reference_df,
    standard_id,
    date_column = "start_time",
    compound_column = "compound_id",
    delta_column = "delta_value",
    verbose = TRUE,
    sample_id = standard_id,
    exclude_injection_ids = NULL,
    exclude_compound_ids = NULL
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

  if (
    length(sample_id) != 1 ||
    is.na(sample_id) ||
    !nzchar(as.character(sample_id))
  ) {
    stop("'sample_id' must be one non-empty value.")
  }
  sample_id <- as.character(sample_id)

  if (is.null(exclude_injection_ids)) {
    exclude_injection_ids <- character(0)
  }
  exclude_injection_ids <- unique(trimws(as.character(
    exclude_injection_ids
  )))
  exclude_injection_ids <- exclude_injection_ids[
    !is.na(exclude_injection_ids) & nzchar(exclude_injection_ids)
  ]
  if (is.null(exclude_compound_ids)) {
    exclude_compound_ids <- character(0)
  }
  exclude_compound_ids <- unique(trimws(as.character(
    exclude_compound_ids
  )))
  exclude_compound_ids <- exclude_compound_ids[
    !is.na(exclude_compound_ids) & nzchar(exclude_compound_ids)
  ]

  data <- canonical_df
  data$correction_date <- as.Date(data[[date_column]])

  is_external <- !is.na(data$sample_id) &
    trimws(as.character(data$sample_id)) == sample_id

  if (!any(is_external)) {
    stop(
      "No GC injections found for external standard sample ID '",
      sample_id,
      "'."
    )
  }

  external_injection_ids <- unique(as.character(
    data$injection_id[is_external]
  ))
  excluded_standard_ids_found <- intersect(
    exclude_injection_ids,
    external_injection_ids
  )
  excluded_standard_ids_missing <- setdiff(
    exclude_injection_ids,
    external_injection_ids
  )
  external_compound_ids <- unique(as.character(
    trimws(as.character(data[[compound_column]][is_external]))
  ))
  excluded_standard_compounds_found <- intersect(
    exclude_compound_ids,
    external_compound_ids
  )
  excluded_standard_compounds_missing <- setdiff(
    exclude_compound_ids,
    external_compound_ids
  )
  if (length(excluded_standard_ids_missing) > 0) {
    warning(
      "Configured excluded external-standard injection ID(s) were not found: ",
      paste(excluded_standard_ids_missing, collapse = ", ")
    )
  }
  if (length(excluded_standard_compounds_missing) > 0) {
    warning(
      "Configured excluded external-standard compound ID(s) were not found: ",
      paste(excluded_standard_compounds_missing, collapse = ", ")
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

  standard_rows <- is_external &
    !as.character(data$injection_id) %in% excluded_standard_ids_found &
    !trimws(as.character(data[[compound_column]])) %in%
      excluded_standard_compounds_found
  standard_measurements <- data[
    standard_rows,
    ,
    drop = FALSE
  ] |>
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
      "Applied external-standard correction using reference '",
      standard_id,
      "' for sample ID '",
      sample_id,
      " on ",
      nrow(daily_offsets),
      " date(s); excluded ", length(excluded_standard_ids_found),
      " configured injection(s) from offset calculations."
    )
  }

  list(
    data = data,
    injection_offsets = injection_offsets,
    daily_offsets = daily_offsets,
    excluded_standard_injection_ids = excluded_standard_ids_found,
    missing_excluded_standard_injection_ids = excluded_standard_ids_missing,
    excluded_standard_compound_ids = excluded_standard_compounds_found,
    missing_excluded_standard_compound_ids =
      excluded_standard_compounds_missing
  )
}


