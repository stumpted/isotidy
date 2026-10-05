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
adapt_gc_data <- function(
    raw_df,
    config = NULL,
    isotope = "C13",
    verbose = TRUE
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
    stop(
      "'isotope' must contain exactly one non-empty value."
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

  if (is.null(gc_schema$column_mapping)) {
    stop(
      "The GC schema does not contain 'column_mapping'."
    )
  }

  mapping <- gc_schema$column_mapping

  element <- isotope_system$element
  isotope_mass <- isotope_system$isotope_mass
  reference_mass <- isotope_system$reference_mass

  if (is.null(element)) {
    stop(
      "Isotope system '",
      isotope,
      "' does not define 'element'."
    )
  }

  if (is.null(isotope_mass)) {
    stop(
      "Isotope system '",
      isotope,
      "' does not define 'isotope_mass'."
    )
  }

  if (is.null(reference_mass)) {
    stop(
      "Isotope system '",
      isotope,
      "' does not define 'reference_mass'."
    )
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
      fixed("{element}"),
      element
    )

    template <- stringr::str_replace_all(
      template,
      fixed("{isotope_mass}"),
      as.character(isotope_mass)
    )

    template <- stringr::str_replace_all(
      template,
      fixed("{reference_mass}"),
      as.character(reference_mass)
    )

    template
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

    isotope_ratio = resolve_template(
      mapping$isotope_ratio
    ),

    retention_time = resolve_template(
      mapping$retention_time
    ),

    peak_amplitude = resolve_template(
      mapping$peak_amplitude
    ),

    peak_area = resolve_template(
      mapping$peak_area
    )
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
      paste(
        missing_mapping,
        collapse = ", "
      )
    )
  }


  mapped_columns <- unlist(
    resolved_mapping,
    use.names = FALSE
  )

  missing_columns <- setdiff(
    mapped_columns,
    names(raw_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "The raw GC-IRMS data are missing mapped columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }


  excel_datetime <- function(x) {

    if (inherits(x, "POSIXt")) {
      return(x)
    }

    x <- as.numeric(x)

    as.POSIXct(
      x * 86400,
      origin = "1899-12-30",
      tz = "UTC"
    )
  }


  data <- raw_df |>
    dplyr::transmute(
      injection_id =
        .data[[resolved_mapping$injection_id]],

      sample_id =
        .data[[resolved_mapping$sample_id]],

      analysis_no =
        .data[[resolved_mapping$analysis_no]],

      comment =
        .data[[resolved_mapping$comment]],

      sample_type =
        .data[[resolved_mapping$sample_type]],

      status =
        .data[[resolved_mapping$status]],

      start_time =
        excel_datetime(
          .data[[resolved_mapping$start_time]]
        ),

      stop_time =
        excel_datetime(
          .data[[resolved_mapping$stop_time]]
        ),

      peak_number =
        .data[[resolved_mapping$peak_number]],

      retention_time =
        as.numeric(
          .data[[resolved_mapping$retention_time]]
        ),

      delta_value =
        as.numeric(
          .data[[resolved_mapping$isotope_ratio]]
        ),

      peak_amplitude =
        as.numeric(
          .data[[resolved_mapping$peak_amplitude]]
        ),

      area_or_voltage =
        as.numeric(
          .data[[resolved_mapping$peak_area]]
        ),

      element =
        element,

      isotope =
        isotope,

      amount =
        NA_real_,

      instrument =
        "GC-IRMS"
    )

  validate_gc_canonical_data(
    data
  )


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
