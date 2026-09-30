# EA adapter: Convert raw EA-IRMS data to the canonical IRMS format

#' Adapt raw EA-IRMS data to the canonical format
#'
#' Uses the EA peripheral schema to map raw instrument columns into the
#' package-wide canonical IRMS data format.
#'
#' @param raw_df Raw EA-IRMS data frame.
#' @param config Experiment configuration loaded from YAML.
#' @param schema EA schema loaded with get_peripheral_schema().
#' @param isotope Isotope system used for processing, e.g. "C13".
#' @param run_id Optional run identifier. Defaults to config$experiment$name.
#' @param instrument Instrument name stored in the canonical data.
#' @param verbose Print progress messages.
#'
#' @return A canonical IRMS data frame with EA-specific columns.
#' @export
adapt_ea_data <- function(
    raw_df,
    config = NULL,
    schema = NULL,
    isotope = "C13",
    run_id = NULL,
    instrument = "EA",
    verbose = TRUE
) {

  # Check raw data
  if (!is.data.frame(raw_df)) {
    stop("'raw_df' must be a data frame.")
  }

  if (nrow(raw_df) == 0) {
    stop("'raw_df' contains no rows.")
  }

  # Require configuration
  if (is.null(config)) {
    stop("An experiment config is required.")
  }

  # Require schema
  if (is.null(schema)) {
    stop("An EA-IRMS schema is required.")
  }

  # Check isotope
  if (
    length(isotope) != 1 ||
    is.na(isotope) ||
    !nzchar(isotope)
  ) {
    stop(
      "'isotope' must contain exactly one non-empty value."
    )
  }

  # Check schema structure
  if (is.null(schema$column_mapping)) {
    stop(
      "The EA-IRMS schema does not contain 'column_mapping'."
    )
  }

  mapping <- schema$column_mapping

  if (!is.list(mapping)) {
    stop(
      "'column_mapping' in the EA-IRMS schema must be a list."
    )
  }

  # Get the element from EA configuration
  element <- get_config_element(config)

  # Resolve isotope-specific mapping
  if (is.null(mapping$isotope_ratio)) {
    stop(
      "The EA-IRMS schema does not contain an ",
      "'isotope_ratio' mapping."
    )
  }

  if (!is.list(mapping$isotope_ratio)) {
    stop(
      "'isotope_ratio' in the EA-IRMS schema must contain ",
      "isotope-specific mappings."
    )
  }

  if (!isotope %in% names(mapping$isotope_ratio)) {
    stop(
      "No isotope ratio mapping found for isotope '",
      isotope,
      "'. Available isotopes: ",
      paste(
        names(mapping$isotope_ratio),
        collapse = ", "
      )
    )
  }

  resolved_mapping <- mapping

  resolved_mapping$isotope_ratio <-
    mapping$isotope_ratio[[isotope]]

  # Required EA mappings needed to create canonical data
  required_mapping <- c(
    "identifier",
    "amount",
    "peak_number",
    "isotope_ratio",
    "area_all"
  )

  missing_mapping <- setdiff(
    required_mapping,
    names(resolved_mapping)
  )

  if (length(missing_mapping) > 0) {
    stop(
      "The EA-IRMS schema is missing required mappings: ",
      paste(
        missing_mapping,
        collapse = ", "
      )
    )
  }

  # Validate that required mappings are single raw column names
  for (field in required_mapping) {

    value <- resolved_mapping[[field]]

    if (
      length(value) != 1 ||
      !is.character(value) ||
      is.na(value) ||
      !nzchar(value)
    ) {
      stop(
        "Mapping '",
        field,
        "' in the EA-IRMS schema must contain ",
        "one raw column name."
      )
    }
  }

  # Optional EA mappings
  optional_mapping <- c(
    "area_44",
    "area_45",
    "area_46",
    "amplitude_44",
    "amplitude_45",
    "amplitude_46",
    "bgd_44",
    "bgd_45",
    "bgd_46",
    "time"
  )

  # Keep only mappings relevant to EA processing
  resolved_mapping <- resolved_mapping[
    names(resolved_mapping) %in%
      c(
        required_mapping,
        optional_mapping
      )
  ]

  # Check that mapped raw columns exist
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
      "The raw EA-IRMS data are missing mapped columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  # Validate EA-specific schema requirements
  validate_ea_raw_data(
    raw_df = raw_df,
    schema = schema,
    mapping = resolved_mapping
  )

  # Determine run ID
  if (is.null(run_id)) {

    if (!is.null(config$experiment_name)) {
      run_id <- config$experiment_name
    } else {
      stop(
        "No 'run_id' was supplied and ",
        "'config$experiment_name' is missing."
      )
    }
  }

  if (
    length(run_id) != 1 ||
    is.na(run_id) ||
    !nzchar(run_id)
  ) {
    stop(
      "'run_id' must contain exactly one non-empty value."
    )
  }

  # Create canonical data
  canonical_df <- data.frame(
    run_id = rep(
      run_id,
      nrow(raw_df)
    ),
    sample_id = raw_df[[resolved_mapping$identifier]],
    element = rep(
      element,
      nrow(raw_df)
    ),
    delta_value = raw_df[[resolved_mapping$isotope_ratio]],
    area_or_voltage = raw_df[[resolved_mapping$area_all]],
    amount = raw_df[[resolved_mapping$amount]],
    instrument = rep(
      instrument,
      nrow(raw_df)
    ),
    is_standard = rep(
      FALSE,
      nrow(raw_df)
    ),
    is_blank = rep(
      FALSE,
      nrow(raw_df)
    ),
    peak_number = raw_df[[resolved_mapping$peak_number]],
    stringsAsFactors = FALSE
  )

  # Add EA-specific columns
  canonical_df <- add_ea_specific_columns(
    canonical_df = canonical_df,
    raw_df = raw_df,
    mapping = resolved_mapping
  )

  # Apply EA-specific filtering rules
  canonical_df <- apply_ea_validation_rules(
    canonical_df = canonical_df,
    schema = schema
  )

  # Validate final EA canonical data
  validate_ea_canonical_data(canonical_df)

  if (verbose) {
    message(
      "Adapted ",
      nrow(canonical_df),
      " EA-IRMS rows ",
      "for element ",
      element,
      " using isotope ",
      isotope,
      "."
    )
  }

  canonical_df
}


#' Get the element specified in the EA experiment configuration
#'
#' @param config Experiment configuration loaded from YAML.
#'
#' @return A single element name.
#' @export
get_config_element <- function(config) {

  if (is.null(config)) {
    stop(
      "An experiment config is required to specify the element."
    )
  }

  if (is.null(config$processing)) {
    stop(
      "The EA experiment config is missing the ",
      "'processing' section."
    )
  }

  if (is.null(config$processing$element)) {
    stop(
      "The EA experiment config is missing ",
      "'processing$element'."
    )
  }

  element <- as.character(
    config$processing$element
  )

  if (
    length(element) != 1 ||
    is.na(element) ||
    !nzchar(element)
  ) {
    stop(
      "'processing$element' must contain exactly one ",
      "non-empty value."
    )
  }

  element
}


#' Validate raw EA data against the EA schema
#'
#' Checks EA-specific required columns and schema validation rules.
#'
#' @param raw_df Raw EA-IRMS data.
#' @param schema EA peripheral schema.
#' @param mapping Resolved EA column mapping.
#'
#' @return Invisibly returns TRUE.
validate_ea_raw_data <- function(
    raw_df,
    schema,
    mapping
) {

  # Check schema validation section
  if (is.null(schema$validation)) {
    return(invisible(TRUE))
  }

  validation <- schema$validation

  # Check schema-defined required columns
  if (!is.null(validation$required_columns)) {

    required_fields <- unlist(
      validation$required_columns,
      use.names = FALSE
    )

    missing_fields <- setdiff(
      required_fields,
      names(mapping)
    )

    if (length(missing_fields) > 0) {
      stop(
        "The EA schema requires mappings for: ",
        paste(
          missing_fields,
          collapse = ", "
        )
      )
    }
  }

  # Check that mapped columns exist in the raw data
  mapped_columns <- unlist(
    mapping,
    use.names = FALSE
  )

  missing_columns <- setdiff(
    mapped_columns,
    names(raw_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "The raw EA data are missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}


#' Apply EA-specific validation and filtering rules
#'
#' Applies rules defined in schema$validation.
#'
#' @param canonical_df Canonical EA data.
#' @param schema EA peripheral schema.
#'
#' @return Filtered canonical EA data.
apply_ea_validation_rules <- function(
    canonical_df,
    schema
) {

  if (is.null(schema$validation)) {
    return(canonical_df)
  }

  validation <- schema$validation

  # Keep only configured peak numbers
  if (!is.null(validation$peak_numbers_to_keep)) {

    peak_numbers <- as.numeric(
      unlist(
        validation$peak_numbers_to_keep,
        use.names = FALSE
      )
    )

    canonical_df <- canonical_df[
      canonical_df$peak_number %in% peak_numbers,
      ,
      drop = FALSE
    ]
  }

  # Exclude configured identifiers
  if (!is.null(validation$exclude_identifiers)) {

    excluded_identifiers <- as.character(
      unlist(
        validation$exclude_identifiers,
        use.names = FALSE
      )
    )

    canonical_df <- canonical_df[
      !canonical_df$sample_id %in% excluded_identifiers,
      ,
      drop = FALSE
    ]
  }

  if (nrow(canonical_df) == 0) {
    stop(
      "No EA data remain after applying schema validation rules."
    )
  }

  canonical_df
}


#' Add EA-specific columns to canonical data
#'
#' @param canonical_df Canonical data frame.
#' @param raw_df Raw EA-IRMS data frame.
#' @param mapping Resolved EA column mapping.
#'
#' @return Canonical data frame with EA-specific columns.
add_ea_specific_columns <- function(
    canonical_df,
    raw_df,
    mapping
) {

  # Area columns
  if (!is.null(mapping$area_44)) {
    canonical_df$area_44 <-
      raw_df[[mapping$area_44]]
  }

  if (!is.null(mapping$area_45)) {
    canonical_df$area_45 <-
      raw_df[[mapping$area_45]]
  }

  if (!is.null(mapping$area_46)) {
    canonical_df$area_46 <-
      raw_df[[mapping$area_46]]
  }

  # Amplitude columns
  if (!is.null(mapping$amplitude_44)) {
    canonical_df$amplitude_44 <-
      raw_df[[mapping$amplitude_44]]
  }

  if (!is.null(mapping$amplitude_45)) {
    canonical_df$amplitude_45 <-
      raw_df[[mapping$amplitude_45]]
  }

  if (!is.null(mapping$amplitude_46)) {
    canonical_df$amplitude_46 <-
      raw_df[[mapping$amplitude_46]]
  }

  # Background columns
  if (!is.null(mapping$bgd_44)) {
    canonical_df$bgd_44 <-
      raw_df[[mapping$bgd_44]]
  }

  if (!is.null(mapping$bgd_45)) {
    canonical_df$bgd_45 <-
      raw_df[[mapping$bgd_45]]
  }

  if (!is.null(mapping$bgd_46)) {
    canonical_df$bgd_46 <-
      raw_df[[mapping$bgd_46]]
  }

  # Time code
  if (!is.null(mapping$time)) {
    canonical_df$time_code <-
      raw_df[[mapping$time]]
  }

  canonical_df
}


#' Validate canonical EA-IRMS data
#'
#' Checks that the EA adapter produced all required canonical and
#' EA-specific fields.
#'
#' @param canonical_df Canonical EA-IRMS data.
#'
#' @return Invisibly returns TRUE.
validate_ea_canonical_data <- function(
    canonical_df
) {

  required_columns <- c(
    "run_id",
    "sample_id",
    "element",
    "delta_value",
    "area_or_voltage",
    "amount",
    "instrument",
    "is_standard",
    "is_blank",
    "peak_number"
  )

  missing_columns <- setdiff(
    required_columns,
    names(canonical_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Canonical EA data are missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  # Check element
  if (
    any(
      is.na(canonical_df$element) |
      canonical_df$element == ""
    )
  ) {
    stop(
      "Canonical EA data contain missing or empty element values."
    )
  }

  # Check sample IDs
  if (any(is.na(canonical_df$sample_id))) {
    warning(
      "Canonical EA data contain missing sample identifiers."
    )
  }

  # Check peak numbers
  if (any(is.na(canonical_df$peak_number))) {
    warning(
      "Canonical EA data contain missing peak numbers."
    )
  }

  # Check required numeric fields
  numeric_columns <- c(
    "delta_value",
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
      "The following canonical EA fields must be numeric: ",
      paste(
        non_numeric,
        collapse = ", "
      )
    )
  }

  # Check logical flags
  if (!is.logical(canonical_df$is_standard)) {
    stop(
      "'is_standard' must be logical."
    )
  }

  if (!is.logical(canonical_df$is_blank)) {
    stop(
      "'is_blank' must be logical."
    )
  }

  invisible(TRUE)
}
