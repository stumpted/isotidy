# GB-IRMS adapter: convert peak-per-row GasBench exports to isotidy's canonical format.

#' Adapt raw GB-IRMS data to the canonical format
#'
#' Maps one row per exported peak into the canonical IRMS format. Labels
#' beginning with a configured standard ID (for example, NaHCO3_1umol) are
#' assigned that base ID for existing standards utilities. The full raw label
#' and the suffix are retained in raw_sample_id and standard_variant.
#'
#' @param raw_df Raw GB-IRMS data frame, such as a table imported from Excel.
#' @param config Experiment config returned by load_IRMS_config().
#' @param schema Optional GB schema. Loaded from package configuration by default.
#' @param isotope Isotope system to adapt, e.g. "C13".
#' @param run_id Optional run identifier. Defaults to the experiment name.
#' @param instrument Instrument label stored in canonical data.
#' @param verbose Print progress information.
#' @return A canonical data frame with one row per raw GB peak.
#' @export
adapt_gb_data <- function(
    raw_df,
    config = NULL,
    schema = NULL,
    isotope = "C13",
    run_id = NULL,
    instrument = "GB-IRMS",
    verbose = TRUE
) {
  if (!is.data.frame(raw_df)) stop("'raw_df' must be a data frame.")
  if (!nrow(raw_df)) stop("'raw_df' contains no rows.")
  if (is.null(config)) stop("An experiment config is required.")

  if (is.null(schema)) {
    schemas <- load_peripheral_schemas()
    schema <- get_peripheral_schema(schemas = schemas, peripheral = "GB")
  }
  mapping <- schema$column_mapping
  if (is.null(mapping)) stop("The GB schema does not contain 'column_mapping'.")

  element <- get_config_element(config)
  isotope_system <- get_isotope_system(isotope = isotope, config = config)
  isotope_mass <- isotope_system$isotope_mass
  reference_isotope_mass <- isotope_system$reference_isotope_mass
  reference_mass <- as.character(isotope_system$reference_mass)

  resolve_template <- function(template) {
    if (length(template) != 1 || is.na(template) || !nzchar(template)) {
      stop("GB column mapping templates must contain one non-empty value.")
    }
    template <- stringr::str_replace_all(as.character(template), stringr::fixed("{element}"), element)
    template <- stringr::str_replace_all(template, stringr::fixed("{isotope_mass}"), as.character(isotope_mass))
    template <- stringr::str_replace_all(template, stringr::fixed("{reference_isotope_mass}"), as.character(reference_isotope_mass))
    stringr::str_replace_all(template, stringr::fixed("{reference_mass}"), reference_mass)
  }

  required_mapping <- c("injection_id", "identifier", "peak_number", "isotope_ratio", "area")
  missing_mapping <- setdiff(required_mapping, names(mapping))
  if (length(missing_mapping)) {
    stop("The GB schema is missing required mappings: ", paste(missing_mapping, collapse = ", "))
  }

  resolve_mass_mapping <- function(field) {
    mass_mapping <- mapping[[field]]
    if (is.null(mass_mapping)) return(NULL)
    if (!is.list(mass_mapping) || !"{mass}" %in% names(mass_mapping)) {
      stop("GB '", field, "' mapping must use '{mass}' as its key.")
    }
    masses <- as.character(unlist(schema$measurement_masses, use.names = FALSE))
    if (!length(masses)) stop("The GB schema must define 'measurement_masses'.")
    resolved <- vapply(masses, function(mass) {
      stringr::str_replace_all(
        resolve_template(mass_mapping[["{mass}"]]),
        stringr::fixed("{mass}"),
        mass
      )
    }, character(1))
    names(resolved) <- masses
    resolved
  }

  resolved <- list(
    injection_id = resolve_template(mapping$injection_id),
    identifier = resolve_template(mapping$identifier),
    amount = if (is.null(mapping$amount)) NULL else resolve_template(mapping$amount),
    amount_unit = if (is.null(mapping$amount_unit)) NULL else resolve_template(mapping$amount_unit),
    peak_number = resolve_template(mapping$peak_number),
    isotope_ratio = resolve_template(mapping$isotope_ratio),
    area = resolve_mass_mapping("area"),
    amplitude = resolve_mass_mapping("amplitude"),
    background = resolve_mass_mapping("background")
  )
  optional_fields <- c(
    sample_dilution = "sample_dilution",
    date = "analysis_date",
    time = "analysis_time",
    retention_time = "retention_time",
    peak_width = "peak_width",
    gas_configuration = "gas_configuration"
  )
  for (field in names(optional_fields)) {
    if (!is.null(mapping[[field]])) resolved[[field]] <- resolve_template(mapping[[field]])
  }

  if (is.null(resolved$area[[reference_mass]])) {
    stop(
      "The GB schema must map area for reference mass ", reference_mass,
      " to populate canonical 'area_or_voltage'."
    )
  }
  raw_columns <- unique(unlist(resolved, use.names = FALSE))
  missing_columns <- setdiff(raw_columns, names(raw_df))
  if (length(missing_columns)) {
    stop("The raw GB data are missing mapped columns: ", paste(missing_columns, collapse = ", "))
  }

  if (is.null(run_id)) {
    run_id <- config$experiment_name
    if (is.null(run_id) && !is.null(config$experiment$name)) run_id <- config$experiment$name
  }
  if (length(run_id) != 1 || is.na(run_id) || !nzchar(run_id)) {
    stop("'run_id' must contain exactly one non-empty value.")
  }

  standard_ids <- as.character(unlist(config$stds_used, use.names = FALSE))
  standard_ids <- standard_ids[!is.na(standard_ids) & nzchar(standard_ids)]
  if (!length(standard_ids)) stop("The experiment config must list standards in 'stds_used'.")

  raw_sample_id <- trimws(as.character(raw_df[[resolved$identifier]]))
  standard_id <- vapply(raw_sample_id, function(label) {
    if (is.na(label) || !nzchar(label)) return(NA_character_)
    matches <- standard_ids[vapply(standard_ids, function(id) {
      identical(label, id) || startsWith(label, paste0(id, "_"))
    }, logical(1))]
    if (!length(matches)) return(NA_character_)
    matches[[which.max(nchar(matches))]]
  }, character(1))
  sample_id <- ifelse(!is.na(standard_id), standard_id, raw_sample_id)
  standard_variant <- ifelse(
    !is.na(standard_id) & nchar(raw_sample_id) > nchar(standard_id),
    substring(raw_sample_id, nchar(standard_id) + 1L),
    NA_character_
  )

  numeric_column <- function(column) suppressWarnings(as.numeric(as.character(column)))
  amount <- if (is.null(resolved$amount)) rep(NA_real_, nrow(raw_df)) else numeric_column(raw_df[[resolved$amount]])
  area_by_mass <- lapply(resolved$area, function(column) numeric_column(raw_df[[column]]))
  names(area_by_mass) <- names(resolved$area)

  canonical_df <- data.frame(
    run_id = rep(as.character(run_id), nrow(raw_df)),
    injection_id = raw_df[[resolved$injection_id]],
    analysis_no = raw_df[[resolved$injection_id]],
    sample_id = sample_id,
    raw_sample_id = raw_sample_id,
    standard_id = standard_id,
    standard_variant = standard_variant,
    element = rep(element, nrow(raw_df)),
    isotope = rep(isotope, nrow(raw_df)),
    delta_value = numeric_column(raw_df[[resolved$isotope_ratio]]),
    area_or_voltage = area_by_mass[[reference_mass]],
    amount = amount,
    amount_unit = if (is.null(resolved$amount_unit)) rep(NA_character_, nrow(raw_df)) else as.character(raw_df[[resolved$amount_unit]]),
    instrument = rep(instrument, nrow(raw_df)),
    is_standard = !is.na(standard_id),
    is_blank = stringr::str_detect(raw_sample_id, stringr::regex(schema$metadata$blank_identifier_pattern %||% "blank", ignore_case = TRUE)),
    peak_number = numeric_column(raw_df[[resolved$peak_number]]),
    stringsAsFactors = FALSE
  )

  mass_fields <- list(area = resolved$area, amplitude = resolved$amplitude, bgd = resolved$background)
  for (prefix in names(mass_fields)) {
    mass_mapping <- mass_fields[[prefix]]
    if (is.null(mass_mapping)) next
    for (mass in names(mass_mapping)) {
      canonical_df[[paste0(prefix, "_", mass)]] <- numeric_column(raw_df[[mass_mapping[[mass]]]])
    }
  }
  for (field in names(optional_fields)) {
    raw_column <- resolved[[field]]
    if (!is.null(raw_column)) canonical_df[[optional_fields[[field]]]] <- raw_df[[raw_column]]
  }

  qc_ids <- as.character(unlist(schema$metadata$quality_control_identifiers, use.names = FALSE))
  canonical_df$is_quality_control <- if (length(qc_ids)) raw_sample_id %in% qc_ids else rep(FALSE, nrow(raw_df))

  if (verbose) {
    message(
      "Adapted ", nrow(canonical_df), " GB-IRMS peak rows for element ",
      element, " using isotope ", isotope, "; matched ",
      sum(canonical_df$is_standard, na.rm = TRUE), " standard peak rows."
    )
  }
  canonical_df
}