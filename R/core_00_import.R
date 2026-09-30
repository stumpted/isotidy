# Category 0: Configuration and canonical data handling
#
# Functions:
#   - load_peripheral_schemas(): Load peripheral schemas
#   - get_peripheral_schema(): Retrieve one peripheral schema
#   - validate_peripheral_schemas(): Validate all peripheral schemas
#   - validate_peripheral_schema(): Validate one peripheral schema
#   - load_lab_standards(): Load laboratory standards
#   - load_IRMS_config(): Load experimental configuration
#   - validate_config(): Validate experimental configuration
#   - create_config_template(): Create a configuration template
#   - normalize_data(): Convert raw data to canonical format
#
# Dependencies: yaml, readr, purrr, tibble
# Used by: All peripheral-specific wrappers and processing pipelines


#' Load peripheral schemas
#'
#' Loads the peripheral schema YAML file from the installed package.
#'
#' @param path Optional path to a schema YAML file.
#' @param config_dir Optional configuration directory. Used when developing
#'   the package locally.
#' @return A named list containing the peripheral schemas.
#' @export
load_peripheral_schemas <- function(
    path = NULL,
    config_dir = NULL
) {

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      "The 'yaml' package is required to load peripheral schemas. ",
      "Install it with install.packages('yaml')."
    )
  }

  if (!is.null(path)) {

    schema_path <- path

  } else if (!is.null(config_dir)) {

    schema_path <- file.path(
      config_dir,
      "peripheral_schemas.yml"
    )

  } else {

    schema_path <- system.file(
      "config",
      "peripheral_schemas.yml",
      package = "isotidy"
    )

    if (schema_path == "") {
      schema_path <- file.path(
        "inst",
        "config",
        "peripheral_schemas.yml"
      )
    }
  }

  if (!file.exists(schema_path)) {
    stop(
      "Peripheral schema file not found: ",
      schema_path
    )
  }

  schemas <- yaml::read_yaml(schema_path)

  if (is.null(schemas) || length(schemas) == 0) {
    stop(
      "Peripheral schema file is empty: ",
      schema_path
    )
  }

  validate_peripheral_schemas(schemas)

  schemas
}


#' Get one peripheral schema
#'
#' Retrieves a schema for a specific peripheral.
#'
#' @param schemas Loaded schemas from load_peripheral_schemas().
#' @param peripheral Name of the peripheral, such as "EA".
#' @return The requested peripheral schema.
#' @export
get_peripheral_schema <- function(
    schemas,
    peripheral
) {

  if (is.null(schemas) || !is.list(schemas)) {
    stop(
      "schemas must be a list returned by load_peripheral_schemas()."
    )
  }

  if (
    length(peripheral) != 1 ||
    is.na(peripheral) ||
    !nzchar(peripheral)
  ) {
    stop(
      "peripheral must be a single non-empty character value."
    )
  }

  available <- if (!is.null(schemas$peripherals)) {
    schemas$peripherals
  } else {
    schemas
  }

  if (!peripheral %in% names(available)) {
    stop(
      "No schema found for peripheral '",
      peripheral,
      "'. Available schemas: ",
      paste(
        names(available),
        collapse = ", "
      )
    )
  }

  schema <- available[[peripheral]]

  validate_peripheral_schema(
    schema,
    peripheral
  )

  schema
}


#' Validate all peripheral schemas
#'
#' Performs generic structural validation of all loaded peripheral schemas.
#'
#' @param schemas Loaded schemas.
#' @return Invisibly returns TRUE.
#' @export
validate_peripheral_schemas <- function(
    schemas
) {

  available <- if (!is.null(schemas$peripherals)) {
    schemas$peripherals
  } else {
    schemas
  }

  if (
    !is.list(available) ||
    length(available) == 0
  ) {
    stop("No peripheral schemas were found.")
  }

  if (is.null(names(available))) {
    stop(
      "Peripheral schemas must have named entries."
    )
  }

  for (peripheral in names(available)) {

    validate_peripheral_schema(
      available[[peripheral]],
      peripheral
    )
  }

  invisible(TRUE)
}


#' Validate one peripheral schema
#'
#' Performs generic structural validation without imposing
#' peripheral-specific fields or processing rules.
#'
#' @param schema A single peripheral schema.
#' @param peripheral Name of the peripheral.
#' @return Invisibly returns TRUE.
#' @export
validate_peripheral_schema <- function(
    schema,
    peripheral = "unknown"
) {

  if (!is.list(schema)) {
    stop(
      "Schema for peripheral '",
      peripheral,
      "' must be a YAML mapping/list."
    )
  }

  if (
    is.null(schema$column_mapping)
  ) {
    stop(
      "Schema for peripheral '",
      peripheral,
      "' is missing 'column_mapping'."
    )
  }

  if (!is.list(schema$column_mapping)) {
    stop(
      "'column_mapping' for peripheral '",
      peripheral,
      "' must be a YAML mapping/list."
    )
  }

  if (length(schema$column_mapping) == 0) {
    stop(
      "'column_mapping' for peripheral '",
      peripheral,
      "' is empty."
    )
  }

  if (is.null(names(schema$column_mapping))) {
    stop(
      "'column_mapping' for peripheral '",
      peripheral,
      "' must contain named mappings."
    )
  }

  # Recursively validate mapping values.
  validate_mapping <- function(
    mapping,
    field_path = "column_mapping"
  ) {

    if (!is.list(mapping)) {

      if (
        !is.character(mapping) ||
        length(mapping) != 1 ||
        is.na(mapping) ||
        !nzchar(mapping)
      ) {
        stop(
          "Mapping '",
          field_path,
          "' in peripheral '",
          peripheral,
          "' must contain one non-empty raw column name."
        )
      }

      return(invisible(TRUE))
    }

    if (length(mapping) == 0) {
      stop(
        "Mapping '",
        field_path,
        "' in peripheral '",
        peripheral,
        "' is empty."
      )
    }

    if (is.null(names(mapping))) {
      stop(
        "Mapping '",
        field_path,
        "' in peripheral '",
        peripheral,
        "' must contain named entries."
      )
    }

    for (field in names(mapping)) {

      if (
        is.null(field) ||
        !nzchar(field)
      ) {
        stop(
          "A mapping in peripheral '",
          peripheral,
          "' has an empty field name."
        )
      }

      validate_mapping(
        mapping[[field]],
        paste0(
          field_path,
          "$",
          field
        )
      )
    }

    invisible(TRUE)
  }

  validate_mapping(
    schema$column_mapping
  )

  invisible(TRUE)
}


#' Load laboratory standards
#'
#' Loads and validates the laboratory standards reference database.
#'
#' @param path Path to the standards YAML or CSV file.
#' @param peripheral Optional peripheral name used to filter applicable standards.
#' @param standard_ids Optional standard IDs to load.
#' @param verbose Print loading information.
#' @return A data frame containing laboratory standards.
#' @export
load_lab_standards <- function(
    path,
    peripheral = NULL,
    standard_ids = NULL,
    verbose = TRUE
) {

  if (!file.exists(path)) {
    stop(
      "Laboratory standards file not found: ",
      path
    )
  }

  if (grepl(
    "\\.ya?ml$",
    path,
    ignore.case = TRUE
  )) {

    standards <- yaml::read_yaml(path)

    if (
      is.null(standards) ||
      length(standards) == 0
    ) {
      stop(
        "Laboratory standards file is empty: ",
        path
      )
    }

    standards <- purrr::map_dfr(
      standards,
      tibble::as_tibble
    )

  } else if (grepl(
    "\\.csv$",
    path,
    ignore.case = TRUE
  )) {

    standards <- readr::read_csv(
      path,
      show_col_types = FALSE
    )

  } else {

    stop(
      "Unsupported laboratory standards file format. ",
      "Use YAML or CSV."
    )
  }

  required_columns <- c(
    "element",
    "delta_value",
    "delta_uncertainty",
    "scale",
    "element_fraction",
    "standard_id",
    "applicable_peripherals"
  )

  missing_columns <- setdiff(
    required_columns,
    names(standards)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Laboratory standards file is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if (!is.null(peripheral)) {

    standards <- standards[
      purrr::map_lgl(
        standards$applicable_peripherals,
        ~ peripheral %in% .x
      ),
      ,
      drop = FALSE
    ]
  }

  if (!is.null(standard_ids)) {

    standard_ids <- as.character(
      standard_ids
    )

    missing_ids <- setdiff(
      standard_ids,
      standards$standard_id
    )

    if (length(missing_ids) > 0) {
      stop(
        "Requested standards not found: ",
        paste(
          missing_ids,
          collapse = ", "
        )
      )
    }

    standards <- standards[
      standards$standard_id %in% standard_ids,
      ,
      drop = FALSE
    ]
  }

  if (anyDuplicated(standards$standard_id)) {
    stop(
      "Duplicate standard IDs found in laboratory standards."
    )
  }

  if (verbose) {
    message(
      "Loaded ",
      nrow(standards),
      " laboratory standards."
    )
  }

  standards
}


#' Load IRMS experimental configuration
#'
#' Loads an experiment configuration from YAML.
#'
#' @param yaml_file_path Path to the experiment configuration YAML file.
#' @param stds_reference_df Optional pre-loaded laboratory standards table.
#' @param verbose Print loading information.
#' @return An object of class `irms_config`.
#' @export
load_IRMS_config <- function(
    yaml_file_path,
    stds_reference_df = NULL,
    verbose = TRUE
) {

  if (!file.exists(yaml_file_path)) {
    stop(
      "Configuration file not found: ",
      yaml_file_path
    )
  }

  config <- yaml::read_yaml(
    yaml_file_path
  )

  # Generic experiment information
  if (is.null(config$experiment)) {
    stop(
      "Configuration must contain an 'experiment' section."
    )
  }

  if (is.null(config$experiment$name)) {
    stop(
      "Configuration must contain experiment$name."
    )
  }

  if (is.null(config$experiment$file_name)) {
    stop(
      "Configuration must contain experiment$file_name."
    )
  }

  # Standards used in the experiment
  if (is.null(config$standards)) {
    stop(
      "Configuration must contain a 'standards' section."
    )
  }

  if (is.null(config$standards$used)) {
    stop(
      "Configuration must contain standards$used."
    )
  }

  stds_used <- unlist(
    config$standards$used,
    use.names = FALSE
  )

  if (length(stds_used) == 0) {
    stop(
      "standards$used must contain at least one standard."
    )
  }

  # Calibration standard selection
  if (
    is.null(
      config$standards$calibration_standards
    )
  ) {
    stop(
      "Configuration must contain ",
      "standards$calibration_standards."
    )
  }

  calibration_standards <- unlist(
    config$standards$calibration_standards
  )

  if (
    length(calibration_standards) == 0 ||
    is.null(names(calibration_standards))
  ) {
    stop(
      "standards$calibration_standards must be a named ",
      "mapping from element to standard ID."
    )
  }

  calibration_standards <- as.character(
    calibration_standards
  )

  if (any(!nzchar(calibration_standards))) {
    stop(
      "Calibration standard IDs cannot be empty."
    )
  }

  if (
    any(
      !calibration_standards %in% stds_used
    )
  ) {

    missing_calibration_standards <-
      setdiff(
        calibration_standards,
        stds_used
      )

    stop(
      "Calibration standards must also appear in ",
      "standards$used: ",
      paste(
        missing_calibration_standards,
        collapse = ", "
      )
    )
  }

  # Generic paths
  experiment_name <- as.character(
    config$experiment$name
  )

  raw_data_dir <- NULL

  output_dir <- NULL

  stds_data_path <- NULL

  if (!is.null(config$paths)) {

    raw_data_dir <- config$paths$raw_data_dir

    output_dir <- config$paths$output_dir

    stds_data_path <- config$paths$stds_data_path
  }

  if (is.null(raw_data_dir)) {

    raw_data_dir <- file.path(
      "data",
      experiment_name,
      "raw_data"
    )
  }

  if (is.null(output_dir)) {

    output_dir <- file.path(
      "data",
      experiment_name,
      "analyzed_data"
    )
  }

  # Locate laboratory standards
  if (is.null(stds_data_path)) {

    stds_data_path <- system.file(
      "config",
      "lab_standards.yml",
      package = "isotidy"
    )

    if (stds_data_path == "") {

      stds_data_path <- file.path(
        "inst",
        "config",
        "lab_standards.yml"
      )
    }
  }

  # Load laboratory standards if they were not supplied
  if (is.null(stds_reference_df)) {

    stds_reference_df <- load_lab_standards(
      path = stds_data_path,
      verbose = FALSE
    )
  }

  # Check calibration standards against reference database
  missing_standards <- setdiff(
    calibration_standards,
    stds_reference_df$standard_id
  )

  if (length(missing_standards) > 0) {

    stop(
      "Calibration standards not found in laboratory standards: ",
      paste(
        missing_standards,
        collapse = ", "
      )
    )
  }

  result <- list(
    experiment_name = experiment_name,
    file_name = as.character(
      config$experiment$file_name
    ),
    raw_data_dir = raw_data_dir,
    output_dir = output_dir,
    stds_data_path = stds_data_path,
    stds_used = stds_used,
    calibration_standards = calibration_standards,
    processing = config$processing,
    stds_reference_df = stds_reference_df,
    yaml_file_path = yaml_file_path
  )

  class(result) <- c(
    "irms_config",
    "list"
  )

  if (verbose) {

    message(
      "Loaded IRMS configuration for experiment: ",
      experiment_name
    )
  }

  result
}


#' Validate an IRMS configuration
#'
#' @param config An object returned by load_IRMS_config().
#' @return Invisibly returns TRUE.
#' @export
validate_config <- function(
    config
) {

  if (!inherits(
    config,
    "irms_config"
  )) {
    stop(
      "config must be an object created by load_IRMS_config()."
    )
  }

  if (
    is.null(config$experiment_name) ||
    !nzchar(config$experiment_name)
  ) {
    stop(
      "Configuration has no valid experiment name."
    )
  }

  if (
    is.null(config$file_name) ||
    !nzchar(config$file_name)
  ) {
    stop(
      "Configuration has no valid file name."
    )
  }

  if (
    is.null(config$stds_data_path) ||
    !file.exists(config$stds_data_path)
  ) {
    stop(
      "Laboratory standards file not found: ",
      config$stds_data_path
    )
  }

  if (
    is.null(config$stds_used) ||
    length(config$stds_used) == 0
  ) {
    stop(
      "No standards are specified in the configuration."
    )
  }

  if (
    is.null(config$calibration_standards) ||
    length(config$calibration_standards) == 0
  ) {
    stop(
      "No calibration standards are specified in the configuration."
    )
  }

  if (!dir.exists(config$output_dir)) {

    warning(
      "Output directory does not exist: ",
      config$output_dir
    )
  }

  invisible(TRUE)
}


#' Create an IRMS configuration template
#'
#' @param output_path Path for the YAML file to create.
#' @param experiment_name Experiment name.
#' @param file_name Experiment file name.
#' @param stds_used Standards used in the experiment.
#' @param calibration_standards Named mapping from element to calibration standard.
#' @param raw_data_dir Optional raw data directory.
#' @param output_dir Optional output directory.
#' @param stds_data_path Optional laboratory standards path.
#' @param overwrite Overwrite an existing file.
#' @return Invisibly returns the output path.
#' @export
create_config_template <- function(
    output_path,
    experiment_name,
    file_name,
    stds_used,
    calibration_standards,
    raw_data_dir = NULL,
    output_dir = NULL,
    stds_data_path = NULL,
    overwrite = FALSE
) {

  if (
    file.exists(output_path) &&
    !overwrite
  ) {
    stop(
      "Configuration file already exists: ",
      output_path
    )
  }

  if (is.null(raw_data_dir)) {

    raw_data_dir <- file.path(
      "data",
      experiment_name,
      "raw_data"
    )
  }

  if (is.null(output_dir)) {

    output_dir <- file.path(
      "data",
      experiment_name,
      "analyzed_data"
    )
  }

  config <- list(
    experiment = list(
      name = experiment_name,
      file_name = file_name
    ),
    standards = list(
      used = stds_used,
      calibration_standards = calibration_standards
    ),
    paths = list(
      raw_data_dir = raw_data_dir,
      output_dir = output_dir,
      stds_data_path = stds_data_path
    ),
    notes = ""
  )

  yaml::write_yaml(
    config,
    output_path
  )

  invisible(output_path)
}


#' Normalize raw peripheral data to the canonical format
#'
#' Applies a peripheral-specific adapter and validates the resulting
#' canonical data structure.
#'
#' @param raw_df Raw peripheral data.
#' @param adapter_function Function that converts raw data to canonical format.
#' @param ... Additional arguments passed to the adapter.
#' @return A canonical IRMS data frame.
#' @export
normalize_data <- function(
    raw_df,
    adapter_function,
    ...
) {

  if (!is.data.frame(raw_df)) {
    stop(
      "raw_df must be a data frame."
    )
  }

  if (!is.function(adapter_function)) {
    stop(
      "adapter_function must be a function."
    )
  }

  canonical_df <- adapter_function(
    raw_df,
    ...
  )

  required_columns <- c(
    "run_id",
    "sample_id",
    "element",
    "delta_value",
    "area_or_voltage",
    "amount",
    "instrument",
    "is_standard",
    "is_blank"
  )

  missing_columns <- setdiff(
    required_columns,
    names(canonical_df)
  )

  if (length(missing_columns) > 0) {

    stop(
      "Adapter did not produce required canonical columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  # Validate canonical numeric fields
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
      "Canonical numeric fields must be numeric: ",
      paste(
        non_numeric,
        collapse = ", "
      )
    )
  }

  # Validate canonical identifiers
  character_columns <- c(
    "run_id",
    "sample_id",
    "element",
    "instrument"
  )

  non_character <- character_columns[
    !vapply(
      canonical_df[character_columns],
      is.character,
      logical(1)
    )
  ]

  if (length(non_character) > 0) {

    stop(
      "Canonical identifier fields must be character: ",
      paste(
        non_character,
        collapse = ", "
      )
    )
  }

  # Validate canonical flags
  logical_columns <- c(
    "is_standard",
    "is_blank"
  )

  non_logical <- logical_columns[
    !vapply(
      canonical_df[logical_columns],
      is.logical,
      logical(1)
    )
  ]

  if (length(non_logical) > 0) {

    stop(
      "Canonical flag fields must be logical: ",
      paste(
        non_logical,
        collapse = ", "
      )
    )
  }

  canonical_df
}
