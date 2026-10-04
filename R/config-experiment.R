# Experiment configuration: load, validate and template the per-experiment YAML config.

#' Load IRMS experimental configuration
#'
#' Loads an experiment configuration from YAML and converts it into
#' the canonical `irms_config` structure used by the processing pipeline.
#'
#' @param yaml_file_path Path to the experiment configuration YAML file.
#' @param lab_config_dir Directory containing laboratory configuration files.
#' @param stds_reference_df Optional pre-loaded laboratory standards table.
#' @param verbose Print loading information.
#' @return An object of class `irms_config`.
#' @export
load_IRMS_config <- function(
    yaml_file_path,
    lab_config_dir = NULL,
    stds_reference_df = NULL,
    verbose = TRUE
) {

  # Check configuration file
  if (!file.exists(yaml_file_path)) {
    stop(
      "Configuration file not found: ",
      yaml_file_path
    )
  }

  config <- yaml::read_yaml(
    yaml_file_path
  )

  # Experiment information
  if (is.null(config$experiment)) {
    stop(
      "Configuration must contain an 'experiment' section."
    )
  }

  if (
    is.null(config$experiment$name) ||
    length(config$experiment$name) != 1 ||
    is.na(config$experiment$name) ||
    !nzchar(config$experiment$name)
  ) {
    stop(
      "Configuration must contain a valid experiment$name."
    )
  }

  if (
    is.null(config$experiment$file_name) ||
    length(config$experiment$file_name) != 1 ||
    is.na(config$experiment$file_name) ||
    !nzchar(config$experiment$file_name)
  ) {
    stop(
      "Configuration must contain a valid experiment$file_name."
    )
  }

  experiment_name <- as.character(
    config$experiment$name
  )

  file_name <- as.character(
    config$experiment$file_name
  )

  # Standards used
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

  stds_used <- as.character(
    stds_used
  )

  if (length(stds_used) == 0) {
    stop(
      "standards$used must contain at least one standard."
    )
  }

  if (any(is.na(stds_used)) || any(!nzchar(stds_used))) {
    stop(
      "standards$used cannot contain empty standard IDs."
    )
  }

  # Calibration standards: named character vector, element -> standard ID (must stay named)
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
    config$standards$calibration_standards,
    use.names = TRUE
  )

  if (length(calibration_standards) == 0) {
    stop(
      "standards$calibration_standards cannot be empty."
    )
  }

  calibration_elements <- names(
    calibration_standards
  )

  if (
    is.null(calibration_elements) ||
    any(!nzchar(calibration_elements))
  ) {
    stop(
      "standards$calibration_standards must be a named ",
      "mapping from element to standard ID."
    )
  }

  if (anyDuplicated(calibration_elements)) {
    stop(
      "Each element may have only one calibration standard. ",
      "Duplicate elements found: ",
      paste(
        unique(
          calibration_elements[
            duplicated(calibration_elements)
          ]
        ),
        collapse = ", "
      )
    )
  }

  calibration_standards <- as.character(
    calibration_standards
  )

  # Restore names dropped by as.character()
  names(calibration_standards) <- calibration_elements

  if (
    any(is.na(calibration_standards)) ||
    any(!nzchar(calibration_standards))
  ) {
    stop(
      "Calibration standard IDs cannot be empty."
    )
  }

  # Every calibration standard must also be listed in standards$used.
  missing_from_run <- setdiff(
    unname(calibration_standards),
    stds_used
  )

  if (length(missing_from_run) > 0) {
    stop(
      "Calibration standards must also appear in ",
      "standards$used: ",
      paste(
        missing_from_run,
        collapse = ", "
      )
    )
  }

  # Paths
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

    if (is.null(lab_config_dir)) {
      stop(
        "No laboratory configuration directory was provided. ",
        "Supply lab_config_dir or specify paths$stds_data_path."
      )
    }

    stds_data_path <- file.path(
      lab_config_dir,
      "lab_standards.yml"
    )
  }

  # Load laboratory standards
  if (is.null(stds_reference_df)) {

    stds_reference_df <- load_lab_standards(
      path = stds_data_path,
      verbose = FALSE
    )
  }

  # Check calibration standards against the reference table
  missing_reference_standards <- setdiff(
    unname(calibration_standards),
    stds_reference_df$standard_id
  )

  if (length(missing_reference_standards) > 0) {

    stop(
      "Calibration standards not found in laboratory standards: ",
      paste(
        missing_reference_standards,
        collapse = ", "
      )
    )
  }

  # Build the irms_config object
  result <- list(
    experiment_name = experiment_name,
    file_name = file_name,
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

  # Validate the finished object
  validate_config(
    result
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
#' Validates the canonical structure of an object returned by
#' load_IRMS_config().
#'
#' @param config An object returned by load_IRMS_config().
#' @return Invisibly returns TRUE.
#' @export
validate_config <- function(
    config
) {

  # Must be an irms_config object
  if (!inherits(
    config,
    "irms_config"
  )) {
    stop(
      "config must be an object created by load_IRMS_config()."
    )
  }

  # Experiment information
  if (
    is.null(config$experiment_name) ||
    length(config$experiment_name) != 1 ||
    is.na(config$experiment_name) ||
    !nzchar(config$experiment_name)
  ) {
    stop(
      "Configuration has no valid experiment name."
    )
  }

  if (
    is.null(config$file_name) ||
    length(config$file_name) != 1 ||
    is.na(config$file_name) ||
    !nzchar(config$file_name)
  ) {
    stop(
      "Configuration has no valid file name."
    )
  }

  # Paths
  if (
    is.null(config$stds_data_path) ||
    length(config$stds_data_path) != 1 ||
    is.na(config$stds_data_path) ||
    !nzchar(config$stds_data_path)
  ) {
    stop(
      "Configuration has no valid laboratory standards path."
    )
  }

  if (!file.exists(config$stds_data_path)) {
    stop(
      "Laboratory standards file not found: ",
      config$stds_data_path
    )
  }

  if (
    is.null(config$output_dir) ||
    length(config$output_dir) != 1 ||
    is.na(config$output_dir) ||
    !nzchar(config$output_dir)
  ) {
    stop(
      "Configuration has no valid output directory."
    )
  }

  # The output directory is allowed not to exist yet.
  if (!dir.exists(config$output_dir)) {

    warning(
      "Output directory does not exist: ",
      config$output_dir
    )
  }

  # Standards used
  if (
    is.null(config$stds_used) ||
    length(config$stds_used) == 0
  ) {
    stop(
      "No standards are specified in the configuration."
    )
  }

  if (
    !is.character(config$stds_used) ||
    any(is.na(config$stds_used)) ||
    any(!nzchar(config$stds_used))
  ) {
    stop(
      "stds_used must contain non-empty standard IDs."
    )
  }

  # Calibration standards: named character vector, element -> standard ID
  if (
    is.null(config$calibration_standards) ||
    length(config$calibration_standards) == 0
  ) {
    stop(
      "No calibration standards are specified in the configuration."
    )
  }

  if (!is.character(config$calibration_standards)) {
    stop(
      "calibration_standards must be a character vector."
    )
  }

  calibration_elements <- names(
    config$calibration_standards
  )

  if (
    is.null(calibration_elements) ||
    any(!nzchar(calibration_elements))
  ) {
    stop(
      "calibration_standards must be a named mapping ",
      "from element to standard ID."
    )
  }

  if (anyDuplicated(calibration_elements)) {
    stop(
      "Each element may have only one calibration standard. ",
      "Duplicate elements found: ",
      paste(
        unique(
          calibration_elements[
            duplicated(calibration_elements)
          ]
        ),
        collapse = ", "
      )
    )
  }

  if (
    any(
      is.na(config$calibration_standards)
    ) ||
    any(
      !nzchar(config$calibration_standards)
    )
  ) {
    stop(
      "Calibration standard IDs cannot be empty."
    )
  }

  # Every calibration standard must be in stds_used
  missing_from_run <- setdiff(
    unname(config$calibration_standards),
    config$stds_used
  )

  if (length(missing_from_run) > 0) {

    stop(
      "Calibration standards must also appear in ",
      "stds_used: ",
      paste(
        missing_from_run,
        collapse = ", "
      )
    )
  }

  # Laboratory standards reference table
  if (
    is.null(config$stds_reference_df) ||
    !is.data.frame(config$stds_reference_df)
  ) {
    stop(
      "Configuration must contain a laboratory standards ",
      "reference data frame."
    )
  }

  if (!"standard_id" %in% names(config$stds_reference_df)) {
    stop(
      "Laboratory standards reference data are missing ",
      "the 'standard_id' column."
    )
  }

  missing_reference_standards <- setdiff(
    unname(config$calibration_standards),
    config$stds_reference_df$standard_id
  )

  if (length(missing_reference_standards) > 0) {

    stop(
      "Calibration standards are missing from the laboratory ",
      "standards reference data: ",
      paste(
        missing_reference_standards,
        collapse = ", "
      )
    )
  }

  # Processing section
  if (
    is.null(config$processing) ||
    !is.list(config$processing)
  ) {
    stop(
      "Configuration must contain a processing section."
    )
  }

  if (
    is.null(config$processing$element) ||
    length(config$processing$element) != 1 ||
    is.na(config$processing$element) ||
    !nzchar(config$processing$element)
  ) {
    stop(
      "Configuration processing section must contain ",
      "a valid element."
    )
  }

  invisible(TRUE)
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
