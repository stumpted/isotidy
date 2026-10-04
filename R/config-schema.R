# Peripheral schemas: load, retrieve and validate the YAML column-mapping schemas.

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
