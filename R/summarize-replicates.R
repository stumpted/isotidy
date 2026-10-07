# Schema-driven replicate summaries for visualization ----------------------

#' Summarize analytical replicates for visualization
#'
#' Uses a peripheral schema's `replicates` settings to derive sample groups and
#' analytical replicate identifiers. Repeated rows within one analytical
#' replicate are first reduced to one value. The result contains group
#' identifiers, replicate counts, and summary statistics for final measurement
#' columns.
#'
#' @param data Processed canonical data containing `sample_id` and the selected
#'   measurement columns.
#' @param peripheral Peripheral name in the schema, such as `EA` or `GC`.
#' @param config Optional experiment configuration. Used to exclude configured
#'   standards and blanks.
#' @param schema Optional peripheral schema. Loaded from the package when
#'   omitted.
#' @param value_columns Measurement columns to summarize. Defaults to
#'   `delta_value_final`.
#' @param include_standards Include standard rows in the summary.
#' @param include_blanks Include blank rows in the summary.
#' @param include_unassigned Include rows without a compound assignment.
#' @param verbose Print a summary of the result.
#' @return A visualization-ready data frame with one row per sample group and
#'   compound where applicable.
#' @export
summarize_replicates <- function(
    data,
    peripheral,
    config = NULL,
    schema = NULL,
    value_columns = NULL,
    include_standards = FALSE,
    include_blanks = FALSE,
    include_unassigned = FALSE,
    verbose = TRUE
) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!"sample_id" %in% names(data)) {
    stop("'data' must contain a 'sample_id' column.")
  }
  if (
    length(peripheral) != 1 || is.na(peripheral) ||
    !nzchar(as.character(peripheral))
  ) {
    stop("'peripheral' must be one non-empty name.")
  }

  peripheral <- as.character(peripheral)
  processing_config <- config[["processing"]]
  if (!is.list(processing_config)) {
    processing_config <- list()
  }
  if (is.null(schema)) {
    schemas <- load_peripheral_schemas()
    schema <- get_peripheral_schema(schemas, peripheral)
  }

  if (
    !"isotope" %in% names(data) &&
    !is.null(processing_config$isotope)
  ) {
    data$isotope <- as.character(processing_config$isotope)
  }

  replicate_settings <- schema$replicates
  if (!is.list(replicate_settings)) {
    stop("Schema for peripheral '", peripheral, "' must define replicates.")
  }

  identify_by <- replicate_settings$identify_by
  if (
    is.null(identify_by) || length(identify_by) != 1 ||
    is.na(identify_by) || !nzchar(identify_by)
  ) {
    stop("Replicate schema for '", peripheral, "' must define identify_by.")
  }
  identify_by <- as.character(identify_by)

  if (is.null(value_columns)) {
    value_columns <- "delta_value_final"
  }

  if (
    length(value_columns) == 0 ||
    anyNA(value_columns) ||
    any(!nzchar(value_columns))
  ) {
    stop("Specify at least one valid measurement column in 'value_columns'.")
  }

  missing_values <- setdiff(value_columns, names(data))
  if (length(missing_values) > 0) {
    stop(
      "Measurement columns are missing from 'data': ",
      paste(missing_values, collapse = ", ")
    )
  }
  non_numeric_values <- value_columns[
    !vapply(data[value_columns], is.numeric, logical(1))
  ]
  if (length(non_numeric_values) > 0) {
    stop(
      "Measurement columns must be numeric: ",
      paste(non_numeric_values, collapse = ", ")
    )
  }

  data$sample_id <- as.character(data$sample_id)
  if (identical(identify_by, "name_pattern")) {
    replicate_pattern <- replicate_settings$pattern
    if (
      is.null(replicate_pattern) || length(replicate_pattern) != 1 ||
      is.na(replicate_pattern) || !nzchar(replicate_pattern)
    ) {
      stop("Name-based replicate schema must define a valid pattern.")
    }

    data$.visual_sample_id <- sub(
      replicate_pattern,
      "\\1",
      data$sample_id,
      perl = TRUE
    )
    data$.analytical_replicate_id <- data$sample_id
  } else {
    group_column <- replicate_settings$group_column
    if (
      is.null(group_column) || length(group_column) != 1 ||
      is.na(group_column) || !nzchar(group_column)
    ) {
      group_column <- "sample_id"
    }
    missing_identifiers <- setdiff(
      c(identify_by, group_column),
      names(data)
    )
    if (length(missing_identifiers) > 0) {
      stop(
        "Replicate identifiers required by the '", peripheral,
        "' schema are missing: ", paste(missing_identifiers, collapse = ", ")
      )
    }

    data$.visual_sample_id <- as.character(data[[group_column]])
    data$.analytical_replicate_id <- as.character(data[[identify_by]])
  }

  keep <- rep(TRUE, nrow(data))
  if (!include_standards) {
    if ("is_standard" %in% names(data)) {
      keep <- keep & (is.na(data$is_standard) | !data$is_standard)
    }
    if ("is_external_standard" %in% names(data)) {
      keep <- keep & (is.na(data$is_external_standard) | !data$is_external_standard)
    }

    standard_ids <- character()
    if (!is.null(config$stds_used)) {
      standard_ids <- c(standard_ids, as.character(config$stds_used))
    }
    if (!is.null(config$calibration_standards)) {
      standard_ids <- c(
        standard_ids,
        as.character(unname(config$calibration_standards))
      )
    }
    standard_ids <- unique(standard_ids[!is.na(standard_ids)])
    if (length(standard_ids) > 0) {
      keep <- keep & !data$sample_id %in% standard_ids
    }
  }

  if (!include_blanks) {
    if ("is_blank" %in% names(data)) {
      keep <- keep & (is.na(data$is_blank) | !data$is_blank)
    }
    if ("sample_type" %in% names(data)) {
      keep <- keep & !grepl(
        "blank",
        as.character(data$sample_type),
        ignore.case = TRUE
      )
    }

    blank_patterns <- character()
    if (!is.null(processing_config$blank_identifiers)) {
      blank_patterns <- c(
        blank_patterns,
        as.character(unlist(processing_config$blank_identifiers))
      )
    }
    schema_blank_pattern <- schema$metadata$blank_identifier_pattern
    if (!is.null(schema_blank_pattern)) {
      blank_patterns <- c(blank_patterns, as.character(schema_blank_pattern))
    }
    blank_patterns <- unique(blank_patterns[!is.na(blank_patterns) & nzchar(blank_patterns)])
    if (length(blank_patterns) > 0) {
      keep <- keep & !grepl(
        paste(blank_patterns, collapse = "|"),
        data$sample_id,
        ignore.case = TRUE
      )
    }
  }

  if (!include_unassigned && "compound_id" %in% names(data)) {
    keep <- keep & !is.na(data$compound_id) & nzchar(as.character(data$compound_id))
  }

  data <- data[keep, , drop = FALSE]

  group_columns <- intersect(
    c(
      "run_id",
      ".visual_sample_id",
      "element",
      "isotope",
      "method_name",
      "compound_id"
    ),
    names(data)
  )

  if (nrow(data) == 0) {
    summary <- data[FALSE, group_columns, drop = FALSE]
    names(summary)[names(summary) == ".visual_sample_id"] <- "sample_id"
    summary$n_replicates <- integer()
    for (column in value_columns) {
      summary[[paste0(column, "_mean")]] <- numeric()
      summary[[paste0(column, "_sd")]] <- numeric()
      summary[[paste0(column, "_se")]] <- numeric()
    }
    return(summary)
  }

  if (anyNA(data$.analytical_replicate_id) || any(!nzchar(data$.analytical_replicate_id))) {
    stop("Analytical replicate identifiers cannot be missing or empty.")
  }

  mean_or_na <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  sd_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) NA_real_ else stats::sd(x)
  }
  se_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) NA_real_ else stats::sd(x) / sqrt(length(x))
  }

  replicate_columns <- c(group_columns, ".analytical_replicate_id")
  replicate_values <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(replicate_columns))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(value_columns), mean_or_na),
      .groups = "drop"
    )

  primary_value <- value_columns[[1]]
  summary <- replicate_values |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_columns))) |>
    dplyr::summarise(
      n_replicates = sum(!is.na(.data[[primary_value]])),
      dplyr::across(
        dplyr::all_of(value_columns),
        list(mean = mean_or_na, sd = sd_or_na, se = se_or_na),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    ) |>
    dplyr::rename(sample_id = .visual_sample_id)

  if (verbose) {
    message(
      "Summarized ", nrow(data), " rows into ", nrow(summary),
      " visualization group(s) for peripheral '", peripheral, "'."
    )
  }

  summary
}
