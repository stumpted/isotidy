# Laboratory standards: load the reference table, flag standards in the data, merge reference values.

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


identify_standards <- function(df, standard_names, verbose = TRUE) {

  if (verbose) cat("Identifying standards\n")

  if (verbose) {
    cat("  Standards to find:", paste(standard_names, collapse = ", "), "\n")
  }

  if (!is.character(standard_names) || length(standard_names) == 0) {
    stop("standard_names must be non-empty character vector")
  }

  df <- df %>%
    dplyr::mutate(is_standard = .data$sample_id %in% standard_names)

  n_stds <- sum(df$is_standard, na.rm = TRUE)

  if (verbose) {
    cat("  Found", n_stds, "standard measurements\n")
    cat("  ✓ Standards identified\n\n")
  }

  return(df)
}


#' Extract standards and merge with reference values
#'
#' Isolates standards from dataset and enriches them with certified
#' isotope values from reference database. This merged dataframe is
#' used to build calibration curves.
#'
#' @param df Canonical dataframe with is_standard column
#' @param stds_reference_df Reference database with certified values
#'   Must have at least: standard_id, delta_value, element_fraction, etc.
#'   Other columns (std.d15N.AIR, standard_type, etc.) preserved
#' @param config Configuration (optional, for element, etc.)
#' @param verbose Print progress (default: TRUE)
#'
#' @return Dataframe containing ONLY standards with reference values merged
#'   Columns include:
#'   - All canonical columns (run_id, sample_id, delta_value, etc.)
#'   - Reference columns (delta_value, element_fraction, etc.)
#'   - Warnings if standards not found in reference DB
#'
#' @details
#' Merge is by:
#' - canonical: sample_id
#' - reference: standard_id
#'
#' A warning is issued for unmatched standards.
#' These unmatched rows will have NA in reference columns.
#' Consider this an error to investigate.
#'
extract_standards <- function(df, stds_reference_df, config = NULL, verbose = TRUE) {

  if (verbose) cat("Extracting standards\n")

  stds_df <- df %>%
    dplyr::filter(.data$is_standard == TRUE)

  if (verbose) {
    cat("  Standards to process:", nrow(stds_df), "measurements\n")
  }

  if (!is.data.frame(stds_reference_df)) {
    stop("stds_reference_df must be a dataframe")
  }

  if (!"standard_id" %in% colnames(stds_reference_df)) {
    stop("stds_reference_df must have 'standard_id' column")
  }

  # Rename so measured and reference deltas stay distinct
  stds_reference_df <- stds_reference_df %>%
    dplyr::rename(delta_value_reference = "delta_value") %>%
    dplyr::select(-dplyr::all_of("element"))

  stds_df <- dplyr::left_join(
    stds_df,
    stds_reference_df,
    by = c("sample_id" = "standard_id"),
    relationship = "many-to-one"
  )

  unmatched <- sum(is.na(stds_df$delta_value_reference))

  if (unmatched > 0) {
    unmatched_names <- stds_df %>%
      dplyr::filter(is.na(.data$delta_value_reference)) %>%
      dplyr::pull("sample_id") %>%
      unique()

    warning(paste(
      unmatched, "standard measurements not found in reference database:",
      paste(unmatched_names, collapse = ", ")
    ))
  }

  if (verbose) {
    cat("  Merged with reference database\n")
    if (unmatched > 0) {
      cat("  ⚠ Warning:", unmatched, "unmatched standards\n")
    }
    cat("  ✓ Standards extracted\n\n")
  }

  return(stds_df)
}
