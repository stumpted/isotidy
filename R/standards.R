# Laboratory standards: load the reference table, flag standards in the data, merge reference values.

# Helper function to handle bulk and compound specific standards
normalize_yaml_standards <- function(records) {
  purrr::imap_dfr(records, function(record, record_name) {
    if (!is.list(record)) {
      stop("Standard '", record_name, "' must be a YAML mapping.")
    }

    standard_id <- as.character(record$standard_id %||% record_name)
    values <- record$compound_values %||% record$compound_delta_values

    if (is.null(values)) {
      compound_ids <- NA_character_
      values <- list(list(
        delta_value = record$delta_value,
        delta_uncertainty = record$delta_uncertainty
      ))
    } else if (is.null(names(values))) {
      # Support a sequence of one-item mappings, e.g.
      # - C14-ME: -29.98
      valid_entries <- vapply(
        values,
        function(x) is.list(x) && length(x) == 1 &&
          !is.null(names(x)) && nzchar(names(x)[[1]]),
        logical(1)
      )

      if (!all(valid_entries)) {
        stop(
          "Compound values for standard '", record_name,
          "' must map each compound ID to a value."
        )
      }

      compound_ids <- vapply(values, function(x) names(x)[[1]], character(1))
      values <- lapply(values, `[[`, 1)
    } else {
      # Support a named mapping, e.g. C14M: {delta_value: -29.98}
      compound_ids <- names(values)
    }

    scalar_standard <-
      length(compound_ids) == 1 &&
      length(values) == 1 &&
      is.na(compound_ids[[1]])

    invalid_compound_ids <-
      length(compound_ids) != length(values) ||
      (!scalar_standard && (
        anyNA(compound_ids) ||
        anyDuplicated(compound_ids)
      ))

    if (invalid_compound_ids) {
      stop("Compound IDs for standard '", record_name, "' must be unique.")
    }
    if (any(!is.na(compound_ids) & !nzchar(compound_ids))) {
      stop("Compound IDs for standard '", record_name, "' cannot be empty.")
    }

    peripheral_values <- as.character(
      unlist(record$applicable_peripherals, use.names = FALSE)
    )

    purrr::map2_dfr(values, compound_ids, function(value, compound_id) {
      value_metadata <- if (is.list(value)) value else list()
      if (is.list(value)) {
        delta_value <- value$delta_value
        delta_uncertainty <- value$delta_uncertainty
      } else {
        delta_value <- value
        delta_uncertainty <- NULL
      }

      delta_value <- suppressWarnings(as.numeric(delta_value))
      delta_uncertainty <- if (is.null(delta_uncertainty)) {
        NA_real_
      } else {
        suppressWarnings(as.numeric(delta_uncertainty))
      }

      if (length(delta_value) != 1 || is.na(delta_value)) {
        stop(
          "Standard '", record_name,
          "' must have a numeric delta value for every entry."
        )
      }

      if (length(delta_uncertainty) != 1) {
        stop(
          "Standard '", record_name,
          "' uncertainty must be a single number or null."
        )
      }

      aa_c_count <- if (is.null(value_metadata$AA_C_count)) {
        NA_real_
      } else {
        suppressWarnings(as.numeric(value_metadata$AA_C_count))
      }
      nacme_c_count <- if (is.null(value_metadata$NACME_C_count)) {
        NA_real_
      } else {
        suppressWarnings(as.numeric(value_metadata$NACME_C_count))
      }
      if (
        length(aa_c_count) != 1 || length(nacme_c_count) != 1 ||
        (!is.null(value_metadata$AA_C_count) && is.na(aa_c_count)) ||
        (!is.null(value_metadata$NACME_C_count) && is.na(nacme_c_count))
      ) {
        stop(
          "Standard '", record_name,
          "' carbon counts must be single numeric values or null."
        )
      }
      if (
        (!is.na(aa_c_count) && (!is.finite(aa_c_count) || aa_c_count <= 0)) ||
        (!is.na(nacme_c_count) &&
          (!is.finite(nacme_c_count) || nacme_c_count <= 0))
      ) {
        stop(
          "Standard '", record_name,
          "' carbon counts must be finite positive values."
        )
      }

      tibble::tibble(
        element = as.character(record$element %||% NA_character_),
        delta_value = delta_value,
        delta_uncertainty = delta_uncertainty,
        scale = as.character(record$scale %||% NA_character_),
        element_fraction = if (is.null(record$element_fraction)) {
          NA_real_
        } else {
          as.numeric(record$element_fraction)
        },
        standard_id = standard_id,
        compound_id = compound_id,
        AA_C_count = aa_c_count,
        NACME_C_count = nacme_c_count,
        applicable_peripherals = list(peripheral_values)
      )
    })
  })
}

#' Load laboratory standards
#'
#' Loads and validates scalar or compound-specific laboratory standards.
#'
#' @param path Path to the standards YAML or CSV file.
#' @param peripheral Optional peripheral name used to filter applicable standards.
#' @param standard_ids Optional standard IDs to load.
#' @param verbose Print loading information.
#' @return A data frame of laboratory standard reference values.
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

  if (grepl("\\.ya?ml$", path, ignore.case = TRUE)) {

    records <- yaml::read_yaml(path)

    if (
      is.null(records) ||
      length(records) == 0
    ) {
      stop(
        "Laboratory standards file is empty: ",
        path
      )
    }

    standards <- normalize_yaml_standards(records)

  } else if (grepl(
    "\\.csv$",
    path,
    ignore.case = TRUE
  )) {

    standards <- readr::read_csv(
      path,
      show_col_types = FALSE
    )

    # Existing CSV reference files describe scalar standards.
    if (!"compound_id" %in% names(standards)) {
      standards$compound_id <- NA_character_
    }

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
    "compound_id",
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

  if (
    anyNA(standards$standard_id) ||
    any(!nzchar(as.character(standards$standard_id)))
  ) {
    stop(
      "Laboratory standards must have non-empty standard IDs."
    )
  }

  if (anyNA(standards$delta_value)) {
    stop(
      "Laboratory standards must have a numeric delta_value ",
      "for every standard or compound."
    )
  }

  # A standard ID may appear on multiple rows when each row has a
  # different compound ID. Duplicate standard/compound pairs are invalid.
  if (anyDuplicated(standards[c("standard_id", "compound_id")])) {
    stop(
      "Duplicate standard_id + compound_id entries found ",
      "in laboratory standards."
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

  if (verbose) {
    message(
      "Loaded ",
      nrow(standards),
      " laboratory standard reference values."
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

