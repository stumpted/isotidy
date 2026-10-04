# Generic entry point for converting raw peripheral data to the canonical format.

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
