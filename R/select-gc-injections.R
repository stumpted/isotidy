#' Select GC injections
#'
#' Selects complete GC injections from canonical GC-IRMS data.
#' All peaks belonging to a selected injection are retained.
#'
#' @param canonical_df Canonical GC-IRMS data frame.
#' @param injection_ids Optional vector of injection IDs to select.
#' @param sample_ids Optional vector of sample IDs to select.
#' @param verbose Print a selection summary.
#'
#' @return A canonical GC-IRMS data frame containing selected injections.
#' @export
select_gc_injections <- function(
    canonical_df,
    injection_ids = NULL,
    sample_ids = NULL,
    verbose = TRUE
) {

  if (!is.data.frame(canonical_df)) {
    stop(
      "'canonical_df' must be a data frame."
    )
  }

  required_columns <- c(
    "injection_id",
    "sample_id",
    "comment"
  )

  missing_columns <- setdiff(
    required_columns,
    names(canonical_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Canonical GC data are missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (
    is.null(injection_ids) &&
    is.null(sample_ids)
  ) {
    stop(
      "At least one of 'injection_ids' or 'sample_ids' ",
      "must be provided."
    )
  }

  selected_injections <- character(0)

  if (!is.null(injection_ids)) {

    matching_injections <- canonical_df$injection_id[
      canonical_df$injection_id %in% injection_ids
    ]

    selected_injections <- c(
      selected_injections,
      as.character(matching_injections)
    )
  }

  if (!is.null(sample_ids)) {

    matching_injections <- canonical_df$injection_id[
      canonical_df$sample_id %in% sample_ids
    ]

    selected_injections <- c(
      selected_injections,
      as.character(matching_injections)
    )
  }

  selected_injections <- unique(
    selected_injections
  )

  if (length(selected_injections) == 0) {

    warning(
      "No injections matched the supplied selection criteria."
    )

    return(
      canonical_df[0, , drop = FALSE]
    )
  }

  selected <- canonical_df[
    as.character(canonical_df$injection_id) %in%
      selected_injections,
    ,
    drop = FALSE
  ]

  if (verbose) {

    n_injections <- length(
      unique(selected$injection_id)
    )

    n_peaks <- nrow(selected)

    cat(
      "Selected ",
      n_injections,
      " injections containing ",
      n_peaks,
      " peaks.\n",
      sep = ""
    )

    cat(
      "\nSelected injections:\n"
    )

    injection_summary <- selected |>
      dplyr::distinct(
        injection_id,
        sample_id,
        comment
      )

    print(
      injection_summary
    )
  }

  selected
}
