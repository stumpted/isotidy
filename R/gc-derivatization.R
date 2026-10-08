# GC-IRMS amino-acid derivatization correction helpers.

#' Estimate date- and compound-specific NACME isotope values
#'
#' Estimate the empirical NACME contribution from externally corrected AA mix
#' injections and the known compound-specific AA isotope values.
#'
#' @param canonical_df GC data after external-standard correction.
#' @param standards_reference_df Laboratory standard references, including the
#'   AA and NACME carbon counts for each compound.
#' @param standard_id Reference standard ID for the AA mix.
#' @param sample_id Sample label used for AA mix injections. Defaults to
#'   `standard_id`.
#' @param external_standard_reference_id Reference ID for the external GC
#'   standard; its compounds are excluded from NACME reference validation.
#' @param date_column Date or timestamp column used as the derivatization
#'   batch key.
#' @param compound_column Compound identifier column.
#' @param delta_column Externally corrected isotope value column.
#'
#' @return A data frame with one empirical NACME estimate per date and compound.
#' @keywords internal
estimate_gc_nacme <- function(
    canonical_df,
    standards_reference_df,
    standard_id,
    sample_id = standard_id,
    external_standard_reference_id = NULL,
    date_column = "start_time",
    compound_column = "compound_id",
    delta_column = "delta_value_external_corrected"
) {
  required_data <- c(
    "injection_id", "sample_id", date_column, compound_column, delta_column
  )
  missing_data <- setdiff(required_data, names(canonical_df))
  if (length(missing_data) > 0) {
    stop("GC data are missing required columns: ",
         paste(missing_data, collapse = ", "))
  }

  required_reference <- c(
    "standard_id", "compound_id", "delta_value", "delta_uncertainty",
    "AA_C_count", "NACME_C_count"
  )
  missing_reference <- setdiff(required_reference, names(standards_reference_df))
  if (length(missing_reference) > 0) {
    stop("Standards reference data are missing required columns: ",
         paste(missing_reference, collapse = ", "))
  }
  if (length(standard_id) != 1 || is.na(standard_id) || !nzchar(standard_id)) {
    stop("'standard_id' must be one non-empty value.")
  }
  if (length(sample_id) != 1 || is.na(sample_id) || !nzchar(sample_id)) {
    stop("'sample_id' must be one non-empty value.")
  }
  if (
    !is.null(external_standard_reference_id) &&
      (length(external_standard_reference_id) != 1 ||
        is.na(external_standard_reference_id) ||
        !nzchar(external_standard_reference_id))
  ) {
    stop("'external_standard_reference_id' must be NULL or one non-empty value.")
  }
  if (
    !is.null(external_standard_reference_id) &&
      as.character(standard_id) ==
        as.character(external_standard_reference_id)
  ) {
    stop(
      "The NACME reference standard ('", standard_id,
      "') is also configured as the external standard. Configure the AA mix ",
      "reference ID separately from the external standard reference ID."
    )
  }

  references <- standards_reference_df |>
    dplyr::filter(
      as.character(.data$standard_id) == as.character(standard_id),
      if (is.null(external_standard_reference_id)) {
        TRUE
      } else {
        as.character(.data$standard_id) !=
          as.character(external_standard_reference_id)
      },
      !is.na(.data$compound_id)
    ) |>
    dplyr::transmute(
      compound_id = as.character(.data$compound_id),
      delta_known = as.numeric(.data$delta_value),
      sigma_known = as.numeric(.data$delta_uncertainty),
      AA_C_count = as.numeric(.data$AA_C_count),
      NACME_C_count = as.numeric(.data$NACME_C_count)
    )
  if (nrow(references) == 0) {
    stop("No compound-specific references found for AA mix standard '",
         standard_id, "'.")
  }
  if (anyDuplicated(references$compound_id)) {
    stop("AA mix reference data contain duplicate compound IDs for '",
         standard_id, "'.")
  }
  invalid_counts <- !is.finite(references$AA_C_count) |
    references$AA_C_count <= 0 |
    !is.finite(references$NACME_C_count) |
    references$NACME_C_count <= 0
  if (any(invalid_counts)) {
    stop("AA mix references are missing valid carbon counts for: ",
         paste(references$compound_id[invalid_counts], collapse = ", "))
  }

  nacme_input <- canonical_df
  if ("is_external_standard" %in% names(nacme_input)) {
    nacme_input <- nacme_input |>
      dplyr::filter(
        is.na(.data$is_external_standard) | !.data$is_external_standard
      )
  }

  standard_rows <- nacme_input |>
    dplyr::filter(
      !is.na(.data$sample_id),
      trimws(as.character(.data$sample_id)) == as.character(sample_id),
      !is.na(.data[[compound_column]]),
      !is.na(.data[[date_column]]),
      is.finite(as.numeric(.data[[delta_column]])),
      !is.na(.data$injection_id)
    ) |>
    dplyr::transmute(
      injection_id = as.character(.data$injection_id),
      correction_date = as.Date(.data[[date_column]]),
      compound_id = as.character(.data[[compound_column]]),
      measured_delta = as.numeric(.data[[delta_column]])
    )
  if (nrow(standard_rows) == 0) {
    stop("No usable AA mix injections found for sample ID '", sample_id, "'.")
  }

  injection_values <- standard_rows |>
    dplyr::group_by(.data$correction_date, .data$injection_id,
                    .data$compound_id) |>
    dplyr::summarise(
      delta_injection = mean(.data$measured_delta),
      .groups = "drop"
    )

  estimates <- injection_values |>
    dplyr::group_by(.data$correction_date, .data$compound_id) |>
    dplyr::summarise(
      delta_derivatized_mean = mean(.data$delta_injection),
      sigma_derivatized = if (dplyr::n() > 1) {
        stats::sd(.data$delta_injection)
      } else {
        NA_real_
      },
      n_standard_injections = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::left_join(references, by = "compound_id",
                     relationship = "many-to-one") |>
    dplyr::mutate(
      total_C_count = .data$AA_C_count + .data$NACME_C_count,
      delta_NACME = (
        .data$total_C_count * .data$delta_derivatized_mean -
          .data$AA_C_count * .data$delta_known
      ) / .data$NACME_C_count,
      sigma_NACME = sqrt(
        (.data$total_C_count / .data$NACME_C_count *
           .data$sigma_derivatized)^2 +
          (.data$AA_C_count / .data$NACME_C_count *
             .data$sigma_known)^2
      )
    ) |>
    dplyr::arrange(.data$correction_date, .data$compound_id)

  estimates
}

# Apply date- and compound-specific NACME estimates to sample peaks.
correct_gc_nacme <- function(
    canonical_df,
    nacme_estimates,
    sample_id,
    date_column = "start_time",
    delta_column = "delta_value_external_corrected"
) {
  required_data <- c("sample_id", date_column, "compound_id", delta_column)
  missing_data <- setdiff(required_data, names(canonical_df))
  if (length(missing_data) > 0) {
    stop("GC data are missing required columns: ",
         paste(missing_data, collapse = ", "))
  }
  required_estimates <- c(
    "correction_date", "compound_id", "delta_NACME", "sigma_NACME"
  )
  missing_estimates <- setdiff(required_estimates, names(nacme_estimates))
  if (length(missing_estimates) > 0) {
    stop("NACME estimates are missing required columns: ",
         paste(missing_estimates, collapse = ", "))
  }
  if (anyDuplicated(nacme_estimates[c("correction_date", "compound_id")])) {
    stop("NACME estimates must be unique by date and compound.")
  }

  data <- canonical_df
  data$correction_date <- as.Date(data[[date_column]])
  data$is_derivatization_standard <- !is.na(data$sample_id) &
    trimws(as.character(data$sample_id)) == as.character(sample_id)
  data$delta_value_nacme_corrected <- NA_real_
  data$delta_NACME <- NA_real_
  data$sigma_NACME <- NA_real_

  row_id <- seq_len(nrow(data))
  eligible <- !data$is_derivatization_standard &
    !is.na(data$compound_id) & !is.na(data$correction_date)
  if (any(eligible)) {
    matched <- data[eligible, , drop = FALSE] |>
      dplyr::mutate(.gc_row_id = row_id[eligible]) |>
      dplyr::left_join(
        nacme_estimates |>
          dplyr::select(
            correction_date, compound_id,
            nacme_estimate = delta_NACME,
            nacme_estimate_sd = sigma_NACME,
            AA_C_count, NACME_C_count
          ),
        by = c("correction_date", "compound_id"),
        relationship = "many-to-one"
      )
    valid <- is.finite(as.numeric(matched[[delta_column]])) &
      is.finite(matched$nacme_estimate)
    corrected <- rep(NA_real_, nrow(matched))
    corrected[valid] <- (
      (matched$AA_C_count[valid] + matched$NACME_C_count[valid]) *
        as.numeric(matched[[delta_column]][valid]) -
        matched$NACME_C_count[valid] * matched$nacme_estimate[valid]
    ) / matched$AA_C_count[valid]
    data$delta_value_nacme_corrected[matched$.gc_row_id] <- corrected
    data$delta_NACME[matched$.gc_row_id] <- matched$nacme_estimate
    data$sigma_NACME[matched$.gc_row_id] <- matched$nacme_estimate_sd
  }

  data$delta_value_final <- as.numeric(data[[delta_column]])
  is_external_standard <- if ("is_external_standard" %in% names(data)) {
    !is.na(data$is_external_standard) & data$is_external_standard
  } else {
    rep(FALSE, nrow(data))
  }
  is_sample <- !data$is_derivatization_standard &
    !is_external_standard
  data$delta_value_final[is_sample] <-
    data$delta_value_nacme_corrected[is_sample]
  data
}
