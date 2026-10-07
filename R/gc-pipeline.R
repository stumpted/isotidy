# GC-IRMS pipeline: adapt peaks, assign compounds, and apply external standard corrections.

#' Process GC-IRMS data through external-standard correction
#'
#' Converts raw GC-IRMS peak data to canonical form, assigns compounds from
#' the selected GC method, and applies date-specific external-standard offsets.
#' This pipeline does not apply derivatization or internal-standard corrections.
#'
#' @param raw_df Raw GC-IRMS data frame.
#' @param config Experiment configuration loaded by load_IRMS_config().
#' @param stds_reference_df Compound-specific standard reference table. Defaults
#'   to config$stds_reference_df.
#' @param gc_methods GC methods YAML path or parsed list.
#' @param method_name GC method key in the methods configuration.
#' @param external_standard_id Sample ID and standard_id for the external
#'   standard, such as "F8".
#' @param isotope Isotope system, such as "C13".
#' @param verbose Print progress messages.
#'
#' @return A list with canonical data, final output, injection offsets, and
#'   daily offsets.
#' @export
process_gc <- function(
    raw_df,
    config,
    stds_reference_df = NULL,
    gc_methods,
    method_name,
    external_standard_id,
    isotope = "C13",
    verbose = TRUE
) {
  if (is.null(stds_reference_df)) {
    stds_reference_df <- config$stds_reference_df
  }

  if (!is.data.frame(stds_reference_df)) {
    stop("'stds_reference_df' must be a data frame.")
  }

  if (
    length(external_standard_id) != 1 ||
    is.na(external_standard_id) ||
    !nzchar(external_standard_id)
  ) {
    stop("'external_standard_id' must be one non-empty value.")
  }

  data <- adapt_gc_data(
    raw_df = raw_df,
    config = config,
    isotope = isotope,
    verbose = verbose,
    gc_methods = gc_methods,
    method_name = method_name
  )

  correction <- correct_gc_external_standard(
    canonical_df = data,
    stds_reference_df = stds_reference_df,
    standard_id = external_standard_id,
    verbose = verbose
  )

  data <- correction$data
  data$run_id <- config$experiment_name
  data$method_name <- method_name
  data$is_external_standard <- !is.na(data$sample_id) &
    trimws(as.character(data$sample_id)) == external_standard_id
  data$delta_value_final <- data$delta_value_external_corrected

  output_cols <- c(
    "run_id",
    "method_name",
    "injection_id",
    "sample_id",
    "analysis_no",
    "start_time",
    "stop_time",
    "element",
    "isotope",
    "compound_id",
    "compound_match_status",
    "is_external_standard",
    "peak_number",
    "retention_time",
    "delta_value",
    "delta_value_external_corrected",
    "delta_value_final",
    "peak_amplitude",
    "area_or_voltage",
    "instrument"
  )
  output_cols <- output_cols[output_cols %in% names(data)]
  output <- data[, output_cols, drop = FALSE]

  if (verbose) {
    message(
      "GC processing complete: ",
      nrow(output),
      " peaks; external standard ",
      external_standard_id,
      "; ",
      nrow(correction$daily_offsets),
      " corrected date(s)."
    )
  }

  list(
    data = data,
    output = output,
    injection_offsets = correction$injection_offsets,
    daily_offsets = correction$daily_offsets,
    external_standard_id = external_standard_id,
    method_name = method_name,
    isotope = isotope
  )
}
