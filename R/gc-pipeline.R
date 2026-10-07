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
#' @param external_standard_id Optional sample ID for the external-standard
#'   injections. Defaults to config$gc_standards$external$sample_id.
#' @param external_reference_id Optional laboratory reference ID. Defaults to
#'   config$gc_standards$external$reference_id, then external_standard_id.
#' @param exclude_standard_injection_ids Injection IDs to omit when calculating
#'   external-standard offsets. Defaults to
#'   config$gc_standards$external$exclude_injection_ids.
#' @param exclude_standard_compound_ids Compound IDs to omit from all external
#'   standard injection and daily offset averages. Defaults to
#'   config$gc_standards$external$exclude_compound_ids.
#' @param exclude_measurements Optional data frame of additional manual
#'   exclusions. It is combined with
#'   `config$processing$gc$exclude_measurements`. Requires `analysis_no`;
#'   optional `peak_number` restricts an exclusion to one peak, while a
#'   missing/NA peak number excludes the entire analysis. Optional `reason`
#'   values are retained for the audit.
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
    external_standard_id = NULL,
    external_reference_id = NULL,
    isotope = "C13",
    verbose = TRUE,
    exclude_standard_injection_ids = NULL,
    exclude_measurements = NULL,
    exclude_standard_compound_ids = NULL
) {
  if (is.null(stds_reference_df)) {
    stds_reference_df <- config$stds_reference_df
  }

  if (!is.data.frame(stds_reference_df)) {
    stop("'stds_reference_df' must be a data frame.")
  }

  external_role <- config$gc_standards$external

  if (is.null(external_standard_id)) {
    external_standard_id <- external_role$sample_id
  }

  if (is.null(external_reference_id)) {
    external_reference_id <- external_role$reference_id %||%
      external_standard_id
  }

  if (is.null(exclude_standard_injection_ids)) {
    exclude_standard_injection_ids <-
      external_role$exclude_injection_ids
  }
  if (is.null(exclude_standard_compound_ids)) {
    exclude_standard_compound_ids <- external_role$exclude_compound_ids
  }

  if (
    length(external_standard_id) != 1 ||
    is.na(external_standard_id) ||
    !nzchar(external_standard_id)
  ) {
    stop(
      "Supply external_standard_id or configure ",
      "standards$gc$external$sample_id."
    )
  }

  if (
    length(external_reference_id) != 1 ||
    is.na(external_reference_id) ||
    !nzchar(external_reference_id)
  ) {
    stop(
      "Supply external_reference_id or configure ",
      "standards$gc$external$reference_id."
    )
  }

  data <- adapt_gc_data(
    raw_df = raw_df,
    config = config,
    isotope = isotope,
    verbose = verbose,
    gc_methods = gc_methods,
    method_name = method_name
  )

  processing_config <- config[["processing"]]
  gc_processing_config <- if (is.list(processing_config)) {
    processing_config[["gc"]]
  } else {
    NULL
  }
  configured_exclusions <- if (is.list(gc_processing_config)) {
    gc_processing_config[["exclude_measurements"]]
  } else {
    NULL
  }

  # Some legacy config objects flatten processing settings. Read this nested
  # GC-only setting from the source YAML when the in-memory structure lacks it.
  yaml_file_path <- config[["yaml_file_path"]]
  if (
    is.null(configured_exclusions) &&
    length(yaml_file_path) == 1 &&
    !is.na(yaml_file_path) &&
    nzchar(yaml_file_path) &&
    file.exists(yaml_file_path)
  ) {
    yaml_config <- yaml::read_yaml(yaml_file_path)
    yaml_processing <- yaml_config[["processing"]]
    yaml_gc_processing <- if (is.list(yaml_processing)) {
      yaml_processing[["gc"]]
    } else {
      NULL
    }
    configured_exclusions <- if (is.list(yaml_gc_processing)) {
      yaml_gc_processing[["exclude_measurements"]]
    } else {
      NULL
    }
  }
  if (is.null(configured_exclusions) || length(configured_exclusions) == 0) {
    configured_exclusions <- data.frame(
      analysis_no = character(),
      peak_number = character(),
      reason = character()
    )
  } else if (!is.data.frame(configured_exclusions)) {
    if (!is.list(configured_exclusions)) {
      stop(
        "config$processing$gc$exclude_measurements must be a list of ",
        "analysis/peak selectors."
      )
    }
    if ("analysis_no" %in% names(configured_exclusions)) {
      configured_exclusions <- list(configured_exclusions)
    }
    configured_exclusions <- lapply(configured_exclusions, function(selector) {
      if (is.data.frame(selector)) {
        return(selector)
      }
      if (!is.list(selector) || is.null(selector$analysis_no)) {
        stop(
          "Each configured GC exclusion must define an 'analysis_no'."
        )
      }
      as.data.frame(selector, stringsAsFactors = FALSE)
    })
    configured_exclusions <- dplyr::bind_rows(configured_exclusions)
  }

  if (!is.null(exclude_measurements) && !is.data.frame(exclude_measurements)) {
    stop("'exclude_measurements' must be a data frame or NULL.")
  }
  if (is.null(exclude_measurements)) {
    exclude_measurements <- configured_exclusions
  } else if (nrow(configured_exclusions) > 0) {
    exclude_measurements <- dplyr::bind_rows(
      configured_exclusions,
      exclude_measurements
    )
  }

  # Apply the user's explicit analysis/peak selection before calculating
  # standard offsets. A missing peak number excludes the whole analysis.
  unmatched_measurement_exclusions <- data.frame()
  excluded_measurement_rows <- integer(0)
  exclusion_reasons <- rep(NA_character_, nrow(data))

  if (!is.null(exclude_measurements)) {
    if (!"analysis_no" %in% names(exclude_measurements)) {
      stop("'exclude_measurements' must contain an 'analysis_no' column.")
    }
    if (!"peak_number" %in% names(exclude_measurements)) {
      exclude_measurements$peak_number <- NA
    }
    if (!"reason" %in% names(exclude_measurements)) {
      exclude_measurements$reason <- "manually excluded"
    }

    selector_analysis_numbers <- trimws(as.character(
      exclude_measurements$analysis_no
    ))
    selector_peaks <- as.character(exclude_measurements$peak_number)
    selector_reasons <- trimws(as.character(exclude_measurements$reason))
    selector_reasons[is.na(selector_reasons) | !nzchar(selector_reasons)] <-
      "manually excluded"
    unmatched_selector_rows <- integer(0)

    for (i in seq_len(nrow(exclude_measurements))) {
      if (
        is.na(selector_analysis_numbers[[i]]) ||
        !nzchar(selector_analysis_numbers[[i]])
      ) {
        stop("Every manual exclusion must have a non-empty analysis_no.")
      }
      matching_rows <- which(
        !is.na(data$analysis_no) &
          trimws(as.character(data$analysis_no)) ==
            selector_analysis_numbers[[i]]
      )
      if (!is.na(selector_peaks[[i]])) {
        matching_rows <- matching_rows[
          !is.na(data$peak_number[matching_rows]) &
            as.character(data$peak_number[matching_rows]) ==
              selector_peaks[[i]]
        ]
      }
      if (length(matching_rows) > 0) {
        excluded_measurement_rows <- union(
          excluded_measurement_rows,
          matching_rows
        )
        exclusion_reasons[matching_rows] <- selector_reasons[[i]]
      } else {
        unmatched_selector_rows <- c(unmatched_selector_rows, i)
      }
    }

    unmatched_measurement_exclusions <- exclude_measurements[
      unmatched_selector_rows,
      ,
      drop = FALSE
    ]
    if (nrow(unmatched_measurement_exclusions) > 0) {
      warning(
        nrow(unmatched_measurement_exclusions),
        " manual analysis exclusion selector(s) did not match GC data."
      )
    }
  }

  excluded_measurements <- data[
    sort(excluded_measurement_rows),
    ,
    drop = FALSE
  ]
  excluded_measurements$exclusion_reason <-
    exclusion_reasons[sort(excluded_measurement_rows)]
  if (length(excluded_measurement_rows) > 0) {
    data <- data[-excluded_measurement_rows, , drop = FALSE]
  }

  correction <- correct_gc_external_standard(
    canonical_df = data,
    stds_reference_df = stds_reference_df,
    standard_id = external_reference_id,
    sample_id = external_standard_id,
    exclude_injection_ids = exclude_standard_injection_ids,
    exclude_compound_ids = exclude_standard_compound_ids,
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

  visualization_data <- summarize_replicates(
    data = output,
    peripheral = "GC",
    config = config,
    verbose = verbose
  )

  if (verbose) {
    message(
      "GC processing complete: ",
      nrow(output),
      " peaks; external standard ",
      external_standard_id,
      " (reference ",
      external_reference_id,
      "); ",
      nrow(correction$daily_offsets),
      " corrected date(s)."
    )
  }

  list(
    data = data,
    output = output,
    visualization_data = visualization_data,
    excluded_measurements = excluded_measurements,
    unmatched_measurement_exclusions = unmatched_measurement_exclusions,
    excluded_standard_injection_ids =
      correction$excluded_standard_injection_ids,
    missing_excluded_standard_injection_ids =
      correction$missing_excluded_standard_injection_ids,
    excluded_standard_compound_ids =
      correction$excluded_standard_compound_ids,
    missing_excluded_standard_compound_ids =
      correction$missing_excluded_standard_compound_ids,
    injection_offsets = correction$injection_offsets,
    daily_offsets = correction$daily_offsets,
    external_standard_id = external_standard_id,
    external_reference_id = external_reference_id,
    method_name = method_name,
    isotope = isotope
  )
}
