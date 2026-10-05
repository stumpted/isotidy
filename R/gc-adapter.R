adapt_gc_data <- function(
    raw_df,
    config = NULL,
    isotope = "C13",
    verbose = TRUE
) {

  excel_datetime <- function(x) {
    x <- as.numeric(x)

    as.POSIXct(
      x * 86400,
      origin = "1899-12-30",
      tz = "UTC"
    )
  }

  isotope_system <- get_isotope_system(
    isotope = isotope,
    config = config
  )

  gc_schema <- get_peripheral_schema(
    peripheral = "GC",
    config = config
  )

  element <- isotope_system$element
  isotope_mass <- isotope_system$isotope_mass
  reference_mass <- isotope_system$reference_mass

  resolve_column <- function(template) {
    template |>
      stringr::str_replace_all(
        fixed("{element}"),
        element
      ) |>
      stringr::str_replace_all(
        fixed("{isotope_mass}"),
        as.character(isotope_mass)
      ) |>
      stringr::str_replace_all(
        fixed("{reference_mass}"),
        as.character(reference_mass)
      )
  }

  column_mapping <- gc_schema$column_mapping

  isotope_ratio_column <- resolve_column(
    column_mapping$isotope_ratio
  )

  retention_time_column <- resolve_column(
    column_mapping$retention_time
  )

  peak_amplitude_column <- resolve_column(
    column_mapping$peak_amplitude
  )

  peak_area_column <- resolve_column(
    column_mapping$peak_area
  )

  required_columns <- c(
    column_mapping$injection_id,
    column_mapping$sample_id,
    column_mapping$analysis_no,
    column_mapping$comment,
    column_mapping$sample_type,
    column_mapping$status,
    column_mapping$start_time,
    column_mapping$stop_time,
    column_mapping$peak_number,
    isotope_ratio_column,
    retention_time_column,
    peak_amplitude_column,
    peak_area_column
  )

  missing_columns <- setdiff(
    required_columns,
    names(raw_df)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing required GC columns for isotope ",
      isotope,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data <- raw_df |>
    dplyr::transmute(
      injection_id = .data[[column_mapping$injection_id]],
      sample_id = .data[[column_mapping$sample_id]],
      analysis_no = .data[[column_mapping$analysis_no]],
      comment = .data[[column_mapping$comment]],
      sample_type = .data[[column_mapping$sample_type]],
      status = .data[[column_mapping$status]],
      start_time = excel_datetime(
        .data[[column_mapping$start_time]]
      ),
      stop_time = excel_datetime(
        .data[[column_mapping$stop_time]]
      ),
      peak_number = .data[[column_mapping$peak_number]],
      retention_time = as.numeric(
        .data[[retention_time_column]]
      ),
      delta_value = as.numeric(
        .data[[isotope_ratio_column]]
      ),
      peak_amplitude = as.numeric(
        .data[[peak_amplitude_column]]
      ),
      area_or_voltage = as.numeric(
        .data[[peak_area_column]]
      ),
      element = element,
      isotope = isotope,
      amount = NA_real_,
      instrument = "GC-IRMS"
    )

  if (verbose) {
    cat(
      "Adapted ",
      nrow(data),
      " GC rows for element ",
      element,
      " using isotope ",
      isotope,
      " and reference mass ",
      reference_mass,
      "\n",
      sep = ""
    )
  }

  data
}
