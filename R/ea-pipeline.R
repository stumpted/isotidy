# EA-IRMS pipeline: orchestrates adapt -> clean -> standards/blanks -> mass -> blank -> linearity -> scale -> output.

process_ea <- function(
    raw_df,
    config,
    schema,
    stds_reference_df,
    isotope = "C13",
    run_id = NULL,
    instrument = "EA",
    apply_mass = TRUE,
    apply_linearity = TRUE,
    verbose = TRUE
) {

  element <- get_config_element(config)

  # 1. Adapt raw EA data to the canonical format
  data <- adapt_ea_data(
    raw_df = raw_df,
    config = config,
    schema = schema,
    isotope = isotope,
    run_id = run_id,
    instrument = instrument,
    verbose = verbose
  )

  # 2. Clean canonical data
  data <- clean_data(
    df = data,
    config = config,
    verbose = verbose
  )

  # 3. Identify standards and blanks
  standard_names <- config$stds_used

  data <- identify_standards(
    df = data,
    standard_names = standard_names,
    verbose = verbose
  )

  data <- identify_blanks(
    df = data,
    blank_identifier = "Blank",
    verbose = verbose
  )

  # 4. Extract standards and merge reference values
  stds_df <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  # 5. Mass calibration
  mass_stds <- prepare_mass_standards(
    stds_df = stds_df,
    element = element,
    config = config,
    verbose = verbose
  )

  mass_model <- build_mass_model(
    stds_df = mass_stds,
    element = element,
    config = config,
    verbose = verbose,
    create_plot = TRUE
  )

  if (apply_mass) {
    data <- apply_mass_correction(
      df = data,
      mass_model = mass_model,
      config = config,
      verbose = verbose
    )
  }

  # 6. Blank correction
  blank_df <- data %>%
    dplyr::filter(.data$is_blank)

  blank_stats <- calculate_blank_stats(
    blank_df = blank_df,
    verbose = verbose
  )

  data <- apply_blank_correction(
    df = data,
    blank_stats = blank_stats,
    method = "ea_area_subtraction",
    verbose = verbose
  )

  # 7. Re-extract standards after blank correction
  stds_df <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  # 8. Linearity calibration
  linearity_model <- build_linearity_model(
    stds_df = stds_df,
    config = config,
    verbose = verbose,
    create_plot = TRUE
  )

  if (apply_linearity) {
    data <- apply_linearity_correction(
      df = data,
      linearity_model = linearity_model,
      verbose = verbose
    )
  }

  # 9. Re-extract standards after linearity correction
  stds_df <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  # 10. Scale calibration (uses the enriched standards df; no separate reference table)
  scale_delta_column <- if (
    apply_linearity &&
    "delta_value_linearity_corrected" %in% colnames(stds_df)
  ) {
    "delta_value_linearity_corrected"
  } else if (
    "delta_value_blank_corrected" %in% colnames(stds_df)
  ) {
    "delta_value_blank_corrected"
  } else {
    "delta_value"
  }

  scale_model <- build_scale_model(
    stds_df = stds_df,
    delta_column = scale_delta_column,
    reference_column = "delta_value_reference",
    config = config,
    verbose = verbose
  )

  # 11. Apply final scale correction
  data <- apply_scale_correction(
    df = data,
    scale_model = scale_model,
    verbose = verbose
  )

  # 12. Re-extract standards after all corrections
  final_standards <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  # 13. Prepare final output ({element}_amount / {element}_percent column names)
  output_cols <- c(
    "run_id",
    "sample_id",
    "element",
    "amount",
    "area_or_voltage",
    paste0(element, "_amount"),
    paste0(element, "_percent"),
    "delta_value",
    "delta_value_blank_corrected",
    "delta_value_linearity_corrected",
    "delta_value_final"
  )

  output_cols <- output_cols[
    output_cols %in% colnames(data)
  ]

  final_output <- prepare_output(
    df = data,
    output_cols = output_cols,
    verbose = verbose
  )

  visualization_data <- summarize_replicates(
    data = final_output,
    peripheral = "EA",
    config = config,
    verbose = verbose
  )

  # 14. Return results
  models <- list(
    mass = mass_model,
    linearity = linearity_model,
    scale = scale_model
  )

  return(
    list(
      data = data,
      output = final_output,
      visualization_data = visualization_data,
      standards = final_standards,
      mass_standards = mass_stds,
      mass_model = mass_model,
      blank_stats = blank_stats,
      linearity_model = linearity_model,
      scale_model = scale_model,
      models = models,
      apply_mass = apply_mass,
      apply_linearity = apply_linearity
    )
  )
}
