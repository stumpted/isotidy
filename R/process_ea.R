# Process EA-IRMS data through the isotidy correction pipeline
#
# This function orchestrates the EA workflow using the core processing
# functions. Optional mass and linearity corrections can be applied
# automatically or decided interactively after reviewing their calibration plots.

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

  # Adapt raw EA data to the canonical format

  data <- adapt_ea_data(
    raw_df = raw_df,
    config = config,
    schema = schema,
    isotope = isotope,
    run_id = run_id,
    instrument = instrument,
    verbose = verbose
  )

  # Clean canonical data

  data <- clean_data(
    df = data,
    config = config,
    verbose = verbose
  )

  # Identify standards

  standard_names <- config$stds_used

  data <- identify_standards(
    df = data,
    standard_names = standard_names
  )

  # Identify blanks

  data <- identify_blanks(
    df = data,
    blank_identifier = "Blank"
  )

  # Extract standards

  stds_df <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  # Prepare standards for mass calibration

  mass_stds <- prepare_mass_standards(
    stds_df = stds_df,
    element = get_config_element(config),
    config = config,
    verbose = verbose
  )

  # Build mass calibration model

  mass_model <- build_mass_model(
    stds_df = mass_stds,
    element = get_config_element(config),
    config = config,
    verbose = verbose,
    create_plot = TRUE
  )

  # Apply mass correction

  if (apply_mass) {
    data <- apply_mass_correction(
      df = data,
      mass_model = mass_model,
      config = config,
      verbose = verbose
    )
  }

  # Calculate blank statistics

  blank_df <- data %>%
    dplyr::filter(is_blank)

  blank_stats <- calculate_blank_stats(
    blank_df = blank_df,
    verbose = verbose
  )

  # Apply blank correction

  data <- apply_blank_correction(
    df = data,
    blank_stats = blank_stats,
    method = "ea_area_subtraction",
    verbose = verbose
  )

  # Re-extract standards after blank correction

  stds_df <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  # Build linearity calibration model

  linearity_model <- build_linearity_model(
    stds_df = stds_df,
    config = config,
    verbose = TRUE,
    create_plot = TRUE
  )

  # Apply linearity correction

  if (apply_linearity) {
    data <- apply_linearity_correction(
      df = data,
      linearity_model = linearity_model,
      verbose = verbose
    )
  }

  # Re-extract standards after linearity correction

  stds_df <- extract_standards(
    df = data,
    stds_reference_df = stds_reference_df,
    config = config,
    verbose = verbose
  )

  return(
    list(
      data = data,
      standards = stds_df,
      mass_standards = mass_stds,
      mass_model = mass_model,
      linearity_model = linearity_model
    )
  )
}
