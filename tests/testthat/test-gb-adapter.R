gb_test_config <- function() {
  list(
    experiment_name = "synthetic_gb_test",
    stds_used = c("NaHCO3", "ETH4"),
    processing = list(
      element = "C",
      sample_peak_number_low = 4,
      sample_peak_number_high = 11
    )
  )
}

gb_test_raw_data <- function() {
  data.frame(
    "Analysis" = 101:106,
    "Identifier 1" = c(
      "NaHCO3_1umol",
      "NaHCO3_5umol",
      "ETH4_5umol",
      "NaHCO3X_1umol",
      "on_off",
      "blank"
    ),
    "Identifier 2" = c(0.0922, 0.492, 0.5084, 0.15, NA, NA),
    "Comment" = c("mg", "mg", "mg", "mg", "", ""),
    "Peak Nr" = c(4, 11, 4, 4, 1, 1),
    "Sample Dilution" = rep(0, 6),
    "d 13C/12C" = c(-2.213, -1.804, -9.442, -3, -12.142, -11.836),
    "Area 44" = c(0.912, 2.402, 2.167, 1.2, 140.788, 147.21),
    "Area 45" = c(0.012, 0.031, 0.028, 0.014, 1.647, 1.64),
    "Area 46" = c(0.004, 0.009, 0.008, 0.004, 0.599, 0.597),
    "Ampl  44" = c(830, 2193, 1976, 1000, 7276, 7268),
    "Ampl  45" = c(900, 2400, 2100, 1100, 8510, 8500),
    "Ampl  46" = c(1000, 2700, 2300, 1200, 10320, 10308),
    "BGD 44" = rep(0.3, 6),
    "BGD 45" = rep(0.2, 6),
    "BGD 46" = rep(2.3, 6),
    "Date" = rep("08/26/26", 6),
    "Time" = rep("17:11:50", 6),
    "Rt" = rep(38.7, 6),
    "Width" = rep(27.5, 6),
    "Gasconfiguration" = rep("CO2", 6),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

gb_test_schema <- function() {
  schemas <- load_peripheral_schemas(
    config_dir = testthat::test_path("..", "..", "inst", "config")
  )
  get_peripheral_schema(schemas = schemas, peripheral = "GB")
}

test_that("GB adapter maps peak rows and retains mass-specific signals", {
  raw <- gb_test_raw_data()
  adapted <- adapt_gb_data(
    raw_df = raw,
    config = gb_test_config(),
    schema = gb_test_schema(),
    isotope = "C13",
    verbose = FALSE
  )

  expect_equal(nrow(adapted), nrow(raw))
  expect_equal(adapted$injection_id, raw[["Analysis"]])
  expect_equal(adapted$analysis_no, raw[["Analysis"]])
  expect_equal(adapted$peak_number, raw[["Peak Nr"]])
  expect_equal(adapted$delta_value, raw[["d 13C/12C"]])
  expect_equal(adapted$area_or_voltage, raw[["Area 44"]])
  expect_equal(adapted$area_45, raw[["Area 45"]])
  expect_equal(adapted$amplitude_44, raw[["Ampl  44"]])
  expect_equal(adapted$bgd_46, raw[["BGD 46"]])
  expect_equal(adapted$amount, raw[["Identifier 2"]])
  expect_equal(adapted$amount_unit, raw[["Comment"]])
})

test_that("GB standard prefixes match configured base IDs without losing raw labels", {
  adapted <- adapt_gb_data(
    raw_df = gb_test_raw_data(),
    config = gb_test_config(),
    schema = gb_test_schema(),
    verbose = FALSE
  )

  expect_equal(
    adapted$sample_id,
    c("NaHCO3", "NaHCO3", "ETH4", "NaHCO3X_1umol", "on_off", "blank")
  )
  expect_equal(
    adapted$standard_id,
    c("NaHCO3", "NaHCO3", "ETH4", NA, NA, NA)
  )
  expect_equal(
    adapted$standard_variant[1:3],
    c("_1umol", "_5umol", "_5umol")
  )
  expect_equal(adapted$raw_sample_id, gb_test_raw_data()[["Identifier 1"]])
  expect_equal(adapted$is_standard, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_true(adapted$is_quality_control[[5]])
  expect_true(adapted$is_blank[[6]])
})

test_that("GB adapter reports missing mapped export columns", {
  raw <- gb_test_raw_data()
  raw[["Area 44"]] <- NULL

  expect_error(
    adapt_gb_data(
      raw_df = raw,
      config = gb_test_config(),
      schema = gb_test_schema(),
      verbose = FALSE
    ),
    "Area 44"
  )
})

test_that("GB adapter requires configured standard IDs", {
  config <- gb_test_config()
  config$stds_used <- character()

  expect_error(
    adapt_gb_data(
      raw_df = gb_test_raw_data(),
      config = config,
      schema = gb_test_schema(),
      verbose = FALSE
    ),
    "must list standards"
  )
})
test_that("process_gb applies configurable name exclusions and keeps GB peaks", {
  raw <- gb_test_raw_data()
  extra <- raw[c(1, 1, 1), , drop = FALSE]
  extra[["Analysis"]] <- 107:109
  extra[["Identifier 1"]] <- c("linearity", "CIT_5umol", "f2_si_26087")
  extra[["Peak Nr"]] <- c(1, 4, 5)
  raw <- rbind(raw, extra)

  config <- gb_test_config()
  config$processing$peak_number <- 3
  config$processing$sample_peak_number_high <- 5

  result <- process_gb(
    raw_df = raw,
    config = config,
    schema = gb_test_schema(),
    verbose = FALSE
  )

  expect_equal(
    sort(result$excluded_data$raw_sample_id),
    sort(c("on_off", "blank", "linearity", "CIT_5umol", "f2_si_26087"))
  )
  expect_equal(nrow(result$data), 3L)
  expect_true(all(result$data$peak_number == 4))
  expect_equal(result$excluded_peaks$peak_number, 11)
  expect_equal(nrow(result$all_peak_data), 4L)
  expect_false(any(result$data$peak_number == 3))
  expect_equal(nrow(result$output), nrow(result$data))
  expect_true(all(nzchar(result$excluded_data$exclusion_reason)))
  expect_false(result$drift_correction$applied)
  expect_equal(result$data$delta_value, result$data$delta_value_raw)

  custom <- process_gb(
    raw_df = raw,
    config = config,
    schema = gb_test_schema(),
    exclude_patterns = "CIT",
    sample_peak_numbers = 1:13,
    verbose = FALSE
  )
  expect_true("CIT_5umol" %in% custom$excluded_data$raw_sample_id)
  expect_true("on_off" %in% custom$data$raw_sample_id)
})

test_that("process_gb can exclude whole analyses by analysis number", {
  result <- process_gb(
    raw_df = gb_test_raw_data(),
    config = gb_test_config(),
    schema = gb_test_schema(),
    exclude_analysis_numbers = 101,
    verbose = FALSE
  )

  expect_false(any(result$data$analysis_no == 101))
  expect_equal(unique(result$excluded_analyses$analysis_no), 101)
  expect_equal(nrow(result$excluded_analyses), 1L)
  expect_true(all(result$excluded_analyses$analysis_exclusion_reason ==
    "analysis number listed for exclusion"))
})

test_that("process_gb corrects linear drift using one mean per ETH4 analysis", {
  raw <- gb_test_raw_data()
  raw[3, "d 13C/12C"] <- -9.4
  extra_standards <- raw[rep(3, 3), , drop = FALSE]
  extra_standards[["Analysis"]] <- c(107, 108, 109)
  extra_standards[["Identifier 1"]] <- "ETH4_5umol"
  extra_standards[["Peak Nr"]] <- 4
  extra_standards[["d 13C/12C"]] <- c(-9.0, -8.9, -8.8)
  raw <- rbind(raw, extra_standards)

  result <- process_gb(
    raw_df = raw,
    config = gb_test_config(),
    schema = gb_test_schema(),
    verbose = FALSE
  )

  expect_true(result$drift_correction$applied)
  expect_equal(result$drift_correction$n_standard_analyses, 4L)
  expect_equal(result$drift_correction$slope_per_analysis, 0.1, tolerance = 1e-10)
  expect_equal(result$drift_correction$r_squared, 1, tolerance = 1e-10)
  expect_true(all(result$data$drift_correction_applied))
  expect_true(all(is.finite(result$data$delta_value_raw)))
  expect_lt(sd(result$data$delta_value[result$data$standard_id == "ETH4"]), 1e-10)
  expect_gt(
    result$data$delta_value[result$data$analysis_no == 101],
    result$data$delta_value_raw[result$data$analysis_no == 101]
  )
})
