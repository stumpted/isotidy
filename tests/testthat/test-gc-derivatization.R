test_that("AA and NACME carbon counts are retained from standards YAML", {
  records <- list(
    AA_mix = list(
      element = "C",
      standard_id = "AA",
      compound_values = list(
        Ala = list(
          delta_value = -19.57429981449309,
          delta_uncertainty = 0.1496319220104934,
          AA_C_count = 3,
          NACME_C_count = 3
        )
      ),
      applicable_peripherals = "GC"
    )
  )

  standards <- isotidy:::normalize_yaml_standards(records)

  expect_equal(standards$AA_C_count, 3)
  expect_equal(standards$NACME_C_count, 3)
})

test_that("NACME estimates and sample correction match synthetic workbook", {
  benchmark_path <- testthat::test_path(
    "..", "fixtures", "DeltaCalculations_example.xlsx"
  )
  benchmark_std <- readxl::read_excel(
    benchmark_path,
    sheet = "stdInfo"
  ) |>
    dplyr::filter(as.character(.data$Date) == "260512",
                  .data$compound == "Ala")
  benchmark_sample <- readxl::read_excel(
    benchmark_path,
    sheet = "derivCorr",
    col_names = FALSE
  )[[2]][[4]]
  benchmark_external <- readxl::read_excel(
    benchmark_path,
    sheet = "esCorr",
    col_names = FALSE
  )

  standard_reference <- data.frame(
    standard_id = "AA",
    compound_id = "Ala",
    delta_value = -19.57429981449309,
    delta_uncertainty = 0.1496319220104934,
    AA_C_count = 3,
    NACME_C_count = 3
  )
  date <- as.Date("2026-05-12")
  canonical <- data.frame(
    injection_id = c("STD_1", "STD_2", "SAMPLE_1", "SAMPLE_2", "F8_1"),
    sample_id = c("STD", "STD", "WT_rep3", "WT_rep3", "F8"),
    start_time = c(date, date, date, as.Date("2026-05-13"), date),
    compound_id = "Ala",
    delta_value_external_corrected = c(
      as.numeric(benchmark_external[[2]][[4]]),
      as.numeric(benchmark_external[[3]][[4]]),
      as.numeric(benchmark_external[[4]][[4]]),
      as.numeric(benchmark_external[[4]][[4]]),
      -30
    ),
    is_external_standard = c(FALSE, FALSE, FALSE, FALSE, TRUE)
  )

  estimates <- isotidy:::estimate_gc_nacme(
    canonical_df = canonical,
    standards_reference_df = standard_reference,
    standard_id = "AA",
    sample_id = "STD"
  )

  expect_equal(nrow(estimates), 1)
  expect_equal(estimates$n_standard_injections, 2L)
  expect_equal(estimates$delta_derivatized_mean,
               benchmark_std$`δderivatized`, tolerance = 1e-12)
  expect_equal(estimates$delta_NACME, benchmark_std$`δreagent`,
               tolerance = 1e-10)
  expect_equal(estimates$sigma_NACME, benchmark_std$`σreagent`,
               tolerance = 1e-10)

  corrected <- isotidy:::correct_gc_nacme(
    canonical_df = canonical,
    nacme_estimates = estimates,
    sample_id = "STD"
  )

  expect_true(corrected$is_derivatization_standard[[1]])
  expect_true(is.na(corrected$delta_value_nacme_corrected[[1]]))
  expect_equal(corrected$delta_value_final[[3]], benchmark_sample,
               tolerance = 1e-10)
  expect_true(is.na(corrected$delta_value_final[[4]]))
  expect_equal(corrected$delta_value_final[[5]], -30)
})

test_that("standards without AA or NACME counts remain valid", {
  records <- list(
    F8 = list(
      element = "C",
      standard_id = "F8",
      compound_values = list(C14M = list(delta_value = -29.98)),
      applicable_peripherals = "GC"
    )
  )

  standards <- isotidy:::normalize_yaml_standards(records)

  expect_true(is.na(standards$AA_C_count[[1]]))
  expect_true(is.na(standards$NACME_C_count[[1]]))
})

test_that("NACME estimation cannot use the external standard as its reference", {
  canonical <- data.frame(
    injection_id = "F8_1",
    sample_id = "F8",
    start_time = as.Date("2026-05-12"),
    compound_id = "C14M",
    delta_value_external_corrected = -30,
    is_external_standard = TRUE
  )
  standard_reference <- data.frame(
    standard_id = "F8",
    compound_id = "C14M",
    delta_value = -29.98,
    delta_uncertainty = NA_real_,
    AA_C_count = NA_real_,
    NACME_C_count = NA_real_
  )

  expect_error(
    isotidy:::estimate_gc_nacme(
      canonical_df = canonical,
      standards_reference_df = standard_reference,
      standard_id = "F8",
      sample_id = "F8",
      external_standard_reference_id = "F8"
    ),
    "also configured as the external standard"
  )
})
