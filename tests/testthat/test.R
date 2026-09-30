test_that("package configuration files are available", {
  config_dir <- system.file(
    "config",
    package = "isotidy"
  )

  expect_true(dir.exists(config_dir))

  expect_true(
    file.exists(
      file.path(config_dir, "experiment_config_TEMPLATE.yaml")
    )
  )

  expect_true(
    file.exists(
      file.path(config_dir, "isotope_systems.yml")
    )
  )

  expect_true(
    file.exists(
      file.path(config_dir, "peripheral_schemas.yml")
    )
  )
})


test_that("peripheral schemas can be loaded", {
  schemas <- load_peripheral_schemas()

  expect_type(schemas, "list")
  expect_true(length(schemas) > 0)
})

test_that("a peripheral schema can be retrieved", {
  schemas <- load_peripheral_schemas()

  expect_true(length(schemas) > 0)

  peripheral_name <- names(schemas)[1]

  schema <- get_peripheral_schema(
    peripheral = peripheral_name,
    schemas = schemas
  )

  expect_type(schema, "list")
  expect_true(length(schema) > 0)
})


test_that("a peripheral schema passes validation", {
  schemas <- load_peripheral_schemas()

  peripheral_name <- names(schemas)[1]

  schema <- get_peripheral_schema(
    peripheral = peripheral_name,
    schemas = schemas
  )

  result <- validate_peripheral_schema(schema)

  expect_true(isTRUE(result))
})

test_that("EA data can be adapted to canonical format", {
  schemas <- load_peripheral_schemas()

  ea_schema <- get_peripheral_schema(
    peripheral = "EA",
    schemas = schemas
  )

  config <- list(
    experiment_name = "test_experiment",
    processing = list(
      element = "C"
    )
  )

  raw_ea <- data.frame(
    `Identifier 1` = c("Sample_A", "Sample_B"),
    Amount = c(1.2, 1.5),
    `Peak Nr` = c(3, 3),
    `d 13C/12C` = c(-25.3, -27.1),
    `Area All` = c(1000, 1200),
    `Area 44` = c(900, 1100),
    `Area 45` = c(50, 60),
    `Area 46` = c(20, 25),
    `Ampl  44` = c(10, 12),
    `Ampl  45` = c(1, 1.2),
    `Ampl  46` = c(0.5, 0.6),
    `BGD 44` = c(2, 2.5),
    `BGD 45` = c(0.2, 0.3),
    `BGD 46` = c(0.1, 0.1),
    `Time Code` = c(100, 200),
    check.names = FALSE
  )

  adapted <- adapt_ea_data(
    raw_df = raw_ea,
    config = config,
    schema = ea_schema,
    isotope = "C13",
    verbose = FALSE
  )

  expect_equal(nrow(adapted), 2)

  expect_equal(
    adapted$run_id,
    c("test_experiment", "test_experiment")
  )

  expect_equal(
    adapted$sample_id,
    c("Sample_A", "Sample_B")
  )

  expect_equal(
    adapted$element,
    c("C", "C")
  )

  expect_equal(
    adapted$delta_value,
    c(-25.3, -27.1)
  )

  expect_equal(
    adapted$area_or_voltage,
    c(1000, 1200)
  )

  expect_equal(
    adapted$peak_number,
    c(3, 3)
  )

  expect_true(all(!adapted$is_standard))
  expect_true(all(!adapted$is_blank))
})

test_that("EA adapter rejects empty raw data", {
  schemas <- load_peripheral_schemas()

  ea_schema <- get_peripheral_schema(
    peripheral = "EA",
    schemas = schemas
  )

  config <- list(
    experiment_name = "test_experiment",
    processing = list(
      element = "C"
    )
  )

  expect_error(
    adapt_ea_data(
      raw_df = data.frame(),
      config = config,
      schema = ea_schema,
      isotope = "C13",
      verbose = FALSE
    ),
    "raw_df.*contains no rows"
  )
})

test_that("EA adapter rejects missing mapped columns", {
  schemas <- load_peripheral_schemas()

  ea_schema <- get_peripheral_schema(
    peripheral = "EA",
    schemas = schemas
  )

  config <- list(
    experiment_name = "test_experiment",
    processing = list(
      element = "C"
    )
  )

  raw_ea <- data.frame(
    `Identifier 1` = "Sample_A",
    Amount = 1.2,
    `Peak Nr` = 3,
    `d 13C/12C` = -25.3,
    check.names = FALSE
  )

  expect_error(
    adapt_ea_data(
      raw_df = raw_ea,
      config = config,
      schema = ea_schema,
      isotope = "C13",
      verbose = FALSE
    ),
    "missing mapped columns"
  )
})

test_that("EA adapter rejects unsupported isotope", {
  schemas <- load_peripheral_schemas()

  ea_schema <- get_peripheral_schema(
    peripheral = "EA",
    schemas = schemas
  )

  config <- list(
    experiment_name = "test_experiment",
    processing = list(
      element = "C"
    )
  )

  raw_ea <- data.frame(
    `Identifier 1` = "Sample_A",
    Amount = 1.2,
    `Peak Nr` = 3,
    `d 13C/12C` = -25.3,
    `Area All` = 1000,
    check.names = FALSE
  )

  expect_error(
    adapt_ea_data(
      raw_df = raw_ea,
      config = config,
      schema = ea_schema,
      isotope = "S34",
      verbose = FALSE
    ),
    "No isotope ratio mapping found"
  )
})


