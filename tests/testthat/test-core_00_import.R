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
