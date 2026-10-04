#' Save processed EA results and a human-readable HTML processing report
#'
#' Saves the final processed EA data as a CSV and renders the permanent
#' EA processing-summary R Markdown template as an HTML report.
#'
#' @param result Result returned by process_ea().
#' @param config Experiment configuration used for processing.
#' @param output_dir Directory where the CSV and HTML report will be saved.
#' @param prefix Optional filename prefix. Defaults to the experiment name.
#' @param verbose Print progress messages.
#'
#' @return A list containing paths to the data and HTML report.
#' @export
save_ea_results <- function(
    result,
    config,
    output_dir,
    prefix = NULL,
    verbose = TRUE
) {

  if (!is.list(result)) {
    stop("'result' must be the list returned by process_ea().")
  }

  if (is.null(result$output)) {
    stop(
      "The process_ea() result does not contain 'output'. ",
      "Make sure process_ea() has completed the final output step."
    )
  }

  if (!is.list(config)) {
    stop("'config' must be the experiment configuration used by process_ea().")
  }

  if (
    length(output_dir) != 1 ||
    is.na(output_dir) ||
    !nzchar(output_dir)
  ) {
    stop("'output_dir' must contain one non-empty path.")
  }

  if (is.null(prefix)) {

    prefix <- config$experiment_name

    if (
      is.null(prefix) ||
      is.na(prefix) ||
      !nzchar(prefix)
    ) {
      stop(
        "No filename prefix was supplied and ",
        "config$experiment_name is missing."
      )
    }
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  data_file <- file.path(
    output_dir,
    paste0(prefix, "_processed.csv")
  )

  report_file <- file.path(
    output_dir,
    paste0(prefix, "_processing_summary.html")
  )

  if (verbose) {
    cat("Saving processed EA results\n")
  }

  utils::write.csv(
    result$output,
    file = data_file,
    row.names = FALSE
  )

  if (verbose) {
    cat(
      "  Data:   ",
      data_file,
      "\n",
      sep = ""
    )
  }

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "The 'rmarkdown' package is required to create the HTML report. ",
      "Install it with install.packages('rmarkdown')."
    )
  }

  # Find the installed package template
  template <- system.file(
    "reports",
    "EA_processing_summary.Rmd",
    package = "IRMS_data_pipeline"
  )

  # Fallback for development with the package source tree
  if (!nzchar(template)) {

    development_template <- file.path(
      "inst",
      "reports",
      "EA_processing_summary.Rmd"
    )

    if (file.exists(development_template)) {
      template <- development_template
    }
  }

  if (
    !nzchar(template) ||
    !file.exists(template)
  ) {
    stop(
      "Could not find the EA processing report template. ",
      "Expected: inst/reports/EA_processing_summary.Rmd"
    )
  }

  if (verbose) {
    cat(
      "  Report: ",
      report_file,
      "\n",
      sep = ""
    )
  }

  rmarkdown::render(
    input = template,
    output_file = basename(report_file),
    output_dir = output_dir,
    params = list(
      result = result,
      config = config,
      data_file = data_file
    ),
    envir = new.env(parent = globalenv()),
    quiet = !verbose
  )

  if (verbose) {
    cat("  ✓ EA results saved\n\n")
  }

  invisible(
    list(
      data_file = normalizePath(
        data_file,
        winslash = "/",
        mustWork = FALSE
      ),
      report_file = normalizePath(
        report_file,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}
