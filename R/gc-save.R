# Saving GC results: corrected peak data, correction diagnostics, and HTML report.

#' Save processed GC results and a processing summary
#'
#' Writes final corrected GC peak data, optional replicate-averaged
#' visualization data, excluded measurements, injection-level and daily
#' external-standard offsets, and an HTML processing summary.
#'
#' @param result Result returned by process_gc().
#' @param config Experiment configuration used for processing.
#' @param output_dir Directory for CSV files and the HTML report.
#' @param prefix Optional filename prefix. Defaults to config$file_name.
#' @param verbose Print progress messages.
#'
#' @return A list of paths to the output files, including the optional
#'   visualization CSV.
#' @export
save_gc_results <- function(
    result,
    config,
    output_dir,
    prefix = NULL,
    verbose = TRUE
) {
  if (!is.list(result) || is.null(result$output)) {
    stop("'result' must be the list returned by process_gc().")
  }

  if (!is.list(config)) {
    stop("'config' must be the experiment configuration used by process_gc().")
  }

  if (
    length(output_dir) != 1 ||
    is.na(output_dir) ||
    !nzchar(output_dir)
  ) {
    stop("'output_dir' must contain one non-empty path.")
  }

  if (is.null(prefix)) {
    prefix <- config$file_name
  }

  if (
    length(prefix) != 1 ||
    is.na(prefix) ||
    !nzchar(prefix)
  ) {
    stop("Supply a non-empty 'prefix' or config$file_name.")
  }

  prefix <- gsub("[^A-Za-z0-9._-]+", "_", as.character(prefix))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  data_file <- file.path(output_dir, paste0(prefix, "_gc_processed.csv"))
  visualization_file <- NULL
  if (!is.null(result$visualization_data)) {
    visualization_file <- file.path(
      output_dir,
      paste0(prefix, "_gc_visualization.csv")
    )
  }
  injection_offsets_file <- file.path(
    output_dir,
    paste0(prefix, "_gc_injection_offsets.csv")
  )
  daily_offsets_file <- file.path(
    output_dir,
    paste0(prefix, "_gc_daily_offsets.csv")
  )
  excluded_measurements_file <- NULL
  if (!is.null(result$excluded_measurements)) {
    excluded_measurements_file <- file.path(
      output_dir,
      paste0(prefix, "_gc_excluded_measurements.csv")
    )
  }
  report_file <- file.path(
    output_dir,
    paste0(prefix, "_gc_processing_summary.html")
  )

  utils::write.csv(result$output, data_file, row.names = FALSE)
  if (!is.null(visualization_file)) {
    utils::write.csv(result$visualization_data, visualization_file, row.names = FALSE)
  }
  utils::write.csv(result$injection_offsets, injection_offsets_file, row.names = FALSE)
  utils::write.csv(result$daily_offsets, daily_offsets_file, row.names = FALSE)
  if (!is.null(excluded_measurements_file)) {
    utils::write.csv(
      result$excluded_measurements,
      excluded_measurements_file,
      row.names = FALSE
    )
  }

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "The 'rmarkdown' package is required to create the HTML report. ",
      "Install it with install.packages('rmarkdown')."
    )
  }

  template <- system.file(
    "reports",
    "GC_processing_summary.Rmd",
    package = "isotidy"
  )

  if (!nzchar(template)) {
    development_template <- file.path(
      "inst",
      "reports",
      "GC_processing_summary.Rmd"
    )
    if (file.exists(development_template)) {
      template <- development_template
    }
  }

  if (!nzchar(template) || !file.exists(template)) {
    stop(
      "Could not find the GC processing report template. ",
      "Expected: inst/reports/GC_processing_summary.Rmd"
    )
  }

  rmarkdown::render(
    input = template,
    output_file = basename(report_file),
    output_dir = output_dir,
    params = list(
      result = result,
      config = config,
      data_file = data_file,
      injection_offsets_file = injection_offsets_file,
      daily_offsets_file = daily_offsets_file,
      excluded_measurements_file = excluded_measurements_file
    ),
    envir = new.env(parent = globalenv()),
    quiet = !verbose
  )

  if (verbose) {
    message("Saved GC results to ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
  }

  invisible(list(
    data_file = normalizePath(data_file, winslash = "/", mustWork = FALSE),
    visualization_file = if (is.null(visualization_file)) {
      NULL
    } else {
      normalizePath(visualization_file, winslash = "/", mustWork = FALSE)
    },
    injection_offsets_file = normalizePath(
      injection_offsets_file,
      winslash = "/",
      mustWork = FALSE
    ),
    daily_offsets_file = normalizePath(
      daily_offsets_file,
      winslash = "/",
      mustWork = FALSE
    ),
    excluded_measurements_file = if (is.null(excluded_measurements_file)) {
      NULL
    } else {
      normalizePath(
        excluded_measurements_file,
        winslash = "/",
        mustWork = FALSE
      )
    },
    report_file = normalizePath(report_file, winslash = "/", mustWork = FALSE)
  ))
}
