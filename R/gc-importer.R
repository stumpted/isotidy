import_gc_data <- function(
    file_path,
    sheet = "New Table"
) {

  raw_df <- readxl::read_excel(
    file_path,
    sheet = sheet,
    col_names = FALSE
  )

  header_names <- apply(
    raw_df[1:3, ],
    2,
    function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & x != ""]

      if (length(x) == 0) {
        NA_character_
      } else {
        x[length(x)]
      }
    }
  )

  units <- as.character(raw_df[4, ])

  data <- raw_df[-c(1, 2, 3, 4), ]

  peak_number_col <- which(
    vapply(
      data,
      function(x) {
        any(
          grepl(
            "^Peak\\s+\\d+$",
            as.character(x)
          ),
          na.rm = TRUE
        )
      },
      logical(1)
    )
  )

  if (length(peak_number_col) != 1) {
    stop(
      "Could not uniquely identify the Peak Number column."
    )
  }

  header_names[peak_number_col] <- "Peak Number"

  missing_names <- is.na(header_names) | header_names == ""

  header_names[missing_names] <- paste0(
    "unnamed_",
    which(missing_names)
  )

  header_names <- make.unique(header_names)

  names(data) <- header_names
  names(units) <- header_names

  attr(data, "units") <- units

  metadata_columns <- header_names[
    seq_len(peak_number_col - 1)
  ]

  data <- data |>
    tidyr::fill(
      dplyr::all_of(metadata_columns),
      .direction = "down"
    )

  data
}


#' Import all GC Excel exports for an experiment
#'
#' Reads `gc$selected_input_set` from the experiment YAML, then uses that entry
#' in `gc$input_sets` to resolve the input directory, filename pattern, and
#' worksheet. Matching files are imported with [import_gc_data()] and combined
#' while exact duplicate observations across files are collapsed. The
#' `source_file` column records all files that contained each observation,
#' joined with `; `.
#'
#' @param config Experiment configuration returned by `load_IRMS_config()`.
#' @param pattern Optional regular-expression override for the selected input
#'   set's `file_pattern`.
#' @param sheet Optional worksheet override for the selected input set's
#'   `sheet`.
#' @param verbose Print a summary of imported files and rows.
#' @return A combined raw GC data frame with one row per distinct observation,
#'   a `source_file` column listing its source files, and the `units` attribute
#'   retained from the imported exports.
#' @export
load_gc_data <- function(
    config,
    pattern = NULL,
    sheet = NULL,
    verbose = TRUE
) {

  if (!inherits(config, "irms_config")) {
    stop("'config' must be returned by load_IRMS_config().")
  }

  yaml_file_path <- config$yaml_file_path
  if (
    is.null(yaml_file_path) || length(yaml_file_path) != 1 ||
    is.na(yaml_file_path) || !nzchar(yaml_file_path) ||
    !file.exists(yaml_file_path)
  ) {
    stop("A valid config$yaml_file_path is required to load GC inputs.")
  }

  yaml_config <- yaml::read_yaml(yaml_file_path)
  gc_config <- yaml_config$gc
  selected_input_set <- gc_config$selected_input_set
  input_sets <- gc_config$input_sets

  if (
    is.null(selected_input_set) || length(selected_input_set) != 1 ||
    is.na(selected_input_set) || !nzchar(as.character(selected_input_set))
  ) {
    stop("Experiment YAML must define gc$selected_input_set.")
  }

  selected_input_set <- as.character(selected_input_set)
  if (!is.list(input_sets) || is.null(names(input_sets))) {
    stop("Experiment YAML must define named gc$input_sets.")
  }

  input_set <- input_sets[[selected_input_set]]
  if (is.null(input_set)) {
    stop(
      "Selected GC input set '", selected_input_set,
      "' was not found in gc$input_sets."
    )
  }

  directory <- input_set$directory
  if (
    is.null(directory) || length(directory) != 1 ||
    is.na(directory) || !nzchar(directory)
  ) {
    stop("Selected GC input set must define a non-empty 'directory'.")
  }

  if (is.null(pattern)) {
    pattern <- input_set$file_pattern
  }
  if (is.null(sheet)) {
    sheet <- input_set$sheet
  }
  if (is.null(sheet)) {
    sheet <- "New Table"
  }
  if (length(sheet) != 1 || is.na(sheet) || !nzchar(sheet)) {
    stop("The GC worksheet name must be one non-empty value.")
  }

  config_dir <- dirname(
    normalizePath(yaml_file_path, winslash = "/", mustWork = TRUE)
  )
  is_absolute_directory <- grepl(
    "^(?:[A-Za-z]:|/|\\\\\\\\)",
    directory,
    perl = TRUE
  )
  input_dir <- if (is_absolute_directory) {
    directory
  } else {
    file.path(config_dir, directory)
  }
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)

  if (
    is.null(input_dir) || length(input_dir) != 1 ||
    is.na(input_dir) || !nzchar(input_dir)
  ) {
    stop("The experiment configuration has no valid raw_data_dir.")
  }

  is_absolute_path <- grepl(
    "^(?:[A-Za-z]:|/|\\\\\\\\)",
    input_dir,
    perl = TRUE
  )

  input_dir_candidates <- input_dir
  if (!is_absolute_path) {
    yaml_file_path <- config$yaml_file_path

    if (
      !is.null(yaml_file_path) &&
      length(yaml_file_path) == 1 &&
      !is.na(yaml_file_path) &&
      nzchar(yaml_file_path) &&
      file.exists(yaml_file_path)
    ) {
      yaml_dir <- dirname(
        normalizePath(yaml_file_path, winslash = "/", mustWork = TRUE)
      )

      repeat {
        input_dir_candidates <- c(
          input_dir_candidates,
          file.path(yaml_dir, input_dir)
        )

        parent_dir <- dirname(yaml_dir)
        if (identical(parent_dir, yaml_dir)) {
          break
        }
        yaml_dir <- parent_dir
      }
    }
  }

  input_dir_candidates <- unique(vapply(
    input_dir_candidates,
    normalizePath,
    character(1),
    winslash = "/",
    mustWork = FALSE
  ))

  existing_input_dirs <- input_dir_candidates[
    dir.exists(input_dir_candidates)
  ]

  if (length(existing_input_dirs) == 0) {
    stop(
      "GC input directory for set '", selected_input_set,
      "' not found. Checked: ",
      paste(input_dir_candidates, collapse = ", ")
    )
  }

  input_dir <- existing_input_dirs[[1]]

  if (
    length(pattern) != 1 || is.na(pattern) || !nzchar(pattern)
  ) {
    stop("'pattern' must be one non-empty regular expression.")
  }

  files <- list.files(
    path = input_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  files <- sort(files)

  if (length(files) == 0) {
    stop(
      "No GC input files matched pattern '", pattern,
      "' in: ", input_dir
    )
  }

  imported <- lapply(files, function(file_path) {
    data <- import_gc_data(
      file_path = file_path,
      sheet = sheet
    )

    if ("source_file" %in% names(data)) {
      stop(
        "The imported GC table already contains a 'source_file' column: ",
        basename(file_path)
      )
    }

    data$source_file <- basename(file_path)
    data
  })

  units <- lapply(imported, attr, which = "units")
  if (!all(vapply(units, identical, logical(1), units[[1]]))) {
    stop("GC input files have inconsistent column units.")
  }

  data <- dplyr::bind_rows(imported)
  attr(data, "units") <- units[[1]]

  # The same observation can appear in overlapping exports. Treat every
  # imported field other than the filename as part of its identity, so
  # differences in sample, injection, peak, or measurement are preserved.
  observation_columns <- setdiff(names(data), "source_file")
  data <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(observation_columns))) |>
    dplyr::summarise(
      source_file = paste(sort(unique(.data$source_file)), collapse = "; "),
      .groups = "drop"
    )
  attr(data, "units") <- units[[1]]

  if (verbose) {
    message(
      "Imported ", length(files), " GC file(s) for input set '",
      selected_input_set, "': ", nrow(data),
      " distinct observation(s) from ", length(files), " file(s) in ",
      input_dir
    )
  }

  data
}
