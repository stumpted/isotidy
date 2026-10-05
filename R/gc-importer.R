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
