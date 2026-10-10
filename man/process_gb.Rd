\name{process_gb}
\alias{process_gb}
\title{Process GB-IRMS data with drift correction}
\usage{
process_gb(
  raw_df,
  config,
  schema = NULL,
  isotope = "C13",
  run_id = NULL,
  exclude_patterns = NULL,
  exclude_analysis_numbers = NULL,
  drift_standard_id = NULL,
  apply_drift_correction = TRUE,
  sample_peak_numbers = NULL,
  verbose = TRUE
)
}
\arguments{
\item{raw_df}{Raw GB-IRMS data frame.}
\item{config}{Experiment configuration returned by \code{load_IRMS_config()}.}
\item{schema}{Optional GB peripheral schema.}
\item{isotope}{Isotope system to process, such as \code{"C13"}.}
\item{run_id}{Optional run identifier. Defaults to the experiment name.}
\item{exclude_patterns}{Name patterns to exclude. Defaults to the GB patterns or config value.}
\item{exclude_analysis_numbers}{Analysis numbers to exclude. Defaults to \code{config$processing$exclude_analysis_numbers} or none.}
\item{drift_standard_id}{Standard used to estimate run drift. Defaults to \code{config$processing$drift_standard_id} or "ETH4".}
\item{apply_drift_correction}{Apply a linear analysis-number correction when at least three drift-standard analyses are available.}
\item{sample_peak_numbers}{Optional explicit peak-number vector overriding the configured low/high range.}
\item{verbose}{Print progress information.}
}
\value{
A list containing selected sample peak rows and output, all retained peak rows,
rows excluded by name or analysis number, peaks outside the configured range,
and drift model diagnostics. Corrected \code{delta_value} and raw
\code{delta_value_raw} are both retained.
}
\description{
Adapts raw GB peak rows, applies configurable exclusions, selects the configured
sample peak range, and estimates linear run drift from one mean per analysis of
the selected standard. The correction is centered at the mean fitted standard
value and applied to rows by instrument analysis number. Other peak rows remain
available for review.
}
