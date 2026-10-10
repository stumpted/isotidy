\name{adapt_gb_data}
\alias{adapt_gb_data}
\title{Adapt GB-IRMS data}
\usage{
adapt_gb_data(
  raw_df,
  config = NULL,
  schema = NULL,
  isotope = "C13",
  run_id = NULL,
  instrument = "GB-IRMS",
  verbose = TRUE
)
}
\arguments{
\item{raw_df}{Raw GB-IRMS data frame with one row per peak.}
\item{config}{Experiment configuration returned by \code{load_IRMS_config()}.}
\item{schema}{Optional GB peripheral schema; loaded from package configuration by default.}
\item{isotope}{Isotope system to adapt, such as \code{"C13"}.}
\item{run_id}{Optional run identifier. Defaults to the experiment name.}
\item{instrument}{Instrument label stored in canonical data.}
\item{verbose}{Print progress information.}
}
\value{
A canonical data frame with one row per raw GB peak. The full raw sample label,
standard suffix, sample amount, and amount unit are retained.
}
\description{
Maps GB-IRMS column names into isotidy's canonical data format. Labels that
start with a configured standard ID followed by an underscore are matched to
that base standard ID while the original label is retained.
}