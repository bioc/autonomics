% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.3_read_diann.R
\name{rm_diann_contaminants}
\alias{rm_diann_contaminants}
\title{Rm contaminants}
\usage{
rm_diann_contaminants(
  object,
  contaminants = read_contaminants(),
  verbose = TRUE
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{contaminants}{uniprots (character vector)}

\item{verbose}{TRUE or FALSE}
}
\value{
SummarizedExperiment
}
\description{
Rm contaminants from DIA-NN SumExp
}
\examples{
file <- download_data('dilution.report.tsv')
object <- read_diann_proteingroups(file)
object \%<>\% rm_diann_contaminants()
}
