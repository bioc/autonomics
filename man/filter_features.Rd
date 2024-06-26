% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_filter.R
\name{filter_features}
\alias{filter_features}
\title{Filter features on condition}
\usage{
filter_features(object, condition, verbose = TRUE)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{condition}{filter condition}

\item{verbose}{logical}
}
\value{
filtered eSet
}
\description{
Filter features on condition
}
\examples{
file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
object <- read_metabolon(file)
filter_features(object, SUPER_PATHWAY == 'Lipid')
}
