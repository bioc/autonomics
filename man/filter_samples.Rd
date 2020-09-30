% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter_samples}
\alias{filter_samples}
\title{Filter samples on condition}
\usage{
filter_samples(object, condition, verbose = FALSE, record = TRUE)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{condition}{filter condition}

\item{verbose}{TRUE or FALSE (default)}

\item{record}{TRUE (default) or FALSE}
}
\value{
filtered SummarizedExperiment
}
\description{
Filter samples on condition
}
\examples{
# GLUTAMINASE
    file <- download_data('glutaminase.metabolon.xlsx')
    object <- read_metabolon(file)
    filter_samples(object,   TIME_POINT=='h10',  verbose = TRUE)
}