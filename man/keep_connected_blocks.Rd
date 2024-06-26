% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_filter.R
\name{keep_connected_blocks}
\alias{keep_connected_blocks}
\title{Keep fully connected blocks}
\usage{
keep_connected_blocks(object, block, verbose = TRUE)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{block}{svar}

\item{verbose}{TRUE or FALSE}
}
\description{
Keep fully connected blocks
}
\examples{
file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
object <- read_metabolon(file)
object \%<>\% keep_connected_blocks(  block = 'Subject')
}
