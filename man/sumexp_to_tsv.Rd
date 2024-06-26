% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_sumexp_to_dt.R
\name{sumexp_to_tsv}
\alias{sumexp_to_tsv}
\title{Write sumexp to tsv}
\usage{
sumexp_to_tsv(object, assay = assayNames(object)[1], file)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{assay}{string}

\item{file}{filename}
}
\description{
Write sumexp to tsv
}
\examples{
file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
object <- read_maxquant_proteingroups(file, fit = 'limma')
tsv <- file.path(tempdir(), 'fukuda20.proteingroups.tsv')
sumexp_to_tsv(object, file = tsv)
}
