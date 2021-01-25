% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_read_rnaseq.R
\name{log2counts}
\alias{log2counts}
\alias{log2counts,SummarizedExperiment-method}
\alias{log2counts<-}
\alias{log2counts<-,SummarizedExperiment,matrix-method}
\alias{log2counts<-,SummarizedExperiment,numeric-method}
\title{Get/Set log2counts}
\usage{
log2counts(object)

\S4method{log2counts}{SummarizedExperiment}(object)

log2counts(object) <- value

\S4method{log2counts}{SummarizedExperiment,matrix}(object) <- value

\S4method{log2counts}{SummarizedExperiment,numeric}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment}

\item{value}{count matrix (features x samples)}
}
\value{
count matrix (get) or updated object (set)
}
\description{
Get / Set log2counts matrix
}
\examples{
file <- download_data('billing16.rnacounts.txt')
object <- read_rnaseq_counts(file)
log2counts(object) <- exprs(object)
log2counts(object)[1:3, 1:3]
}