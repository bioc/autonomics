% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.1_read_rnaseq.R
\name{log2countsratios}
\alias{log2countsratios}
\alias{log2countsratios,SummarizedExperiment-method}
\alias{log2countsratios<-}
\alias{log2countsratios<-,SummarizedExperiment,matrix-method}
\alias{log2countsratios<-,SummarizedExperiment,numeric-method}
\title{Get/Set log2countsratios}
\usage{
log2countsratios(object)

\S4method{log2countsratios}{SummarizedExperiment}(object)

log2countsratios(object) <- value

\S4method{log2countsratios}{SummarizedExperiment,matrix}(object) <- value

\S4method{log2countsratios}{SummarizedExperiment,numeric}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment}

\item{value}{log2countsratios matrix (features x samples)}
}
\value{
log2countsratios matrix (get) or updated object (set)
}
\description{
Get / Set log2countsratios matrix
}
\examples{
file <- download_data('billing19.rnacounts.txt')
object <- read_rnaseq_counts(file, plot=FALSE)
log2countsratios(object) <- values(object)
log2countsratios(object)[1:3, 1:3]
}
