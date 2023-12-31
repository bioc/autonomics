% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.1_read_rnaseq.R
\name{log2tpmratios}
\alias{log2tpmratios}
\alias{log2tpmratios,SummarizedExperiment-method}
\alias{log2tpmratios<-}
\alias{log2tpmratios<-,SummarizedExperiment,matrix-method}
\alias{log2tpmratios<-,SummarizedExperiment,numeric-method}
\title{Get/Set log2tpmratios}
\usage{
log2tpmratios(object)

\S4method{log2tpmratios}{SummarizedExperiment}(object)

log2tpmratios(object) <- value

\S4method{log2tpmratios}{SummarizedExperiment,matrix}(object) <- value

\S4method{log2tpmratios}{SummarizedExperiment,numeric}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment}

\item{value}{log2tpmratios matrix (features x samples)}
}
\value{
log2tpmratios matrix (get) or updated object (set)
}
\description{
Get / Set log2tpmratios matrix
}
\examples{
file <- download_data('billing19.rnacounts.txt')
object <- read_rnaseq_counts(file, plot=FALSE)
log2tpmratios(object) <- values(object)
log2tpmratios(object)[1:3, 1:3]
}
