% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.1_read_rnaseq.R
\name{log2cpmratios}
\alias{log2cpmratios}
\alias{log2cpmratios,SummarizedExperiment-method}
\alias{log2cpmratios<-}
\alias{log2cpmratios<-,SummarizedExperiment,matrix-method}
\alias{log2cpmratios<-,SummarizedExperiment,numeric-method}
\title{Get/Set log2cpmratios}
\usage{
log2cpmratios(object)

\S4method{log2cpmratios}{SummarizedExperiment}(object)

log2cpmratios(object) <- value

\S4method{log2cpmratios}{SummarizedExperiment,matrix}(object) <- value

\S4method{log2cpmratios}{SummarizedExperiment,numeric}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment}

\item{value}{log2cpmratios matrix (features x samples)}
}
\value{
log2cpmratios matrix (get) or updated object (set)
}
\description{
Get / Set log2cpmratios matrix
}
\examples{
file <- download_data('billing19.rnacounts.txt')
object <- read_rnaseq_counts(file, plot=FALSE)
log2cpmratios(object) <- values(object)
log2cpmratios(object)[1:3, 1:3]
}
