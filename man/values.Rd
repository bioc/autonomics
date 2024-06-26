% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_get_set.R
\name{values}
\alias{values}
\alias{values,SummarizedExperiment-method}
\alias{values<-}
\alias{values<-,SummarizedExperiment,matrix-method}
\alias{values<-,SummarizedExperiment,numeric-method}
\title{Get/Set expr values}
\usage{
values(object)

\S4method{values}{SummarizedExperiment}(object)

values(object) <- value

\S4method{values}{SummarizedExperiment,matrix}(object) <- value

\S4method{values}{SummarizedExperiment,numeric}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment}

\item{value}{ratio matrix (features x samples)}
}
\value{
value matrix (get) or updated object (set)
}
\description{
Get/Set value matrix
}
\examples{
file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
object <- read_maxquant_proteingroups(file)
values(object)[1:3, 1:3]
values(object) <- 0
values(object)[1:3, 1:3]
}
