% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_get_set.R
\name{fnames}
\alias{fnames}
\alias{fnames,SummarizedExperiment-method}
\alias{fnames<-}
\alias{fnames<-,SummarizedExperiment,character-method}
\title{Get/Set fnames}
\usage{
fnames(object)

\S4method{fnames}{SummarizedExperiment}(object)

fnames(object) <- value

\S4method{fnames}{SummarizedExperiment,character}(object) <- value
}
\arguments{
\item{object}{SummarizedExperiment, eSet, or EList}

\item{value}{character vector with feature names}
}
\value{
feature name vector (get) or updated object (set)
}
\description{
Get/Set feature names
}
\examples{
file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
object <- read_maxquant_proteingroups(file)
fnames(object) \%<>\% paste0('protein_', .)
object
}
