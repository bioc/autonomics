% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_fit_limma.R
\name{varlevels_dont_clash}
\alias{varlevels_dont_clash}
\alias{varlevels_dont_clash.data.table}
\alias{varlevels_dont_clash.SummarizedExperiment}
\title{Are varlevels unique}
\usage{
varlevels_dont_clash(object, ...)

\method{varlevels_dont_clash}{data.table}(object, vars = names(object), ...)

\method{varlevels_dont_clash}{SummarizedExperiment}(object, vars = svars(object), ...)
}
\arguments{
\item{object}{SummarizedExperiment or data.table}

\item{...}{required for s3 dispatch}

\item{vars}{character vector}
}
\value{
TRUE or FALSE
}
\description{
Are varlevels unique
}
\examples{
require(data.table)
object1 <- data.table(expand.grid(genome = c('WT', 'MUT'), treat = c('control', 'drug')))
object2 <- data.table(expand.grid(mutant = c('YES', 'NO'), treated = c('YES', 'NO')))
varlevels_dont_clash(object1)
varlevels_dont_clash(object2)
}
