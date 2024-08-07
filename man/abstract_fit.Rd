% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_enrichment.R
\name{abstract_fit}
\alias{abstract_fit}
\title{Abstract model fit}
\usage{
abstract_fit(
  object,
  sep = guess_fitsep(fdt(object)),
  fit = fits(object),
  coef = coefs(object, fit = fit),
  significancevar = "p",
  significance = 0.05
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{sep}{string}

\item{fit}{character vector}

\item{coef}{character vector}

\item{significancevar}{'p' or 'fdr'}

\item{significance}{fraction : pvalue cutoff}
}
\value{
SummarizedExperiment
}
\description{
Abstract model fit
}
\examples{
file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
object <- read_metabolon(file, fit = 'limma', coef = 't3-t0')
fdt(object)
fdt(abstract_fit(object))
}
