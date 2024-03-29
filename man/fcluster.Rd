% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_plot.R
\name{fcluster}
\alias{fcluster}
\title{Feature cluster}
\usage{
fcluster(
  object,
  method = c("pamk", "hclust", "apcluster"),
  k = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{verbose}{TRUE or FALSE}

\item{clusterer}{clustering engine}
}
\value{
SummarizedExperiment. Two new fvars: `cluster` and  `clustorder`
}
\description{
Feature cluster
}
\examples{
file <- download_data('atkin.metabolon.xlsx')
object <- read_metabolon(file)
fdt(object) \%<>\% extract(, 1)
if (require(fpc))        fdt(fcluster(object, method = 'pamk'         ))
if (require(stats))      fdt(fcluster(object, method = 'hclust', k = 3))
if (require(apcluster))  fdt(fcluster(object, method = 'apcluster'    ))
}
