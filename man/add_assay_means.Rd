% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_plot_volcano.R
\name{add_assay_means}
\alias{add_assay_means}
\title{Add assay means}
\usage{
add_assay_means(object, assay = assayNames(object)[1], bin = TRUE)
}
\arguments{
\item{object}{SummarizedExperiment or NULL}

\item{assay}{string}

\item{bin}{TRUE or FALSE}
}
\value{
SummarizedExperiment
}
\description{
Add assay means
}
\examples{
file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
object <- read_maxquant_proteingroups(file)
fdt(object) \%<>\% extract(, 1:2)
fdt(object)
object \%<>\% add_assay_means(SummarizedExperiment::assayNames(.))
fdt(object)
}
