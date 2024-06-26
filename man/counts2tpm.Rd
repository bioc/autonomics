% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.1_read_rnaseq.R
\name{counts2tpm}
\alias{counts2tpm}
\title{counts to tpm}
\usage{
counts2tpm(x, genesize)
}
\arguments{
\item{x}{count matrix}

\item{genesize}{genesize vector (kilobase)}
}
\value{
tpm matrix
}
\description{
counts to tpm
}
\examples{
file <- system.file('extdata/billing19.rnacounts.txt', package = 'autonomics')
object <- read_rnaseq_counts(file)
counts(object)[1:3, 1:3]
counts2tpm(counts(object), genesize = 1)[1:3, 1:3]
}
