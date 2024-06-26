% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_transform.R
\name{subtract_baseline}
\alias{subtract_baseline}
\alias{subtract_pairs}
\alias{subtract_differences}
\title{Subtract baseline}
\usage{
subtract_baseline(
  object,
  subgroupvar,
  subgroupctr = slevels(object, subgroupvar)[1],
  block = NULL,
  assaynames = setdiff(assayNames(object), c("weights", "pepcounts")),
  verbose = TRUE
)

subtract_pairs(
  object,
  subgroupvar = "subgroup",
  subgroupctr = slevels(object, subgroupvar)[1],
  block,
  assaynames = assayNames(object)[1],
  verbose = TRUE
)

subtract_differences(object, block, subgroupvar, verbose = TRUE)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{subgroupvar}{subgroup svar}

\item{subgroupctr}{control subgroup}

\item{block}{block svar (within which subtraction is performed)}

\item{assaynames}{which assays to subtract for}

\item{verbose}{TRUE/FALSE}
}
\value{
SummarizedExperiment
}
\description{
Subtract baseline level within block
}
\details{
\code{subtract_baseline} subtracts baseline levels within block, using the 
medoid baseline sample if multiple exist. \cr

\code{subtract_pairs} also subtracts baseline level within block. 
It cannot handle multiple baseline samples, but has instead been optimized
for many blocks \cr

\code{subtract_differences} subtracts differences between subsequent levels, 
again within block
}
\examples{
# read 
    file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
    object0 <- read_metabolon(file)
    pca(object0, plot = TRUE, color = 'Time')

# subtract_baseline: takes medoid of baseline samples if multiple
    object <- subtract_baseline(object0, block = 'Subject', subgroupvar = 'Time')
    pca(object, plot = TRUE, color = 'Time')

# subtract_pairs: optimized for many blocks
    object <- subtract_pairs(object0, block = 'Subject', subgroupvar = 'Time')
    pca(object, plot = TRUE, color = 'Time')

# subtract_differences
    object <- subtract_differences(object0, block = 'Subject', subgroupvar = 'Time')
    values(object) \%<>\% na_to_zero()
    pca(object, plot = TRUE, color = 'Time')
}
