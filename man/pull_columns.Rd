% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_utilities.R
\name{pull_columns}
\alias{pull_columns}
\title{Pull columns in a dataframe to the front}
\usage{
pull_columns(df, first_cols, verbose = TRUE)
}
\arguments{
\item{df}{data.frame}

\item{first_cols}{character vector: columns to be pulled to the front}

\item{verbose}{TRUE (default) or FALSE}
}
\value{
dataframe with re-ordered columns
}
\description{
Pull columns in a dataframe to the front
}
\examples{
df <- data.frame(
   symbol = c('A1BG', 'A2M'),
   id     = c('1',    '2'),
   name   = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'),
   type   = c('proteincoding', 'proteincoding'))
first_cols <- c('id', 'symbol', 'location', 'uniprot')
pull_columns(df, first_cols)
}
