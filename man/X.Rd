% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_fit_limma.R
\name{X}
\alias{X}
\alias{beta}
\title{Model based prediction}
\usage{
X(
  object,
  formula = default_formula(object),
  drop = varlevels_dont_clash(object, all.vars(formula)),
  codingfun = contr.treatment.explicit
)

beta(object, fit = fits(object)[1])
}
\arguments{
\item{object}{SummarizedExperiment or data.frame}

\item{formula}{formula}

\item{drop}{TRUE or FALSE}

\item{codingfun}{function}

\item{fit}{'limma', 'lm', 'lme', 'wilcoxon'}
}
\value{
beta matrix (nlevel x nfeature)
}
\description{
Model based prediction
}
\examples{
file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
object <- read_metabolon(file)
object \%<>\% fit_limma(block = 'Subject')
beta(object)                    #    betas : nlevel x nfeature
   X(object)                    #   design : nlevel x nlevel
   X(object) \%*\% beta(object)   # response : nlevel x nfeature
}
