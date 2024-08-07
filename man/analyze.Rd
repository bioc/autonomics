% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_analyze.R
\name{analyze}
\alias{analyze}
\title{Analyze}
\usage{
analyze(
  object,
  pca = TRUE,
  pls = TRUE,
  fit = "limma",
  formula = default_formula(object),
  drop = varlevels_dont_clash(object, all.vars(formula)),
  codingfun = contr.treatment.explicit,
  contrasts = NULL,
  coefs = colnames(create_design(object, formula = formula, drop = drop)),
  block = NULL,
  weightvar = if ("weights" \%in\% assayNames(object)) "weights" else NULL,
  plot = pca & !is.null(fit),
  label = "feature_id",
  palette = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{pca}{TRUE / FALSE: perform pca ?}

\item{pls}{TRUE / FALSE: perform pls ?}

\item{fit}{linmod engine: 'limma', 'lm', 'lme(r)', 'lmer', 'wilcoxon'}

\item{formula}{model formula}

\item{drop}{TRUE / FALSE : drop varname in designmat ?}

\item{codingfun}{factor coding function
\itemize{
    \item contr.treatment:          intercept = y0,     coefi = yi - y0
    \item contr.treatment.explicit: intercept = y0,     coefi = yi - y0
    \item code_control:             intercept = ymean,  coefi = yi - y0
    \item contr.diff:               intercept = y0,     coefi = yi - y(i-1)
    \item code_diff:                intercept = ymean,  coefi = yi - y(i-1)
    \item code_diff_forward:        intercept = ymean,  coefi = yi - y(i+)
    \item code_deviation:           intercept = ymean,  coefi = yi - ymean (drop last)
    \item code_deviation_first:     intercept = ymean,  coefi = yi - ymean (drop first)
    \item code_helmert:             intercept = ymean,  coefi = yi - mean(y0:(yi-1))
    \item code_helmert_forward:     intercept = ymean,  coefi = yi - mean(y(i+1):yp)
}}

\item{contrasts}{model coefficient contrasts of interest: string vector or NULL}

\item{coefs}{model coefficients          of interest: string vector or NULL}

\item{block}{model blockvar}

\item{weightvar}{NULL or name of weight matrix in assays(object)}

\item{plot}{TRUE / FALSE}

\item{label}{fvar}

\item{palette}{NULL or colorvector}

\item{verbose}{TRUE / FALSE: message?}
}
\value{
SummarizedExperiment
}
\description{
Analyze
}
\examples{
file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
object <- read_metabolon(file)
object \%<>\% analyze()
}
