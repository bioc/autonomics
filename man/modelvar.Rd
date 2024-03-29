% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8_fit.R
\name{modelvar}
\alias{modelvar}
\alias{effectvar}
\alias{tvar}
\alias{pvar}
\alias{fdrvar}
\alias{abstractvar}
\alias{modelvec}
\alias{effectvec}
\alias{tvec}
\alias{pvec}
\alias{fdrvec}
\alias{modeldt}
\alias{effectdt}
\alias{tdt}
\alias{pdt}
\alias{modelmat}
\alias{effectmat}
\alias{effectsizemat}
\alias{tmat}
\alias{pmat}
\alias{fdrmat}
\alias{modelfeatures}
\alias{upfeatures}
\alias{downfeatures}
\title{Get model variable}
\usage{
modelvar(
  featuredt,
  quantity,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

effectvar(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

tvar(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

pvar(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

fdrvar(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

abstractvar(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

modelvec(
  featuredt,
  quantity,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id"
)

effectvec(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id"
)

tvec(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id"
)

pvec(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id"
)

fdrvec(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id"
)

modeldt(
  featuredt,
  quantity,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

effectdt(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

tdt(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

pdt(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

modelmat(
  featuredt,
  quantity,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

effectmat(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

effectsizemat(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

tmat(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

pmat(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

fdrmat(
  featuredt,
  fit = fits(featuredt),
  coef = default_coefs(featuredt, fit = fit)
)

modelfeatures(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id",
  significancevar = "p",
  significance = 0.05,
  effectdirection = "<>",
  effectsize = 0
)

upfeatures(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id",
  significancevar = "p",
  significance = 0.05,
  effectsize = 0
)

downfeatures(
  featuredt,
  fit = fits(featuredt)[1],
  coef = default_coefs(featuredt, fit = fit)[1],
  fvar = "feature_id",
  significancevar = "p",
  significance = 0.05,
  effectsize = 0
)
}
\arguments{
\item{quantity}{'p', 'effect', 'fdr', 't', or 'se'}

\item{fit}{string (vector)}

\item{coef}{string (vector)}

\item{fvar}{'feature_id' or other fvar for values (pvec) or names (upfeatures)}

\item{significancevar}{'p' or 'fdr'}

\item{significance}{p or fdr cutoff (fractional number)}

\item{effectdirection}{'<>', '<' or '>'}

\item{effectsize}{effectsize cutoff (positive number)}

\item{object}{SummarizedExperiment}
}
\value{
string (tvar), matrix (tmat), numeric vector (tvec), character vector (tfeatures)
}
\description{
Get model variable
}
\examples{
# Read
    file <- download_data('atkin.metabolon.xlsx')
    object <- read_metabolon(file)
    object \%<>\% fit_limma()
    object \%<>\% fit_lm()
    featuredt <- fdt(object)
# modelvar
    modelvar(featuredt, 'p');                                          pvar(featuredt)
    modelvar(featuredt, 'effect');                                effectvar(featuredt)
    modelvar(featuredt, 'fdr');                                      fdrvar(featuredt)
# modelvec
    modelvec(featuredt, 'p'     )[1:3];                                 pvec(featuredt)[1:3]
    modelvec(featuredt, 'effect')[1:3];                            effectvec(featuredt)[1:3]
    modelvec(featuredt, 'fdr'   )[1:3];                               fdrvec(featuredt)[1:3]
# modelmatrix
    modelmat(featuredt, 'p'     )[1:3, 1:3];                            pmat(featuredt)[1:3, 1:3]
    modelmat(featuredt, 'effect')[1:3, 1:3];                       effectmat(featuredt)[1:3, 1:3]
    modelmat(featuredt, 'fdr'   )[1:3, 1:3];                          fdrmat(featuredt)[1:3, 1:3]
# modelfeatures
    modelfeatures(featuredt      )[1:3]
    modelfeatures(featuredt, effectdirection = '<' )[1:3];      downfeatures(featuredt)[1:3]
    modelfeatures(featuredt, effectdirection = '>' )[1:3];        upfeatures(featuredt)[1:3]
}
