% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_plot.R
\name{plot_exprs_per_coef}
\alias{plot_exprs_per_coef}
\title{Plot exprs per coef}
\usage{
plot_exprs_per_coef(
  object,
  fit = fits(object)[1],
  coefs = default_coefs(object, fit = fit),
  x = default_x(object),
  block = NULL,
  geom = default_geom(object, x, block = block),
  orderbyp = FALSE,
  title = x,
  subtitle = default_subtitle(fit, x, coefs),
  n = 1,
  nrow = 1,
  ncol = NULL,
  theme = ggplot2::theme(legend.position = "bottom", legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
)
}
\arguments{
\item{object}{SummarizedExperiment}

\item{fit}{'limma', 'lm', 'lme', 'lmer', 'wilcoxon'}

\item{coefs}{subset of coefs(object) to consider in selecting top}

\item{x}{x svar}

\item{block}{group svar}

\item{geom}{'boxplot' or 'point'}

\item{orderbyp}{TRUE or FALSE}

\item{title}{string}

\item{subtitle}{string}

\item{n}{number}

\item{nrow}{number of rows in faceted plot}

\item{ncol}{number of cols in faceted plot}

\item{theme}{ggplot2::theme(...) or NULL}
}
\value{
ggplot object
}
\description{
Plot exprs per coef
}
\examples{
file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
object <- read_metabolon(file)
object \%<>\% fit_limma()
object \%<>\% pls(by = 'subgroup')
object \%<>\% pls(by = 'Diabetes')
object \%<>\% pls(by = 'Subject')
plot_exprs_per_coef(object)
plot_exprs_per_coef(object, orderbyp = TRUE)
plot_exprs_per_coef(object, fit = 'pls1', block = 'Subject')
}
\seealso{
\code{\link{plot_sample_densities}},
         \code{\link{plot_sample_violins}}
}
