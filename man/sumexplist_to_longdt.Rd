% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.0_read_omics.R
\name{sumexplist_to_longdt}
\alias{sumexplist_to_longdt}
\title{SummarizedExperiment list to long data.table}
\usage{
sumexplist_to_longdt(
  sumexplist,
  svars = intersect("subgroup", autonomics::svars(sumexplist[[1]])),
  fvars = intersect("gene", autonomics::fvars(sumexplist[[1]])),
  setvarname = "set"
)
}
\arguments{
\item{sumexplist}{list of SummarizedExperiments}

\item{svars}{character vector}

\item{fvars}{character vector}

\item{setvarname}{string}
}
\value{
data.table
}
\description{
SummarizedExperiment list to long data.table
}
\examples{
# RNA
    rnafile <- download_data('billing19.rnacounts.txt')
    rna <- read_rnaseq_counts(rnafile)
    fdt(rna)$gene <- fdt(rna)$gene_name
# PRO/FOS
    profile <- download_data('billing19.proteingroups.txt')
    fosfile <- download_data('billing19.phosphosites.txt')
    subgroups <- paste0(c('E00', 'E01', 'E02', 'E05', 'E15', 'E30', 'M00'), '_STD')
    pro <- read_maxquant_proteingroups(file = profile, subgroups = subgroups)
    fos <- read_maxquant_phosphosites(
              phosphofile = fosfile, proteinfile = profile, subgroups = subgroups)
    pro$subgroup \%<>\% stringi::stri_replace_first_fixed('_STD', '')
    fos$subgroup \%<>\% stringi::stri_replace_first_fixed('_STD', '')
# sumexplist to longdt
    sumexplist <- list(rna = rna, pro = pro, fos = fos)
    dt <- sumexplist_to_longdt(sumexplist, setvarname = 'platform')
    dt \%<>\% extract(gene \%in\% c('TNMD', 'TSPAN6'))
}
