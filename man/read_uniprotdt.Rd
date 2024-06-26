% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3.2_read_annotations.R
\name{read_uniprotdt}
\alias{read_uniprotdt}
\alias{parse_maxquant_hdrs}
\alias{read_contaminantdt}
\title{Read fasta hdrs}
\usage{
read_uniprotdt(fastafile, fastafields = FASTAFIELDS, verbose = TRUE)

parse_maxquant_hdrs(fastahdrs)

read_contaminantdt(force = FALSE, verbose = TRUE)
}
\arguments{
\item{fastafile}{string (or charactervector)}

\item{fastafields}{charactervector : which fastahdr fields to extract ?}

\item{verbose}{bool}

\item{fastahdrs}{character vector}

\item{force}{whether to overwrite existing file}
}
\value{
data.table(uniprot, protein, gene, uniprot, reviewed, existence)
}
\description{
Read fasta hdrs
}
\note{
existence values are always those of the canonical isoform
      (no isoform-level resolution for this field)
}
\examples{
# uniprot hdrs
     fastafile <- system.file('extdata/uniprot_hsa_20140515.fasta', package = 'autonomics')
     read_uniprotdt(fastafile)
     
# maxquant hdrs
    file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
    dt <- .read_maxquant_proteingroups(file)
    parse_maxquant_hdrs(dt$`Fasta headers`)

    profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
    fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics' )
    prodt <- .read_maxquant_proteingroups(profile)
    fosdt <- .read_maxquant_phosphosites(fosfile, profile)
    parse_maxquant_hdrs(prodt$`Fasta headers`)
    parse_maxquant_hdrs(fosdt$`Fasta headers`)
    
# contaminant hdrs
    read_contaminantdt()

}
