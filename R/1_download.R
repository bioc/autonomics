
#' Data used in examples/vignette/tests/longtests
#' @examples 
#' AUTONOMICS_DATASETS
#' @export
AUTONOMICS_DATASETS <- c('atkin.somascan.adat',    # new versions
                         'atkin.metabolon.xlsx',
                         'atkin18.somascan.adat',  # old versions - backward compatibility
                         'atkin18.metabolon.xlsx',
                         'billing16.bam.zip',
                         'billing16.rnacounts.txt',
                         'billing16.proteingroups.txt',
                         'billing16.somascan.adat',
                         'billing19.rnacounts.txt',
                         'billing19.proteingroups.txt',
                         'billing19.phosphosites.txt',
                         'contaminants.tsv',
                         'dilution.report.tsv',
                         'fukuda20.proteingroups.txt',
                         'halama18.metabolon.xlsx',
                         'integer64.proteinGroups.txt',
                         'multiorganism.combined_protein.tsv',
                         'ddglucose.proteingroups.txt',
                         'uniprot_hsa_20140515.fasta')

#' @rdname download_data
#' @export
DATADIR <- file.path(R_user_dir('autonomics', 'cache'), 'datasets')

#' Download autonomics example data
#' 
#' @param filename  file name
#' \tabular{rlr}{
#'         'atkin.somascan.adat'      \tab   \href{https://pubmed.ncbi.nlm.nih.gov/30525282/}{Halama, 2018}   \tab  effects of hypoglycemia     \cr
#'         'atkin.metabolon.xlsx'     \tab                                                                    \tab                              \cr
#'     'billing16.bam.zip'            \tab   \href{https://pubmed.ncbi.nlm.nih.gov/26857143/}{Billing, 2016}  \tab  stemcell comparison         \cr
#'     'billing16.rnacounts.txt'      \tab                                                                    \tab                              \cr
#'     'billing16.somascan.adat'      \tab                                                                    \tab                              \cr
#'     'billing16.proteingroups.txt'  \tab                                                                    \tab                              \cr
#'     'billing19.rnacounts.txt'      \tab   \href{https://pubmed.ncbi.nlm.nih.gov/31332097/}{Billing, 2016}  \tab  stemcell differentiation    \cr
#'     'billing19.proteingroups.txt'  \tab                                                                    \tab                              \cr
#'     'billing19.phosphosites.txt'   \tab                                                                    \tab                              \cr
#'     'ddglucose.proteingroups.txt'  \tab   Omics Module                                                     \tab  glycolysis inhibitor        \cr
#'      'fukuda20.proteingroups.txt'  \tab   \href{https://pubmed.ncbi.nlm.nih.gov/32648304/}{Fukuda, 2020}   \tab  zebrafish development       \cr 
#'      'halama18.metabolon.xlsx'     \tab   \href{https://pubmed.ncbi.nlm.nih.gov/29777783/}{Halama, 2018}   \tab  glutaminase inhibitor
#' }
#' @param localdir local dir to save file to
#' @param verbose  TRUE / FALSE
#' @param force    TRUE / FALSE
#' @return local file path
#' @examples
#' # Show available datasets
#'     download_data()
#'     
#' # atkin 2018 - hypoglycemia - pubmed 30525282
#'     # download_data('atkin.somascan.adat')            # somascan  intensities
#'     # download_data('atkin.metabolon.xlsx')           # metabolon intensities
#'
#' # billing16 - stemcell characterization - pubmed 26857143
#'     # download_data('billing16.proteingroups.txt')      # proteingroup ratios
#'     # download_data('billing16.somascan.adat')          # somascan  intensities
#'     # download_data('billing16.rnacounts.txt')          # rnaseq    counts
#'     # download_data('billing16.bam.zip')                # rnaseq    alignments
#'
#' # billing19 - stemcell differentiation - pubmed 31332097
#'     # download_data('billing19.proteingroups.txt')    # proteingroup ratios
#'     # download_data('billing19.phosphosites.txt')     # phosphosite  ratios
#'     # download_data('billing19.rnacounts.txt')        # rnaseq       counts
#'
#' # fukuda20 - heart regeneration - pubmed PXD016235
#'     # download_data('fukuda20.proteingroups.txt')       # proteingroup LFQ
#'
#' # halama18 - glutaminase inhibition - pubmed 30525282
#'     # download_data('halama18.metabolon.xlsx')          # metabolon intensities
#' @export
download_data <- function(
    filename = NULL,
    localdir = file.path(DATADIR, split_extract_fixed(filename, '.', 1)),
    verbose  = TRUE, 
    force    = FALSE
){
# Assert/Prepare
    . <- NULL
    if (is.null(filename))  return(AUTONOMICS_DATASETS)
    assert_is_subset(filename, AUTONOMICS_DATASETS)
    dir.create(localdir, recursive = TRUE, showWarnings = FALSE)
    filepath <- file.path(localdir, filename)
# Download
    if (force | !file.exists(filepath)) {
        if(verbose) message( "Downloading ", filename)
        url <- "https://bitbucket.org/graumannlabtools/autonomics/downloads"
        url %<>% paste0('/', filename)
        for (i in 1:10){  # retry to deal with possible server timeout (due to coinciding requests)
            try({    Sys.sleep(abs(rnorm(1, sd = 2)))           # wait a few seconds
                     download.file(url, filepath, mode = 'wb')  # prevent corrupt xlsx ('wb')
                     break   })                                 # dont repeat if successful
        }
    }
# Unzip
    if (file_ext(filename)=='zip'){
        if (verbose)  message('unzip')
        unzip(filepath, exdir = substr(filepath, 1, nchar(filepath)-4))
        filepath %<>% substr(1, nchar(.)-4)
    }
# Return
    filepath
}


#' Download mcclain21 data
#' @param counts_or_samples 'counts' or 'samples'
#' @param localdir           dirname
#' @param force              TRUE or FALSE
#' @details \href{https://pubmed.ncbi.nlm.nih.gov/33597532/}{ Mc clain 2021: COVID19 transcriptomics:} 
#' @examples
#' download_mcclain21('counts')
#' download_mcclain21('samples')
#' @export
download_mcclain21 <- function(                  # Integrating into `download_data` not succesful
    counts_or_samples = 'counts',                #    - getGEOSuppFiles <-> download.file
    localdir = file.path(DATADIR, 'mcclain21'),  #    - GSE161731 not intuitive as filename
    force = FALSE                                #    - two files : counts + sampledt
){
    if (!dir.exists(localdir) | force){
        dir.create(localdir)
        GEOquery::getGEOSuppFiles("GSE161731", baseDir = localdir, makeDirectory = FALSE)
    }
    file  <- file.path(localdir, 'GSE161731_counts.csv.gz')
    sfile <- file.path(localdir, 'GSE161731_counts_key.csv.gz')
    out <- switch(counts_or_samples, counts = file, samples = sfile)
    return(out)
}
