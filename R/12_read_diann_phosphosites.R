
#' @rdname read_diann_phosphodiffs
#' @export
read_diann_pgmatrix <- function(dir){
    file <- file.path(dir, 'report.pg_matrix.tsv')
    assert_all_are_existing_files(file)
    object <- fread(file)
    cols <- c('Protein.Group', 'Protein.Names', 'Genes', 'First.Protein.Description')
    fdt0 <- object[, cols, with = FALSE]
    fdt0 %<>% pull_columns('Protein.Names')
    setnames(fdt0, 'Protein.Names',             'feature_id')
    setnames(fdt0, 'Protein.Group',             'uniprot'   )
    setnames(fdt0, 'Genes',                     'gene'      )
    setnames(fdt0, 'First.Protein.Description', 'description')
    object[, (cols) := NULL ]
    object <- cbind(fdt0[, .(feature_id)], object)
    object %<>% dt2mat()
    object %<>% matrix2sumexp()
    assayNames(object) <- 'pgmatrix'
    object %<>% merge_fdt(fdt0)
    assays(object)$pgmatrix %<>% zero_to_na()
    assays(object)$pgmatrix %<>%  log2()
    assayNames(object) %<>% paste0('log2.', .)
    object
    
}


#' @rdname read_diann_phosphodiffs
#' @export
read_diann_phosphosites <- function(dir){
    file <- file.path(dir, 'report.phosphosites_90.tsv')
    assert_all_are_existing_files(file)
    object <- fread(file)
    cols <- c('Protein', 'Protein.Names', 'Gene.Names', 'Residue', 'Site', 'Sequence')
    fdt0 <- object[, cols, with = FALSE]
    Protein.Names <- Residue <- Site <- NULL
    fdt0[, feature_id := paste0(Protein.Names, '.', Residue, Site) ]
    fdt0 %<>% pull_columns('feature_id')
    setnames(fdt0, 'Protein',       'uniprot' )
    setnames(fdt0, 'Protein.Names', 'protein' )
    setnames(fdt0, 'Gene.Names',    'gene'    )
    setnames(fdt0, 'Residue',       'residue' )
    setnames(fdt0, 'Site',          'site'    )
    setnames(fdt0, 'Sequence',      'sequence')
    object[, (cols) := NULL ]
    object <- cbind(fdt0[, .(feature_id)], object)
    object %<>% dt2mat()
    object %<>% matrix2sumexp()
    assayNames(object) <- 'phosphosites90'
    object %<>% merge_fdt(fdt0)
    assays(object)$phosphosites90 %<>% zero_to_na()
    assays(object)$phosphosites90 %<>%  log2()
    assayNames(object) %<>% paste0('log2.', .)
    object
}


#' Read diann phosphosites
#' @param dir  directory with 'report_pgmatrix' and 'report.phosphosites_90.tsv'
#' @return SummarizedExperiment
#' @export
read_diann_phosphodiffs <- function(dir){
    
# Read proteingroups
    proteingroups <- read_diann_pgmatrix(dir)
    phosphosites  <- read_diann_phosphosites(dir)

# Common samples
    commonsamples <- intersect(snames(proteingroups), snames(phosphosites  ))
    proteinsamples <-  setdiff(snames(proteingroups), snames(phosphosites  ))
    phosphosamples <-  setdiff(snames(phosphosites),  snames(proteingroups ))
    if (length(phosphosamples)>0)  cmessage('Drop %d/%d phospho samples with no protein counterpart', length(phosphosamples), ncol(phosphosites ))
    proteingroups %<>% extract(, commonsamples)
    phosphosites  %<>% extract(, commonsamples)
    
# Common features
    commonfeatures <- intersect(fdt(phosphosites)$uniprot, fdt(proteingroups)$uniprot)
    proteinfeatures <- setdiff(fdt(proteingroups)$uniprot, fdt(phosphosites )$uniprot)
    phosphofeatures <- setdiff(fdt(phosphosites )$uniprot, fdt(proteingroups)$uniprot)
    if (length(phosphofeatures)>0)  cmessage('Retain %d/%d phospho features with protein counterpart', nrow(phosphosites ) - length(phosphofeatures), nrow(phosphosites ))
    uniprot <- NULL
    proteingroups %<>% filter_features(uniprot %in% commonfeatures, verbose = FALSE)
    phosphosites  %<>% filter_features(uniprot %in% commonfeatures, verbose = FALSE)
    assert_has_no_duplicates(proteingroups$uniprot)
    idx <- match(fdt(phosphosites )$uniprot, fdt(proteingroups)$uniprot)
    proteingroups %<>% extract(idx, )
    fnames(proteingroups) <- fdt(proteingroups)$feature_id <- fnames(phosphosites)
    assays(phosphosites)$log2.pgmatrix <- assays(proteingroups)$log2.pgmatrix
    
# Differences
    assays(phosphosites)$log2.phosphodiffs90 <- 
    assays(phosphosites)$log2.phosphosites90 - 
    assays(phosphosites)$log2.pgmatrix
    assays(phosphosites) %<>% extract(c('log2.phosphodiffs90', 'log2.phosphosites90', 'log2.pgmatrix'))
    phosphosites
}

