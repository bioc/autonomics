
#' proteingroup to isoforms
#' @param x proteingroups string vector
#' @param unique whether to remove duplicates
#' @return string vector
#' @examples 
#'  (x <- c('Q96JP5;Q96JP5-2', 'Q96JP5', 'Q96JP5-2;P86791'))
#'  pg_to_isoforms(x)
#'  pg_to_canonical(x)
#'  pg_to_isoforms( x, unique = FALSE)
#'  pg_to_canonical(x, unique = FALSE)
#' # .pg_to_isoforms(x[1])   # unexported dot functions
#' # .pg_to_canonical(x[1])  # operate on scalars
#' @export
pg_to_canonical <- function(x, unique = TRUE){
    assert_is_character(x)
    unname(vapply(x, .pg_to_canonical, character(1), unique = unique))
}

.pg_to_canonical <- function(x, unique = TRUE){
    z <- unlist(stri_split_fixed(x, ';'))
    z %<>% split_extract_fixed('-', 1)
    if ({{unique}}) z %<>% unique()
    paste0(z, collapse = ';')
}

#' @rdname pg_to_canonical
#' @export
pg_to_isoforms <- function(x, unique = TRUE){
    assert_is_character(x)
    unname(vapply(x, .pg_to_isoforms, character(1), unique = unique))
}

.pg_to_isoforms <- function(x, unique = TRUE){
    z <- unlist(stri_split_fixed(x, ';'))
    z %<>% split_extract_fixed('-', 2)
    z[z=='NA'] <- '0'
        # Sometimes the canonical isoform is NOT isoform-1 !
        # https://www.uniprot.org/uniprotkb/Q9H0P0/entry#sequences
    if ({{unique}}) z %<>% unique()
    paste0(z, collapse = ',')
}

#' Extract common substring
#' @param a first string
#' @param b second string
#' @return  string
#' @examples
#' a <- "heart-specific Fatty acid binding protein"
#' b <- "Fatty acid binding protein, isoform 3"
#' extract_common_substr(a, b)
#'
#' a <- "Small nuclear ribonucleoprotein-associated proteins B and B'"
#' b <- "Small nuclear ribonucleoprotein-associated protein N"
#' extract_common_substr(a, b)
#' @references https://stackoverflow.com/questions/28261825
#' @noRd
extract_common_substr <- function(a, b){

    tt <- drop(attr(adist(a, b, counts=TRUE), "trafos"))

    # Nothing in common
    if (!stri_detect_regex(tt, 'M+')) return('')

    # Something in common
    aa  <-  stri_sub(tt, stri_locate_all_regex(tt, '[DM]+')[[1]]) %>%
            paste0(collapse = '') %>% trimws()
            # paste is required because multiple substrings can be found
    #bb <- tt %>%
    # stri_sub(stri_locate_all_regex(tt, '[IM]+')[[1]]) %>%
    # paste0(collapse = '')

    stri_sub(a, stri_locate_all_regex(aa, 'M+')[[1]]) %>%
    paste0(collapse = '') %>%
    trimws()

    # different  = c(a %>%
    #    stri_sub(stri_locate_all_regex(aa, 'D+')[[1]]) %>%
    #    trimws(),
    # b %>%
    #    stri_sub(stri_locate_all_regex(bb, 'I+')[[1]])) %>%
    #    trimws())

}

#' Commonify strings
#' @param x character vector
#' @examples
#' # NO DIFFERENCES
#'    x <- c( 'Retrotransposon Gag-like protein 8B',
#'            'Retrotransposon Gag-like protein 8B')
#'    commonify_strings(x)
#' # TAILS DIFFER
#'    x <- c( 'Histone H2B type 1-K',
#'            'Histone H2B type 1-C/E/F/G/I')
#'    commonify_strings(x)
#'    x <- c("Small nuclear ribonucleoprotein-associated proteins B and B'",
#'           "Small nuclear ribonucleoprotein-associated protein N")
#'    commonify_strings(x)
#' # MORE COMPLEX DIFFERENCES
#'    x <- c( 'Fatty acid binding protein, isoform 3',
#'            'Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein, isoform 3')
#'    commonify_strings(x)
#' # NOTHING IN COMMON
#'    x <- c('ABC1', 'DEF2')
#'    commonify_strings(x)
#' @noRd
commonify_strings <- function(x){
    . <- NULL
    common <- Reduce(extract_common_substr, x)
    alternate  <- if (common==''){  x
                } else {            stri_replace_first_fixed(x, common, '') %>%
                                    stri_replace_first_fixed(', ', '')      %>%
                                    trimws()
                }
    if (all(alternate == '')) return(common)

    alternate                          %>%
    unique()                           %>%
    (function(s){s[s==''] <- '.'; s})  %>%
    sort()                             %>%
    #magrittr::extract(.!='')          %>%
    paste0(collapse='|')             %>%
    paste0('(', ., ')')              %>%
    paste0(common, .)
}


#' diann precursor quantity
#' @export
PRECURSOR_QUANTITY <- 'Precursor.Quantity'

excelcols <- function(file, range) colnames(read_excel(file, range = range, n_max = 0))
cols <- function(file)  names(fread(file, nrows = 0))
col1 <- function(file)  cols(file)[1]
col2 <- function(file)  cols(file)[2]
col3 <- function(file)  cols(file)[3]


# x <- 'Q15149;Q15149-3;Q15149-4;Q15149-8'
# uniprot2isoforms(x)
uniprot2isoforms <- function(x){
    x %<>% stri_split_fixed(';') 
    x %<>% unlist() 
    x %<>% split_extract_fixed('-', 2)
    x[x=='NA'] <- '0'
    #x %<>% unique()  #doesnt work with the way diann organises its proteingroups
    x %<>% sort()
    paste0(x, collapse = ',')
}

#' @rdname read_diann_proteingroups
#' @importFrom arrow read_parquet
#' @export
.read_diann_precursors <- function(
    file,
    Global.Q             = 0.01, 
    Q                    = 0.01,
    Global.PG.Q          = 0.01,
    PG.Q                 = 0.05,
    Global.Peptidoform.Q = 0.01,
    Peptidoform.Q        = 0.01,
    Lib.Q                = 0.01,
    Lib.PG.Q             = 0.01,
    Lib.Peptidoform.Q    = 0.01, 
    format               = c("tsv", "parquet")[1],
    verbose              = TRUE)
{
# Assert
    assert_is_fraction(Global.Q)
    assert_is_fraction(Q)
    assert_is_fraction(Global.PG.Q)
    assert_is_fraction(PG.Q)
    assert_is_fraction(Global.Peptidoform.Q)
    assert_is_fraction(Peptidoform.Q)
    assert_is_fraction(Lib.Q)
    assert_is_fraction(Lib.PG.Q)
    assert_is_fraction(Lib.Peptidoform.Q)
    assert_is_a_string(format)
    switch(format,
      'tsv'     = assert_diann_report(file),
      'parquet' = assert_diann_parquet_report(file),
      stop("Not implemented DIA-NN output format: ", format))
    assert_is_a_bool(verbose)
    iprecursor <- isoform    <- NULL
    log2maxlfq <- maxlfq     <- organism     <- pepcounts <- NULL
    precounts  <- precursor  <- preintensity <- protein   <- NULL
    run        <- top1       <- top3         <- total     <- NULL
    uniprot    <- NULL
# Read
    if (format == 'tsv')
    {
      anncols <- c('Run', 'Genes', 'Protein.Names', 'Protein.Group',
                 'Precursor.Id', 'Q.Value', 'Lib.PG.Q.Value',
                 'Stripped.Sequence')
      numcols <- c('Precursor.Quantity', 'PG.Quantity', 'PG.MaxLFQ')
    } else if (format == 'parquet')
    {
      anncols <- c('Run', 'Genes', 'Protein.Names', 'Protein.Group',
                 'Precursor.Id', 'Global.Q.Value', 'Q.Value', 
                 'Global.PG.Q.Value', 'PG.Q.Value',
                 'Global.Peptidoform.Q.Value', 'Peptidoform.Q.Value',
                 'Lib.Q.Value', 'Lib.PG.Q.Value', 'Lib.Peptidoform.Q.Value', 
                 'Stripped.Sequence')
      numcols <- c('Precursor.Quantity', 'PG.TopN', 'PG.MaxLFQ')
    } else {
      stop("Not implemented DIA-NN output format: ", format)
    }
    cols <- c(anncols, numcols)
    if (format == 'tsv')
    {
      dt <- fread(file, select = cols)                  # 1977.16 but 1,35E+11
      for (col in numcols){ dt[, (col) := stri_replace_first_fixed(get(col), ',', '.') ] 
        dt[, (col) := as.numeric(get(col))  ] }
    } else if (format == 'parquet') {
      dt <- read_parquet(file, col_select = cols) %>%
        as.data.table()
    } else stop("Not implemented DIA-NN output format: ", format)
    setnames(dt, 'Run',                'run')
    setnames(dt, 'Genes',              'gene')
    setnames(dt, 'Protein.Names',      'protein')
    setnames(dt, 'Protein.Group',      'uniprot')
    setnames(dt, 'Precursor.Id',       'precursor')
    if (format == 'parquet')
    {
      setnames(dt, 'Global.Q.Value',             'Global.Q')
      setnames(dt, 'Q.Value',                    'Q')
      setnames(dt, 'Global.PG.Q.Value',          'Global.PG.Q')
      setnames(dt, 'PG.Q.Value',                 'PG.Q')
      setnames(dt, 'Global.Peptidoform.Q.Value', 'Global.Peptidoform.Q')
      setnames(dt, 'Peptidoform.Q.Value',        'Peptidoform.Q')
      setnames(dt, 'Lib.Q.Value',                'Lib.Q')
      # setnames(dt, 'Lib.PG.Q.Value',             'Lib.PG.Q')
      setnames(dt, 'Lib.Peptidoform.Q.Value',    'Lib.Peptidoform.Q')
    }
    setnames(dt, 'Lib.PG.Q.Value',     'Lib.PG.Q')
    setnames(dt, 'Stripped.Sequence',  'sequence')
    setnames(dt, 'PG.MaxLFQ',          'maxlfq')
    setnames(
      dt,
      switch(format,
        'tsv'     = 'PG.Quantity',
        'parquet' = 'PG.TopN',
        stop("Not implemented DIA-NN output format: ", format)),
      'intensity')
    setnames(dt, 'Precursor.Quantity', 'preintensity')
# Filter
    if (format == 'parquet') dt %<>% .filter_dianne_proteingroups(
      Global.Q, Q, Global.PG.Q, PG.Q, Global.Peptidoform.Q, Peptidoform.Q,
      Lib.Q, Lib.Peptidoform.Q, verbose = verbose)
    dt %<>% .filter_dianne_proteingroups(Lib.PG.Q)
# Order precursors
    dt <- dt[, .SD[rev(order(preintensity))], by = c('uniprot', 'run')]
    dt[, iprecursor := seq_len(.N),                      by = c('uniprot', 'run')]
    dt[, precounts  := length(unique(precursor)),        by = c('uniprot', 'run')]
    dt[, pepcounts  := length(unique(sequence)),         by = c('uniprot', 'run')]
# Order proteingroups
    pgdt <- dt[, .(uniprot, run, pepcounts, precounts, maxlfq)]
    pgdt %<>% unique()
    pgdt %<>% extract(, .(pepcounts = sum(pepcounts), 
                          precounts = sum(precounts), 
                         log2maxlfq = log2(sum(maxlfq, na.rm = TRUE))), by = 'uniprot')
    pgdt %<>% extract(order(-pepcounts, -precounts, -log2maxlfq))
    dt[, uniprot := factor(uniprot, pgdt$uniprot)]
    dt %<>% extract(order(uniprot, run, iprecursor))
    dt[, uniprot := as.character(uniprot)]
# Intuify protein
    pgdt <- unique(dt[, .(uniprot, protein)])
    pgdt[protein=='', protein := uniprot]
    pgdt %<>% uncollapse(protein, sep = ';')                                     #     uncollapse
    pgdt[, organism := split_extract_fixed(protein, '_', 2)]                     #     drop organism
    pgdt[, protein  := split_extract_fixed(protein, '_', 1)]                     # 
    pgdt[, protein := commonify_strings(protein), by = c('uniprot', 'organism')] #     commonify  (within proteingroup/organism)
    pgdt %<>% recollapse(by = c('uniprot', 'organism'), sep = ';')               #     recollapse (within proteingroup/organism)
    pgdt[, protein := paste0(protein, '_', organism)]                            #     add organism
    pgdt %<>% recollapse(by = 'uniprot')                                         #     recollapse (within proteingroup)
# Add feature_id
    pgdt[, isoform := uniprot2isoforms(uniprot), by = 'uniprot']
    pgdt[, feature_id := paste0(protein, '-', isoform)]
    assert_has_no_duplicates(pgdt$feature_id)
    pgdt[, c('isoform') := NULL]
    # pgdt[, feature_name := forge_pg_descriptions(uniprot, protein, fastadt)]     #     add feature_name
    dt %<>% .merge(pgdt, by = 'uniprot')
    dt %<>% pull_columns(c('gene', 'protein', 'organism', 'feature_id', 'uniprot', 
                'run', 'pepcounts', 'precounts', 'iprecursor', 'precursor', 'sequence'))
# Summarize
    dt[, top1  :=     rev(sort(preintensity))[1],                  by = c('uniprot', 'run')]
    dt[, top3  := sum(rev(sort(preintensity))[1:3], na.rm = TRUE), by = c('uniprot', 'run')]
    dt[, total := sum(         preintensity,        na.rm = TRUE), by = c('uniprot', 'run')]
    dt[]
}

#' @importFrom rlang dots_list
.filter_dianne_proteingroups <- function(dt, ..., verbose = TRUE)
{
  filters <- dots_list(...,  .named = TRUE)
  assert_is_subset(c(names(filters), 'uniprot'), colnames(dt))
  for (fl in names(filters))
  {
    n0 <- length(unique(dt$uniprot))
    dt %<>% extract(dt[[fl]] < filters[[fl]])
    n1 <- length(unique(dt$uniprot))
    if (verbose)  message(
      '\t\tRetain ', n1, '/', n0, ' proteingroups: ', fl, ' < ', filters[[fl]])
  }
  dt
}

#' @rdname read_diann_proteingroups
#' @export
.read_diann_proteingroups <- function(
    file,
    format               = c("tsv", "parquet")[1],
    Global.Q             = 0.01, 
    Q                    = 0.01,
    Global.PG.Q          = 0.01,
    PG.Q                 = 0.05,
    Global.Peptidoform.Q = 0.01,
    Peptidoform.Q        = 0.01,
    Lib.Q                = 0.01,
    Lib.PG.Q             = 0.01,
    Lib.Peptidoform.Q    = 0.01,
    verbose              = TRUE
){
    dt <- .read_diann_precursors(
      file,
      format               = format,
      Global.Q             = Global.Q, 
      Q                    = Q,
      Global.PG.Q          = Global.PG.Q,
      PG.Q                 = PG.Q,
      Global.Peptidoform.Q = Global.Peptidoform.Q,
      Peptidoform.Q        = Peptidoform.Q,
      Lib.Q                = Lib.Q,
      Lib.PG.Q             = Lib.PG.Q,
      Lib.Peptidoform.Q    = Lib.Peptidoform.Q,
      verbose              = verbose)
    dt[, sequence := sequence[1], by = c('uniprot', 'run')]
    cols <- c('gene', 'feature_id', 'protein', 'organism', 'uniprot', 'run',
              'pepcounts', 'precounts', 'sequence',
              'intensity', 'top1', 'top3', 'total', 'maxlfq', 
              'Lib.PG.Q')
    dt %<>% extract(, cols, with = FALSE )
    dt %<>% unique()
    assert_is_identical_to_true(all(dt[, .N, by = c('run', 'feature_id')]$N==1))  # single row per run/protein - yes!
    dt
}

#' Read diann
#'
#' @param file                    DIA-NN report file
#' @param format                  Format of the report ('tsv' DIA-NN < v.2.0; 'parquet')
#' @param Q                       Q cutoff
#' @param Lib.Q                   Lib.Q cutoff
#' @param Global.Q                Global.Q cutoff
#' @param Lib.PG.Q                Lib.PG.Q cutoff
#' @param Global.PG.Q             Global.PG.Q cutoff
#' @param Lib.Peptidoform.Q       Lib.Peptidoform.Q cutoff
#' @param Global.Peptidoform.Q    Global.Peptidoform.Q cutoff
#' @param Peptidoform.Q           Peptidoform.Q cutoff
#' @param PG.Q                    PG.Q cutoff
#' @param simplify_snames         TRUE or FALSE: simplify (drop common parts in) samplenames ?
#' @param rm_contaminants         TRUE or FALSE: rm contaminants ?
#' @param impute                  TRUE or FALSE: impute group-specific NA values ?
#' @param plot                    TRUE or FALSE
#' @param pca                     TRUE or FALSE: run pca ?
#' @param pls                     TRUE or FALSE: run pls ?
#' @param fit                     model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param formula                 model formula
#' @param block                   model blockvar: string or NULL
#' @param coefs                   model coefficients    of interest: character vector or NULL
#' @param contrasts               coefficient contrasts of interest: character vector or NULL
#' @param palette                 color palette: named string vector
#' @param verbose                 TRUE or FALSE
#' @param ...                     used to maintain deprecated functions
#' @return  data.table or SummarizedExperiment
#' @details
#' Defaults for various Q value cutoffs corresppond to recommendations by the
#' DIA-NN teen for DIA-NN v.2 (as of 03.2025). Of these, the reader of the
#' legacy file format (flat tab seperated values, pre-DIA-NN v.2) only utilizes
#' Lib.PG.Q.
#' @examples
#' # Read
#'    file <- download_data('dilution.report.tsv')
#'    .read_diann_precursors(file)         #    precursors longdt
#'    .read_diann_proteingroups(file)      # proteingroups longdt
#'    fdt(read_diann_proteingroups(file))  # proteingroups sumexp
#' # Compare
#'     PR <- .read_diann_precursors(file)
#'     PG <- .read_diann_proteingroups(file)
#'     PG[intensity==top1] # matches      : 24975 (85%) proteingroups
#'     PG[intensity!=top1] # doesnt match :  4531 (15%) proteingroups
#'     RUN <- 'IPT_HeLa_1_DIAstd_Slot1-40_1_9997'
#'     PR[uniprot=='Q96JP5;Q96JP5-2' & run == RUN, 1:6] #    match:    8884 ==   8884
#'     PR[uniprot=='P36578'          & run == RUN, 1:6] # no match:  650887 != 407978
#'     PR[intensity != top1][feature_id == unique(feature_id)[1]][run == unique(run)[1]][1:2, 1:6]
#'     PR[intensity != top1][feature_id == unique(feature_id)[2]][run == unique(run)[1]][1:2, 1:6]
#'     PR[intensity != top1][feature_id == unique(feature_id)[3]][run == unique(run)[1]][1:3, 1:6]
#' @export
read_diann_proteingroups <- function(
                    file,
                  format = .guess_diann_format(file),
                Global.Q = 0.01, 
                       Q = 0.01,
             Global.PG.Q = 0.01,
                    PG.Q = 0.05,
    Global.Peptidoform.Q = 0.01,
           Peptidoform.Q = 0.01,
                   Lib.Q = 0.01,
                Lib.PG.Q = 0.01,
       Lib.Peptidoform.Q = 0.01,
         simplify_snames = TRUE,
         rm_contaminants = TRUE, 
                  impute = FALSE, 
                    plot = FALSE, 
                     pca = plot, 
                     pls = plot, 
                     fit = if (plot) 'limma' else NULL,
                 formula = as.formula('~ subgroup'),
                   block = NULL,
                   coefs = NULL,
               contrasts = NULL,
                 palette = NULL,
                 verbose = TRUE
){
# SumExp
    dt <- .read_diann_proteingroups(
      file,
      format               = format,
      Global.Q             = Global.Q, 
      Q                    = Q,
      Global.PG.Q          = Global.PG.Q,
      PG.Q                 = PG.Q,
      Global.Peptidoform.Q = Global.Peptidoform.Q,
      Peptidoform.Q        = Peptidoform.Q,
      Lib.Q                = Lib.Q,
      Lib.PG.Q             = Lib.PG.Q,
      Lib.Peptidoform.Q    = Lib.Peptidoform.Q,
      verbose              = verbose)
    assert_is_identical_to_true(length(unique(dt$run)) > 1) # SumExp generation fails on single run case
    object <- SummarizedExperiment(list(
        log2maxlfq    = dcast_diann(dt, 'maxlfq',    fill = NA, log2 = TRUE),
        log2intensity = dcast_diann(dt, 'intensity', fill = NA, log2 = TRUE),
        log2top1      = dcast_diann(dt, 'top1',      fill = NA, log2 = TRUE),
        log2top3      = dcast_diann(dt, 'top3',      fill = NA, log2 = TRUE),
        log2total     = dcast_diann(dt, 'total',     fill = NA, log2 = TRUE),
        pepcounts     = dcast_diann(dt, 'pepcounts', fill = 0              ),
        precounts     = dcast_diann(dt, 'precounts', fill = 0              ),
        sequence      = dcast_diann(dt, 'sequence',  fill = '')))
    sdt(object)$sample_id  <- snames(object)
    fdt(object)$feature_id <- fnames(object)
    analysis(object)$nfeatures <- nrow(object)
# fdt
    cols <- c('maxlfq', 'intensity', 'top1', 'top3', 'total', 
              'sequence', 'run', 'pepcounts', 'precounts')
    dt[, (cols) := NULL]
    dt %<>% unique()
    assert_are_identical(nrow(dt), nrow(object)) # if not more fields need to be NULLed
    object %<>% merge_fdt(dt)
    for (assay in assayNames(object))  object %<>% add_assay_means(assay)
# sdt
    snames(object) <- colnames(object)
    if (simplify_snames)  snames(object) %<>% simplify_snames()
    object %<>% add_subgroup()
# Filter. Impute. Analyze
    if (rm_contaminants)  object %<>% rm_diann_contaminants(verbose = verbose)
    object %<>% rm_missing_in_all_samples(verbose = verbose)
    object %<>% extract(order(rowVars(values(.), na.rm = TRUE)), )
  # object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose)
    if ({{impute}})   object %<>% impute()  # above breaks when all subgroups are singletons
    object %<>% analyze( pca         = pca,           pls         = pls,
                         fit         = fit,           formula     = formula,
                         block       = block,         coefs       = coefs,
                         contrasts   = contrasts,     verbose     = verbose,
                         plot        = plot,          palette     = palette )
    object
}

.guess_diann_format <- function (x)
{
  assert_is_a_string(x)
  assert_all_are_existing_files(x)
  
  parquet_magic <- charToRaw("PAR1")
  con_bin <- file(x, "rb")
  on.exit(close(con_bin), add = TRUE)
  header <- readBin(con_bin, what = "raw", n = 4)
  
  if (identical(header, parquet_magic)) {
    seek(con_bin, where = -4, origin = "end")
    footer <- readBin(con_bin, what = "raw", n = 4)
    if (identical(footer, parquet_magic)) return("parquet")
  }
  
  con_text <- file(x, "r")
  lines <- readLines(con_text, n = 5)
  close(con_text)
  if (length(lines) == 0) stop("File empty or cannot be read: ", x)

  tc <- textConnection(lines)
  tab_counts <- count.fields(tc, sep = "\t")
  close(tc)
  if (all(tab_counts >= 2) && length(unique(tab_counts)) == 1) return("tsv")
  
  stop("Not a file in a supported DIA-NN format: ", x)
}

#' @rdname read_diann_proteingroups
#' @export
read_diann <- function(...){
    .Deprecated('read_diann_proteingroups')
    read_diann_proteingroups(...)
}

# file <- download_data('dilution.report.tsv')
# dt <- .read_diann_proteingroups(file)
# dcast_diann(dt, quantity = 'maxlfq',    fill = NA, log2 = TRUE )[1:3, 1:3]
# dcast_diann(dt, quantity = 'precounts', fill = 0               )[1:3, 1:3]
# dcast_diann(dt, quantity = 'pepcounts', fill = 0               )[1:3, 1:3]
# dcast_diann(dt, quantity = 'sequence',  fill = ''              )[1:3, 1:3]
dcast_diann <- function(dt, quantity, fill, log2 = FALSE){
    mat <- data.table::dcast(dt, feature_id ~ run, value.var = quantity, fill = fill)
    mat %<>% dt2mat()
    mat %<>% extract(unique(dt$feature_id), ) # preserve original order
    if (is.na(fill))  mat %<>% zero_to_na() %>% nan_to_na()
    if (log2)         mat %<>% log2()
    mat
}

# No longer publicly accessibly (14.01.2025 )
CONTAMINANTSURL <- paste0(
    'https://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf', 
    '/7994124a4298328fc125748d0048fee2/$FILE/',
    'contaminants.fasta')



#' Rm contaminants
#'
#' Rm contaminants from DIA-NN SumExp
#' @param object         SummarizedExperiment
#' @param verbose        TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('dilution.report.tsv')
#' object <- read_diann_proteingroups(file)
#' object %<>% rm_diann_contaminants()
#' @export
rm_diann_contaminants <- function(object, verbose = TRUE){
# Assert
    assert_is_valid_sumexp(object)
    contaminant <- uniprot <- NULL
# contaminants
    contaminantdt <- read_contaminantdt()
    contaminantdt %<>% extract(uniprot != '')
    contaminants <- contaminantdt$uniprot
# Rm
    fdt0 <- fdt(object)
    fdt0 %<>% separate_rows(uniprot, sep = ';') %>% data.table()
    fdt0[, contaminant := FALSE]
    fdt0[uniprot %in% contaminants, contaminant := TRUE]
    fdt0[, uniprot := NULL]
    fdt0 %<>% extract(, .(contaminant = any(contaminant)), by = 'feature_id')
    object %<>% merge_fdt(fdt0)
    object %<>% filter_features(!contaminant, verbose = verbose)
# Return
    object
}

has_one_level <- function(x) length(unique(x))==1
x <- paste0('pi_exp_', c('wt_r1', 'wt_r2', 'kd_r1', 'kd_r2'))
simplify_snames <- function(x){
    sep <- guess_sep(x)
    if (sep == 'NOSEP')  return(x)
    x <- data.table(x = x)
    x <- x[, tstrsplit(x, split = sep)]
    idx <- !vapply(x, has_one_level, logical(1))
    x <- x[, idx, with = FALSE]
    x <- Reduce(function(a,b) paste(a, b, sep = sep), x)
    x
}

# sampleids <- paste0('pi_exp_', c('wt_r1', 'wt_r2', 'kd_r1', 'kd_r2'))
infer_subgroup <- function(sampleids){
    sep <- guess_sep(sampleids)
    # NoSep: group0
        if (sep == 'NOSEP')  return(rep('group0', length(sampleids)))
    # NoRep: group0
        n <- nfactors(sampleids)
        subgroups <- sampleids %>% split_extract_fixed(sep, 1:(n-1))
        if (!any(duplicated(subgroups)))  return('group0')
    # SepRep: subgroups
        subgroups %<>% factor()
        levels(subgroups) %<>% make.names()  # limma requires!
        return(subgroups)
}

