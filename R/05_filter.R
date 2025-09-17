
#=================
# FILTER FEATURES
#=================


#' Filter features on condition
#' @param object SummarizedExperiment
#' @param condition filter condition
#' @param verbose logical
#' @return filtered eSet
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' filter_features(object, SUPER_PATHWAY == 'Lipid')
#' @export
filter_features <- function(object, condition, verbose = TRUE){
    . <- NULL
    condition <- enquo(condition)
    idx <- eval_tidy(condition, fdata(object))
    idx <- idx & !is.na(idx)
    if (verbose & sum(idx) < length(idx)){
        cmessage('%sRetain %d/%d features: %s', spaces(14), sum(idx), length(idx),
                expr_text(condition) %>% substr(1, min(120, nchar(.))))}
    object %<>% extract(idx,)
    fdata(object) %<>% droplevels()
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>%
            c(structure(sum(idx), names = quo_name(condition)))
    }
    object
}


#' @rdname rm_missing_in_some_samples
#' @export
rm_missing_in_all_samples <- function(object, verbose = TRUE){
    # . != 0 needed due to stupid behaviour of rowAnys
    # https://github.com/HenrikBengtsson/matrixStats/issues/89
    selector <- rowAnys(values(object) != 0, na.rm = TRUE)
    if (verbose && sum(selector)<length(selector)){
        cmessage('%sRetain %d/%d features: non-zero, non-NA, and non-NaN for some sample', 
                 spaces(14), sum(selector), length(selector))
        object %<>% extract(selector, )
        if (!is.null(analysis(object))) {
            analysis(object)$nfeatures %<>% c(structure(sum(selector),
                names = "non-zero, non-NA, and non-NaN for some sample"))
        }
    }
    object
}

is_available_in_all_samples <- function(object)  rowAlls(!is.na(values(object)))


#' Rm features missing in some samples
#' @param object SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return updated object
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' rm_missing_in_all_samples( object)
#' rm_missing_in_some_samples(object)
#' @export
rm_missing_in_some_samples <- function(object, verbose = TRUE){

    # Restrict to available values
    selector <- is_available_in_all_samples(object)
    if (verbose)  message('\t\t\tUse ', sum(selector), '/', length(selector),
                            ' features with available value for each sample')
    object %<>% extract(selector, )
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(sum(selector),
            names = "available value for each sample"))
    }
    object
}


#==================
# FILTER EXPRS
#==================

#' Filter features with replicated expression in some subgroup
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @param assay        string
#' @param comparator   '>' or '!='
#' @param lod          number: limit of detection
#' @param nsample      number
#' @param nsubgroup    number
#' @param verbose      TRUE or FALSE
#' @return Filtered SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% filter_exprs_replicated_in_some_subgroup()
#' filter_exprs_replicated_in_some_subgroup(object, character(0))
#' filter_exprs_replicated_in_some_subgroup(object, NULL)
#' @export
filter_exprs_replicated_in_some_subgroup <- function(
    object, subgroupvar = 'subgroup', assay = assayNames(object)[1],
    comparator = if (contains_ratios(object)) '!=' else '>',
    lod = 0, nsample = 2, nsubgroup = 1, verbose = TRUE
){
# Assert
    assert_is_subset(subgroupvar, svars(object))
    replicated <- NULL
# Datatablify
    replicated_in_its_subgroup <- replicated_in_any_subgroup <- value <- NULL
    dt <- sumexp_to_longdt(object, svars = subgroupvar, assay = assay)
# Find replicated features
    exceeds_lod <- if (comparator == '>'){ function(value, lod) value >  lod
            } else if (comparator == '!=') function(value, lod) value != lod
    
    condition <- sprintf('value %s %s', comparator, lod)
    repfeatures <- dt[,  .(replicated = .(sum(eval(parse(text = condition)), na.rm = TRUE) >= nsample)) , by = c('feature_id', subgroupvar)]
    repfeatures %<>% extract( , .(replicated = sum(as.numeric(replicated)) >= 1), by = 'feature_id')
    repfeatures %<>% extract(replicated == TRUE)
    repfeatures %<>% extract2('feature_id')
    repfeatures %<>% as.character()
# Keep only replicated features
    idx <- fid_values(object) %in% repfeatures
    if (verbose)  if (any(!idx))  cmessage('\t\tFilter %d/%d features: %s %s %s for at least %d samples in %d %s', sum(idx), 
            length(idx), assay, comparator, as.character(lod), nsample, nsubgroup, subgroupvar)
    object %<>% extract(idx, )
# Update analysis log
    if (!is.null(analysis(object))) {
        analysis(object)$nfeatures %<>% c(structure(
                sum(idx),
                names = sprintf(
                    "expr %s %s, for at least two samples in some %s",
                    comparator, as.character(lod), subgroupvar)))
    }
    object
}


#' Are coefs/pvalues estimable
#' @param formula formula
#' @param data    data.table
#' @examples
#' # Onevar design
#' # -------------
#'     # Design not full rank, coefficients/pvalues not estimable
#'          (dt <- data.table( time = factor(c('t0', 't1', 't2', 't3'  ) ), 
#'                            value =        c(  0,    1,    2,   NA   ) ))
#'             coefs_estimable(~time, data = dt)
#'           pvalues_estimable(~time, data = dt)
#'             summary(lm(value~time, data = dt))
#' 
#'     # Design full rank, coefficients estimable.
#'     # No residual dof, pvalues not estimable.
#'          (dt <- data.table( time = factor(c('t0', 't1', 't2', 't3'  ) ), 
#'                            value =        c(  0,    1,    2,    3   ) ))
#'             coefs_estimable(~time, data = dt)
#'           pvalues_estimable(~time, data = dt)
#'             summary(lm(value~time, data = dt))
#' 
#'     # Design full rank, coefficients estimable
#'     # Residual dof, pvalues estimable
#'          (dt <- data.table( time = factor(c('t0', 't1', 't2', 't3', 't3' ) ), 
#'                            value =        c(  0,    1,    2,    3,    3.1) ))
#'             coefs_estimable(~time, data = dt)
#'           pvalues_estimable(~time, data = dt)
#'           summary(lm(value~time, data = dt))
#'
#' # Twovar design
#' # -------------
#'     # Design not full rank, coefficients/pvalues not estimable.
#'          (dt <- data.table( time = factor(c(  't0',  't1', 't2', 't2','t3',  't3',  't0', 't1', 't2', 't3' ) ), 
#'                         diabetes = factor(c(   'C',  'C',  'C',  'C',  'C',   'C',   'D',  'D',  'D',  'D'  ) ),
#'                            value =        c(    0,    1,    2,    2.1,  3,    3.1,   NA,   NA,   NA,   NA  ) ))
#'             coefs_estimable(~time+diabetes, data = dt)
#'           pvalues_estimable(~time+diabetes, data = dt)
#'           # summary(lm(value~time+diabetes, data = dt))
#' 
#'     # Design full rank, coefficients estimable
#'     # No residual dof, pvalues not estimable
#'          (dt <- data.table( time = factor(c( 't0', 't1', 't2', 't3',  't0', 't1', 't2', 't3' ) ), 
#'                         diabetes = factor(c(  'C', 'C',  'C',  'C',   'D',  'D',  'D',  'D'  ) ),
#'                            value =        c(   0,   1,    2,   3,      0.5,  NA,   NA,   NA  ) ))
#'             coefs_estimable(~time+diabetes, data = dt)
#'           pvalues_estimable(~time+diabetes, data = dt)
#'             summary(lm(value~time+diabetes, data = dt))
#' 
#'     # Design full rank, coefficients estimable
#'     # Residual dof, pvalues estimable
#'          (dt <- data.table( time = factor(c( 't0', 't1', 't2', 't3',  't0', 't1', 't2', 't3' ) ), 
#'                         diabetes = factor(c(  'C', 'C',  'C',  'C',   'D',  'D',  'D',  'D'  ) ),
#'                            value =        c(   0,    1,    2,   3,     0.5,  1.6,  NA,   NA  ) ))
#'             coefs_estimable(~time+diabetes, data = dt)
#'           pvalues_estimable(~time+diabetes, data = dt)
#'             summary(lm(value~time+diabetes, data = dt))
#' @export
pvalues_estimable <- function(formula, data){
    data %<>% extract(!is.na(value))
    design <- model.matrix(formula, data = data)
    qr(design)$rank == ncol(design) &   # Full rank design (coefficient estimation)
    qr(design)$rank <  nrow(design)     #        Residual dof (variance estimation)
}


#' @rdname pvalues_estimable
#' @export
coefs_estimable <- function(formula, data){
    data %<>% extract(!is.na(value))
    design <- model.matrix(formula, data = data)
    qr(design)$rank == ncol(design)   # Full rank desing (coefficient estimation)
}

#' Keep estimable features
#' @param object  SummarizedExperiment
#' @param formula model formula
#' @param coding  coding function name 
#' @param verbose TRUE or FALSE
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' keep_estimable_features(object, ~ subgroup)
#' @export
keep_estimable_features <- function(object, formula = ~1, coding = 'code_control', verbose = TRUE){
    sdt(object) %<>% code(coding = coding, vars = all.vars(formula))
    estimdt <- sumexp_to_longdt(object, svars = all.vars(formula))
    estimdt <- estimdt[, .(estimable = pvalues_estimable(formula, .SD)), by = 'feature_id']
    object %<>% merge_fdt(estimdt)
    object %<>% filter_features(estimable == TRUE, verbose = verbose)
    object
}#


#' Keep fully connected blocks
#' @param object  SummarizedExperiment
#' @param block   svar
#' @param verbose TRUE or FALSE
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' keep_connected_blocks(object, block = 'Subject')          # autonomics format
#' keep_connected_blocks(object, block = list(Subject = ~1)) # lme format
#' @export
keep_connected_blocks <- function(object, block, verbose = TRUE){
    if (is.null(block))  return(object)
    if (is.list(block)) block <- names(block)  # linmod_lme(object, ~ Diabetes + Time, block = list(Subject = ~1))
    all_blocks <- unique(object[[block]])
    full_blocks <- sdt(object)[, .N, by = block][N==max(N)][[block]]
    idx <- object[[block]] %in% full_blocks
    if (sum(idx) < length(idx)){
        if (verbose)  cmessage('%sKeep %d/%d fully connected blocks with %d/%d samples',
                         spaces(14), length(full_blocks), length(all_blocks), sum(idx), length(idx))
        object %<>% extract(, idx)
    }
    object
}


#' Keep features with n+ connected blocks
#' @param object   SummarizedExperiment
#' @param block    svar
#' @param n        number
#' @param verbose  TRUE or FALSE
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' keep_connected_features(object, block = 'Subject')
#' keep_connected_features(object, block = list(Subject = ~1))
#' @export
keep_connected_features <- function(object, block, n = 2, verbose = TRUE){
    if (is.null(block))  return(object)
    if (is.list(block)) block <- names(block)  # linmod_lme(object, ~ Diabetes + Time, block = list(Subject = ~1))
    dt <- sumexp_to_longdt(object, svars = block)
    nperblock <- dt[, .N, by = c('feature_id', block)][, max(N)]
    n0 <- length(unique(dt$feature_id))

    idx  <- dt[, .N, by = c('feature_id', block)]                    #   nobs per feature per block
    idx %<>% extract(, .(N = sum(N==nperblock)), by = 'feature_id')  #   n completeblocks per feature 
    idx %<>% extract(N>=n)                                           #   2 completeblocks per feature
    idx %<>% extract(, feature_id)

    if (length(idx) < length(unique(dt$feature_id))){
        if (verbose)  cmessage('\t\t\tRetain %d/%d features: 2+ fully connected blocks',
                                length(idx), length(unique(dt$feature_id)))
        object %<>% extract(feature_id %in% idx, )
    }
    object
    # This earlier approach fails in ~ Time / Diabetes
    # Because a subject cannot be both diabetic AND control
    # obj %<>% extract_connected_features(formula = formula, blockvars = block, verbose = verbose) # doesnt work for complex models
}


#' Tag features
#' 
#' @param object    SummarizedExperiment
#' @param keyvar    string : intersection fvar
#' @param sep       string : keyvar collapse separator
#' @param features  character vector : intersection set
#' @param tagvar    string : 
#' @param verbose   TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
#' object <- read_somascan(file)
#' features <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = 'SYMBOL')
#' object %<>% tag_features(keyvar = 'EntrezGeneSymbol', sep = ' ', features)
#' table(fdt(object)$features)
#' @export
tag_features <- function(
    object, keyvar, sep, features, tagvar = get_name_in_parent(features), verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(keyvar, fvars(object))
    assert_is_a_string(sep)
    assert_is_character(features)
    assert_is_a_string(tagvar)
    assert_is_a_bool(verbose)
# Intersect
    cols <- unique(c('feature_id', keyvar))
    fdt0 <- fdt(object)[, cols, with = FALSE ]
    fdt0 %<>% uncollapse(tidyselect::all_of(keyvar), sep = sep)
    fdt0 %<>% extract(get(keyvar) %in% features)
# Filter
    idx <- fdt(object)$feature_id %in% unique(fdt0$feature_id)
    fdt(object)[[tagvar]]      <- FALSE
    fdt(object)[[tagvar]][idx] <- TRUE
    object
}


#=======================
# FILTER SAMPLES
#=======================


#' Filter samples on condition
#' @param object    SummarizedExperiment
#' @param condition filter condition
#' @param verbose   TRUE/FALSE 
#' @param record    TRUE/FALSE 
#' @param drop      TRUE/FALSE : whether to drop levels
#' @return filtered SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' filter_samples(object, subgroup != 't0', verbose = TRUE)
#' @export
filter_samples <- function(object, condition, verbose = TRUE, record = TRUE, drop = TRUE){
    . <- NULL
    condition <- enquo(condition)
    idx <- eval_tidy(condition, sdt(object))
    idx <- idx & !is.na(idx)
    if (verbose & sum(idx)<length(idx)){
        cmessage('%sRetain %d/%d samples: %s', spaces(14), sum(idx), length(idx), 
                expr_text(condition) %>% substr(1, min(120, nchar(.))))}
    object %<>% extract(, idx)
    if (drop)  sdt(object) %<>% droplevels()
    if (record && !is.null(analysis(object))) {
        analysis(object)$nsamples %<>%  
            c(structure(sum(idx), names = expr_text(condition)))
    }
    object
}


#' Filter samples available for some feature
#' @param object SummarizedExperiment
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @noRd
filter_samples_available_for_some_feature <- function(object, verbose = FALSE){
    subsetter <- is_available_for_some_feature(object)
    if (any(!subsetter)){
        if (verbose)  message('\t\t\tRetain ', sum(subsetter), '/', 
                            length(subsetter),
                            ' samples with a value available for some feature')
    object %<>% extract(, subsetter)
    }
    object
}

is_available_for_some_feature <- function(object){
    subsetter <- (!is.na(values(object))) & (values(object) != 0)
    set_names(colAnys(subsetter), snames(object))
}



