#=============================================================================
#
#                   guess_sep
#                   guess_sep.SummarizedExperiment
#                   guess_sep.factor
#                   guess_sep.character
#                       has_identical_values
#                       is_max
#                           cequals
#
#=============================================================================

#' Convenient equals operator
#'
#' Performs x == y, but returns FALSE rather than NA for NA elements of x.
#' @param x numeric vector or scalar
#' @param y numeric scalar
#' @return logical vector
#' @examples
#' x <- c(A = 1, B = 3, C = 2, D = 3, E = NA)
#' y <- 3
#' cequals(x, y)
#' @noRd
cequals <- function(x,y){
    result <- rep(FALSE, length(x)) %>% set_names(names(x))
    if (is.na(y)){
        result[ is.na(x)] <- TRUE
        result[!is.na(x)] <- FALSE
    } else {
        result[ is.na(x)] <- FALSE
        result[!is.na(x)] <- x[!is.na(x)] == y
    }
    result
}


#' Is maximal
#' @param x numeric vector
#' @return logical vector
#' @examples
#' x <- c(A = 1,B = 3,C = 2,D = 3, E = NA)
#' is_max(x)
#' @noRd
is_max <- function(x) cequals(x, max(x, na.rm = TRUE))


#' All elements of vector are identical
#' @param x vector
#' @return TRUE or FALSE
#' @examples
#' x <- c(2,2,1,2)
#' has_identical_values(x)
#' @noRd
has_identical_values <- function(x) length(unique(x))==1

#' Guess separator
#' @param x          character vector or SummarizedExperiment
#' @param var        svar or fvar
#' @param separators character vector: possible separators to look for
#' @param verbose    TRUE or FALSE
#' @param ...        used for proper S3 method dispatch
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' # charactervector
#'    guess_sep(c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]'))
#'    guess_sep(c('WT_untreated_1', 'WT_untreated_2'))
#'    guess_sep(c('group1', 'group2.R1'))
#' # SummarizedExperiment
#'    file <- download_data('atkin.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    guess_sep(object)
#' @export
guess_sep <- function (x, ...)  UseMethod("guess_sep", x)


#' @rdname guess_sep
#' @export
guess_sep.numeric <- function(x, ...) NULL


#' @rdname guess_sep
#' @export
guess_sep.character <- function(
    x, separators = c('.', '_'), verbose = FALSE, ...
){
# Initialize
    . <- NULL
    sep_freqs <-Map(function(y) stri_split_fixed(x, y), separators) %>%
                lapply(function(y) vapply(y, length, integer(1)))            %>%
                extract( vapply(., has_identical_values, logical(1)))        %>%
                vapply(unique, integer(1))
# No separator detected - return NULL
    if (all(sep_freqs==1)){
        if (verbose) message(x[1],': no (consistent) separator. Returning NULL')
        return('NOSEP')   # no separator detected
    }
# Find best separator
    best_sep <- sep_freqs %>%
                extract(.!=1)  %>%
                extract(is_max(vapply(., extract, integer(1), 1)))   %>%
                names()
# Ambiguous separator - take first from tail
    if (length(best_sep)>1){
        pattern <- best_sep %>% paste0(collapse='') %>% paste0('[', ., ']')
        best_sep <- x[1] %>% stri_extract_last_regex(pattern)
    }
# Separator identified - return
    if (verbose) message("\t\tGuess sep: '", best_sep, "'")
    return(best_sep)
}


#' @rdname guess_sep
#' @export
guess_sep.factor <- function(x, ...)  guess_sep.character(levels(x))


#' @rdname guess_sep
#' @export
guess_sep.SummarizedExperiment <- function(
    x, var = 'sample_id', separators =  c('.', '_'), verbose = FALSE, ...
){
    assert_is_subset(var, c(svars(x), fvars(x)))
    (if (var %in% svars(x)) slevels(x, var) else flevels(x, var)) %>%
    guess_sep(separators = separators, verbose = verbose)
}



#=============================================================================
#
#                       nfactors
#                       split_extract
#
#=============================================================================

#' @export
#' @rdname split_extract_fixed
nfactors <- function(x, sep = guess_sep(x)){
    length(unlist(stri_split_fixed(x[1], sep)))
}

#' stri_split and extract
#' @param x    character vector
#' @param sep  string
#' @param i    integer
#' @return character vector
#' @examples
#' # Read
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     x <- object$sample_id[1:5]
#'     nfactors(x)
#' # Split
#'     split_extract_fixed(x, '_', 1:2)
#'     split_extract_fixed(x, '_', seq_len(nfactors(x)-1))
#'     split_extract_fixed(x, '_', nfactors(x))
#'     split_extract_fixed(fdt(object)$PUBCHEM, ';', 1) # with NA values
#' @export
split_extract_fixed <- function(x, sep, i){
    if (is.factor(x))  x %<>% as.character()
    assert_is_a_string(sep)
    assert_is_numeric(i)
    factors <- stri_split_fixed(x, sep)
    extracted <- rep(NA_character_, length(factors))
    idx <- !is.na(factors)
    extracted[idx] <- vapply(factors[idx], 
                     function(y) paste0(y[i], collapse=sep), character(1))
    extracted
}

#' @export
#' @rdname split_extract_fixed
split_extract_regex <- function(x, sep, i){
    if (is.factor(x))  x %<>% as.character()
    assert_is_a_string(sep)
    assert_is_numeric(i)
    factors <- stri_split_regex(x, sep)
    extracted <- rep(NA_character_, length(factors))
    idx <- !is.na(factors)
    extracted[idx] <- vapply(factors[idx], 
                             function(y) paste0(y[i], collapse=sep), character(1))
    extracted
}

#' @export
#' @rdname split_extract_fixed
split_extract <- function(x, i, sep=guess_sep(x)){
    .Deprecated('split_extract_fixed')
    factors <- stri_split_fixed(x, sep)
    extracted <- rep(NA_character_, length(factors))
    idx <- !is.na(factors)
    extracted[idx] <- vapply(
                        factors[idx], 
                        function(y) paste0(y[i], collapse=sep), character(1))
    extracted
}


#=============================================================================
#
#                   default_sfile
#
#=============================================================================

#' Default sfile
#' @param file data file
#' @return sample file
#' @examples
#' file <- download_data('billing19.proteingroups.txt')
#' default_sfile(file)
#' @export
default_sfile <- function(file){
    sfile <- tools::file_path_sans_ext(file)
    sfile %<>% paste0('.samples.txt')
    sfile
}


#=============================================================================
#
#               subgroup_matrix
#                   split_subgroup_levels
#                       split_subgroup_values
#                           split_values
#
#=============================================================================

split_values <- function(x){
    sep <- guess_sep(x)
    dt <- data.table::data.table(x = x)
    dt[, data.table::tstrsplit(x, sep) ]
}

split_subgroup_values <- function(object, subgroupvar){
    subgroupvalues <- svalues(object, subgroupvar)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}

split_subgroup_levels <- function(object, subgroupvar){
    subgrouplevels <- slevels(object, subgroupvar)
    cbind(subgroup = subgrouplevels, split_values(subgrouplevels))
}


#' @rdname subgroup_matrix
#' @export
subgroup_array <- function(object, subgroupvar){
    . <- NULL
    x <- slevels(object, subgroupvar)
    sep <- guess_sep(object)
    #x %<>% sort()
    dt <- data.table(subgroup = factor(x, x))
    components <- dt[, tstrsplit(subgroup, sep, fixed=TRUE)]
    for (i in seq_len(ncol(components)))   components[[i]] %<>%
                                    factor(., levels=unique(.))
    dt %<>% cbind(components)
    data.table::setorderv(dt, rev(names(components)))
    levels  <- dt[, -1] %>% lapply(unique)
    #levels[1:2] %<>% rev()
    nlevels <- levels %>% vapply(length, integer(1))
    array(dt$subgroup, dim = nlevels, dimnames = levels)
}


#' Get subgroup matrix
#'
#' Arrange (subgroup)levels in matrix
#'
#' @param object SummarizedExperiment
#' @param subgroupvar subgroup svar
#' @return matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' subgroup_matrix(object, 'subgroup')
#' @export
subgroup_matrix <- function(object, subgroupvar){
    . <- NULL
    subgroup_array <- subgroup_array(object, subgroupvar)
    if (length(dim(subgroup_array)) == 1){
        return(matrix(subgroup_array, 
                      byrow = TRUE, 
                      nrow = 1, 
                      dimnames = list(NULL, subgroup_array)))  
    }
    otherdims <- names(dim(subgroup_array)) %>% setdiff('V1')
    ncol1   <- Reduce('*', dim(subgroup_array)[otherdims])
    colnames1 <- dimnames(subgroup_array)[otherdims] %>%
                 expand.grid()                       %>%
                 apply(1, paste0, collapse = '.')
    subgroupmat <- matrix(subgroup_array,
                        nrow = nrow(subgroup_array), ncol = ncol1,
                        dimnames=list(rownames(subgroup_array), colnames1))
    subgroupmat %>% extract(nrow(.):1, )
    #dt <- split_subgroup_levels(object)
    #subgroupmat <- as.matrix(data.table::dcast(
    #    dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    #subgroupmat %>% extract(rev(order(rownames(.))), order(colnames(.)))
}

#------------------------------------------------------------------------------
#
#                   modelvar
#
#------------------------------------------------------------------------------


#' Get model variable
#' @param object          SummarizedExperiment
#' @param quantity        'p', 'effect', 'fdr', 't', or 'se'
#' @param fit             string (vector)
#' @param coef            string (vector)
#' @param fvar            'feature_id' or other fvar for values (pvec) or names (upfeatures)
#' @param significancevar 'p' or 'fdr'
#' @param significance     p or fdr cutoff (fractional number)
#' @param effectdirection  '<>', '<' or '>'
#' @param effectsize      effectsize cutoff (positive number)
#' @return string (tvar), matrix (tmat), numeric vector (tvec), character vector (tfeatures)
#' @examples 
#' # Read
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     object %<>% fit_limma()
#'     object %<>% fit_lm()
#'     featuredt <- fdt(object)
#' # modelvar
#'     modelvar(featuredt, 'p');                                          pvar(featuredt)
#'     modelvar(featuredt, 'effect');                                effectvar(featuredt)
#'     modelvar(featuredt, 'fdr');                                      fdrvar(featuredt)
#' # modelvec
#'     modelvec(featuredt, 'p'     )[1:3];                                 pvec(featuredt)[1:3]
#'     modelvec(featuredt, 'effect')[1:3];                            effectvec(featuredt)[1:3]
#'     modelvec(featuredt, 'fdr'   )[1:3];                               fdrvec(featuredt)[1:3]
#' # modelmatrix
#'     modelmat(featuredt, 'p'     )[1:3, 1:3];                            pmat(featuredt)[1:3, 1:3]
#'     modelmat(featuredt, 'effect')[1:3, 1:3];                       effectmat(featuredt)[1:3, 1:3]
#'     modelmat(featuredt, 'fdr'   )[1:3, 1:3];                          fdrmat(featuredt)[1:3, 1:3]
#' # modelfeatures
#'     modelfeatures(featuredt      )[1:3]
#'     modelfeatures(featuredt, effectdirection = '<' )[1:3];      downfeatures(featuredt)[1:3]
#'     modelfeatures(featuredt, effectdirection = '>' )[1:3];        upfeatures(featuredt)[1:3]
#' @export
modelvar <- function(
    featuredt, 
     quantity, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit)
){
# Assert
    assert_is_data.table(featuredt)
    assert_is_subset(quantity, c('fdr', 'p', 't', 'effect', 'se', 'abstract'))
    assert_is_subset(fit,   fits(featuredt))
    assert_is_subset(coef, coefs(featuredt, fit = fit))
# Return
    sep <- guess_fitsep(featuredt)
    x <- expand.grid(quantity = quantity, fit = fit, coef = coef)
    x <-  paste(x$quantity, x$coef, x$fit, sep = sep)
    x %<>% intersect(names(featuredt))        # fits dont always contain same coefs: 
    if (length(x)==0)  x <- NULL           # `limma(contrasts)` mostly without intercept
    x   # NULL[1] and c('a', NULL) work!   # `lm(coefs)` mostly with intercept               
}
    
#' @rdname modelvar
#' @export
effectvar <- function( 
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit) 
){
    modelvar(featuredt, quantity = 'effect', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
tvar <- function(
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit) 
){
    modelvar(featuredt, quantity = 't', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
pvar <- function( 
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit) 
){
    modelvar(featuredt, quantity = 'p', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
fdrvar <- function( 
    featuredt, 
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    modelvar(featuredt, quantity = 'fdr', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
abstractvar <- function(
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit) 
){
    # cant use modelvar because its
    sep <- guess_fitsep(featuredt)
    y <- paste(coef, fit, sep = sep)        #           t1~limma
    if (y %in% names(featuredt)) y else NULL   #     not p~t1~limma
}

#------------------------------------------------------------------------------
#
#                   modelvec
#
#------------------------------------------------------------------------------

#' @rdname modelvar
#' @export
modelvec <- function(
    featuredt, 
     quantity, 
          fit = fits(featuredt)[1], 
         coef = default_coefs(featuredt, fit = fit)[1], 
         fvar = 'feature_id'
){
    valuevar <- modelvar(featuredt, quantity = quantity, fit = fit, coef = coef)
    if (is.null(valuevar))  return(NULL)
    y        <- featuredt  %>%  extract2(valuevar)
    names(y) <- featuredt  %>%  extract2(fvar)
    y
}

#' @rdname modelvar
#' @export
effectvec <- function(
    featuredt, 
          fit  = fits(featuredt)[1], 
          coef = default_coefs(featuredt, fit = fit)[1], 
          fvar = 'feature_id'
){
    modelvec(featuredt, quantity = 'effect', fit = fit, coef = coef, fvar = fvar)
}

#' @rdname modelvar
#' @export
tvec <- function(
    featuredt, 
          fit = fits(featuredt)[1],
         coef = default_coefs(featuredt, fit = fit)[1], 
         fvar = 'feature_id'
){
    modelvec(featuredt, quantity = 't', fit = fit, coef = coef, fvar = fvar)
}
    
#' @rdname modelvar
#' @export
pvec <- function(
    featuredt, 
          fit = fits(featuredt)[1], 
         coef = default_coefs(featuredt, fit = fit)[1], 
         fvar = 'feature_id'
){
    modelvec(featuredt, quantity = 'p', fit = fit, coef = coef, fvar = fvar)
}

#' @rdname modelvar
#' @export
fdrvec <- function(
     featuredt, 
           fit = fits(featuredt)[1],
          coef = default_coefs(featuredt, fit = fit)[1], 
          fvar = 'feature_id'
){
    modelvec(featuredt, quantity = 'fdr', fit = fit, coef = coef, fvar = fvar)
}

#------------------------------------------------------------------------------
#
#                   modeldt
#
#------------------------------------------------------------------------------


#' @rdname modelvar
#' @export
modeldt <- function(
    featuredt, 
     quantity, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit)
){
    var <- modelvar(featuredt, quantity, coef = coef, fit = fit)
    if (is.null(var))  return(NULL)
    dt <- featuredt[, c('feature_id', var), with = FALSE]
    #names(dt) %<>% stri_replace_first_fixed(paste0(var, sep), '')
    dt
}


#' @rdname modelvar
#' @export
effectdt <- function(
    featuredt, 
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    modeldt(featuredt, quantity = 'effect', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
tdt <- function(
    featuredt, 
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    modeldt(featuredt, quantity = 't', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
pdt <- function(
    featuredt, 
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    modeldt(featuredt, quantity = 'p', fit = fit, coef = coef)
}



#------------------------------------------------------------------------------
#
#                   modelmat
#
#------------------------------------------------------------------------------


#' @rdname modelvar
#' @export
modelmat <- function(
    featuredt, 
     quantity,
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    dt <- modeldt(featuredt, quantity, coef = coef, fit = fit)
    mat <- dt2mat(dt)
    mat
}

#' @rdname modelvar
#' @export
effectmat <- function(
    featuredt, 
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    modelmat(featuredt, quantity = 'effect', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
effectsizemat <- function(
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit)
){
    abs(modelmat(featuredt, quantity = 'effect', fit = fit, coef = coef))
} # dont rm: used in ..extract_statistic_features : 
  # getFromNamespace(sprintf('%smat', statistic), 'autonomics')

#' @rdname modelvar
#' @export
tmat <- function(
    featuredt, 
          fit = fits(featuredt), 
         coef = default_coefs(featuredt, fit = fit)
){
    modelmat(featuredt, quantity = 't', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
pmat <- function(
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit)
){
    modelmat(featuredt, quantity = 'p', fit = fit, coef = coef)
}

#' @rdname modelvar
#' @export
fdrmat <- function(
    featuredt, 
          fit = fits(featuredt),
         coef = default_coefs(featuredt, fit = fit)
){
    modelmat(featuredt, quantity = 'fdr', fit = fit, coef = coef)
}


#==============================================================================
#
#                   modelfeatures
#
#==============================================================================


#' @rdname modelvar
#' @export
modelfeatures <- function(
          featuredt,
                fit = fits(featuredt)[1],
               coef = default_coefs(featuredt, fit = fit)[1], 
               fvar = 'feature_id', 
    significancevar = 'p',
       significance = 0.05,
    effectdirection = '<>',
         effectsize = 0
){
    significancevalues <- modelvec(featuredt, quantity =  significancevar, fit = fit, coef = coef, fvar = fvar)
          effectvalues <- modelvec(featuredt, quantity = 'effect',         fit = fit, coef = coef, fvar = fvar)
    idx <- switch(effectdirection, 
                  `<>` = (effectvalues < -effectsize) | (effectvalues > +effectsize), 
                  `<`  = (effectvalues < -effectsize), 
                  `>`  = (effectvalues > +effectsize) )
    idx <- idx & (significancevalues < significance)
    y <- featuredt[[fvar]][idx]
    unique(y)
}

#' @rdname modelvar
#' @export
upfeatures <- function(
          featuredt, 
                fit = fits(featuredt)[1], 
               coef = default_coefs(featuredt, fit = fit)[1], 
               fvar = 'feature_id',
    significancevar = 'p',
       significance = 0.05,
         effectsize = 0
){
    modelfeatures(    featuredt = featuredt, 
                            fit = fit,
                           coef = coef,
                           fvar = fvar, 
                significancevar = significancevar,
                   significance = significance, 
                effectdirection = '>',
                     effectsize = effectsize )
}

#' @rdname modelvar
#' @export
downfeatures <- function(
          featuredt,
                fit = fits(featuredt)[1], 
               coef = default_coefs(featuredt, fit = fit)[1], 
               fvar = 'feature_id',
    significancevar = 'p',
       significance = 0.05,
         effectsize = 0
){
    modelfeatures(  featuredt = featuredt,
                          fit = fit,
                         coef = coef,
                         fvar = fvar, 
              significancevar = significancevar,
                 significance = significance, 
              effectdirection = '>', 
                   effectsize = effectsize )
}


#----------------------------------------------------------------------------------------


# dont rm - its the lower-level function used by fits() and coefs() !
.effectvars <- function(featuredt){
    . <- NULL
    sep <- guess_fitsep(featuredt)
    names(featuredt) %>% extract(stri_startswith_fixed(., paste0('effect', sep)))
}

#' Get fit models
#' 
#' @param featuredt feature data.table
#' @return  character vector
#' @examples 
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma')
#' fits(fdt(object))
#' @export
fits <- function(featuredt){
    sep <- guess_fitsep(featuredt)
    if (is.null(sep))  return(NULL)
    x <- .effectvars(featuredt)
    x %<>% split_extract_fixed(sep, 3)
    x %<>% unique()
    x
}

#' Get coefs
#' 
#' @param featuredt  data.table or factor
#' @param fit     string: 'limma', 'lm', 'lme', 'lmer'
#' @param svars   NULL or charactervector (svar for which to return coefs)
#' @param ...     required for s3 dispatch
#' @return  character vector
#' @examples
#' # Factor
#'     x <- factor(c('A', 'B', 'C'))
#'     coefs(x)
#'     coefs(code(x, contr.treatment.explicit))
#'     coefs(code(x, code_control))
#'     
#' # SummarizedExperiment
#'     file <- download_data('atkin.metabolon.xlsx')
#'     object <- read_metabolon(file, fit = 'limma')
#'     coefs(fdt(object))
#' @export
coefs <- function(featuredt, ...)  UseMethod('coefs')

#' @rdname coefs
#' @export
coefs.factor <- function(featuredt, ...)   colnames(contrasts(featuredt))

#' @rdname coefs
#' @export
coefs.data.table <- function(
    featuredt, fit = fits(featuredt), svars = NULL, ...
){
    sep <- guess_fitsep(featuredt)
    if (is.null(sep))  return(NULL)
    if (is.null(fit))  return(NULL)
    . <- NULL
    coefs0 <- split_extract_fixed(.effectvars(featuredt), sep, 2)
    fits0  <- split_extract_fixed(.effectvars(featuredt), sep, 3)
    coefs0 %<>% extract(fits0 %in% fit)
    coefs0 %<>% unique()
    #if (!is.null(svars))  coefs0 %<>% extract(Reduce('|', lapply(svars, grepl, .)))
    coefs0 
}

#============================================================================
#
#                   is_sig  .is_sig
#                   plot_contrast_venn
#
#============================================================================

.is_sig <- function(object, fit, contrast, quantity = 'fdr'){
    sigfun <- get(paste0(quantity, 'mat'))
    issig  <- sigfun(object) < 0.05
    isdown <- issig & effectmat(fdt(object)) < 0
    isup   <- issig & effectmat(fdt(object)) > 0
    isdown[is.na(isdown)] <- FALSE
    isup[  is.na(isup)  ] <- FALSE
    testmat <- matrix(0, nrow(issig), ncol(issig), dimnames=dimnames(issig))
    testmat[isdown] <- -1
    testmat[isup]   <-  1
    col <- modelvar(object, fit = fit, coef = contrast, quantity = quantity)
    testmat[, col, drop = FALSE]
}


#' Is significant?
#' @param object    SummarizedExperiment
#' @param fit       subset of autonomics::TESTS
#' @param contrast  subset of colnames(metadata(object)[[fit]])
#' @param quantity  value in dimnames(metadata(object)[[fit]])[3]
#' @return matrix: -1 (downregulated), +1 (upregulatd), 0 (not fdr significant)
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' object %<>% fit_lm()
#' object %<>% fit_limma()
#' issig <- is_sig(object, fit = c('lm','limma'), contrast = 'Adult')
#' plot_contrast_venn(issig)
#' @export
is_sig <- function(
      object,
         fit = fits( fdt(object))[1],
    contrast = coefs(fdt(object)),
    quantity = 'fdr'
){
# Assert
    . <- NULL
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_character(fit)
    assert_is_subset(fit, fits(fdt(object)))
    if (is.character(contrast))  for (fi in fit)  assert_is_subset(contrast, coefs(fdt(object), fit = fi))
# Run across models
    res <-  mapply(.is_sig, fit, 
                    MoreArgs = list(object = object, contrast = contrast, quantity = quantity),
                    SIMPLIFY = FALSE)
    #add_model_names <- function(isfdrmat, model){
    #                    colnames(isfdrmat) %<>% paste(model, sep='.')
    #                    isfdrmat }
    #if (length(fit)>1){
    #    res %<>% mapply(add_model_names , ., names(.), SIMPLIFY=FALSE)
    #}
    res %<>% do.call(cbind, .)
    res
}


