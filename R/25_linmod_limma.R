
#==============================================================================
#
#                       create_design
#                           singlelevel
#                           multilevel
#                               nlevels
#
#==============================================================================

nlevels <- function(object, svar){
    if (!svar %in% svars(object))  return(0)
    length(unique(object[[svar]]))
}

singlelevel <- function(object, svar)   nlevels(object, svar) ==1
multilevel  <- function(object, svar)   nlevels(object, svar) > 1


#' Does object contain ratio values?
#' @param object SummarizedExperiment
#' @return logical
#' @examples
#' file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' contains_ratios(object)
#'
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' contains_ratios(object)
#' @noRd
contains_ratios <- function(object)  any(grepl('[Rr]atio', assayNames(object)))


#' Create default formula
#' @param object SummarizedExperiment
#' @return formula
#' @examples 
#' # Abundances
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#'     default_formula(object)
# # Ratios
#'     file <- download_data('billing16.proteingroups.txt')
#'     object <- read_maxquant_proteingroups(file)
#'     default_formula(object)
#' @export 
default_formula <- function(object){
    if ('subgroup' %in% svars(object)){
        if (singlelevel(object, 'subgroup')){    return(~1)
        } else if (contains_ratios(object)){      return(~0+subgroup)
        } else {                                  return(~subgroup)
        }
    } else {                                      return(~1)
    }
}

character2factor <- function(x)  if (is.character(x)) factor(x) else x


#' Create design matrix
#'
#' Create design matrix for statistical analysis
#'
#' @param object   SummarizedExperiment or data.frame
#' @param formula  formula with svars
#' @param drop     whether to drop predictor names
#' @param coding   string: codingfunname
#' \itemize{
#'     \item contr.treatment:          intercept = y0,     coefi = yi - y0
#'     \item contr.treatment.explicit: intercept = y0,     coefi = yi - y0
#'     \item code_control:             intercept = ymean,  coefi = yi - y0
#'     \item contr.diff:               intercept = y0,     coefi = yi - y(i-1)
#'     \item code_diff:                intercept = ymean,  coefi = yi - y(i-1)
#'     \item code_diff_forward:        intercept = ymean,  coefi = yi - y(i+)
#'     \item code_deviation:           intercept = ymean,  coefi = yi - ymean (drop last)
#'     \item code_deviation_first:     intercept = ymean,  coefi = yi - ymean (drop first)
#'     \item code_helmert:             intercept = ymean,  coefi = yi - mean(y0:(yi-1))
#'     \item code_helmert_forward:     intercept = ymean,  coefi = yi - mean(y(i+1):yp)
#' }
#' @param verbose      whether to message
#' @param ...          required to s3ify
#' @return design matrix
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' unique(create_design(object))
#' unique(create_design(object, ~ Time))
#' unique(create_design(object, ~ Time, coding = 'code_control'))
#' unique(create_design(object, ~ Time, coding = 'code_diff'))
#' unique(create_design(object, ~ Time + Diabetes))
#' unique(create_design(object, ~ Time / Diabetes))
#' unique(create_design(object, ~ Time * Diabetes))
#' @export
create_design <- function(object, ...) UseMethod('create_design')


#' @rdname create_design
#' @export
create_design.SummarizedExperiment <- function(
    object, 
    formula = default_formula(object),
    drop    = varlevels_dont_clash(object, all.vars(formula)), 
    coding  = 'code_control',
    verbose = TRUE, 
    ...
){
    create_design.data.table(sdt(object), 
                            formula = formula,
                            coding  = coding,
                            drop    = drop,
                            verbose = verbose)
}

#' @rdname create_design
#' @export
create_design.data.table <- function(
    object, 
    formula  = default_formula(object),
    drop     = varlevels_dont_clash(object, all.vars(formula)), 
    coding   = 'code_control',
    verbose  = TRUE, 
    ...
){
# Assert
    assert_is_subset(all.vars(formula), names(object))
    . <- NULL
# Contrast Code Factors
    object %<>% code(coding = coding, vars = all.vars(formula), verbose = verbose)
# Create design matrix
    #if (verbose)   message('\t\tDesign: ', formula2str(formula))
    object %<>% data.frame(row.names = .$sample_id)
    myDesign <- model.matrix(formula, data = object)
    colnames(myDesign) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    is_factor_var <- function(x, object) is.factor(object[[x]])
    if (drop){
        for (predictor in all.vars(formula)){
            if (is.factor(object[[predictor]]))  colnames(myDesign) %<>% 
                        stri_replace_first_fixed(predictor, '') }
            # Fails for e.g. Diabetes = YES/NO: a meaningless column "YES" is created
            # For other cases it works wonderfully, so keep for now.
            # If it gives too many issues, roll back to dropping only for "subgroup" levels:
            # colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
    }
# Return
    return(myDesign)
}



#' Model based prediction
#'
#' @param object   SummarizedExperiment or data.frame
#' @param fit     'limma', 'lm', 'lme', 'wilcoxon'
#' @param formula  formula
#' @param drop     TRUE or FALSE
#' @param coding   string: codingfunname
#' @return beta matrix (nlevel x nfeature)
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% linmod_limma(block = 'Subject', coefs = model_coefs(object)) # intercept required!
#' beta(object)                    #    betas : nlevel x nfeature
#'    X(object)                    #   design : nlevel x nlevel
#'    X(object) %*% beta(object)   # response : nlevel x nfeature
#' @export
X <- function(
    object, 
    formula = default_formula(object),
    drop    = varlevels_dont_clash(object, all.vars(formula)), 
    coding  = 'code_control'
){
    design <- create_design(object, formula = formula, drop = drop, coding = coding, verbose = FALSE)
    X <- unique(design)
    cols <- unique(c('sample_id', all.vars(formula)))
    rownamedt <- data.table(sample_id = rownames(X))
    rownamedt %<>% merge(sdt(object)[ , cols, with = FALSE], by = 'sample_id', sort = FALSE)
    rownames(X) <- rownamedt[ , do.call(function(...) paste(..., sep = '.'), .SD), 
                               .SDcols = names(rownamedt)[-1] ]
    X
}


#' @rdname X
#' @export
beta <- function( object, fit = fits(object)[1] ){
    betas <- effectmat(object, fit = fit, coef = coefs(object, intercept = TRUE))
    colnames(betas) %<>% split_extract_fixed('~', 1)
    if ('Intercept' %in% colnames(betas))  betas[ , 'Intercept' ] <- 0
    betas[ pmat(object, fit = fit) > 0.05 ] <- 0
    betas[ is.na(betas) ] <- 0
    betas %<>% t()
    betas
}



#' Contrast Code Factor
#' 
#' Contrast Code Factor for General Linear Model
#'
#' @param object  factor vector
#' @param vars    svars
#' @param coding  string: codingfunname
#' \itemize{
#'     \item contr.treatment:          intercept = y0,     coefi = yi - y0
#'     \item contr.treatment.explicit: intercept = y0,     coefi = yi - y0
#'     \item code_control:             intercept = ymean,  coefi = yi - y0
#'     \item contr.diff:               intercept = y0,     coefi = yi - y(i-1)
#'     \item code_diff:                intercept = ymean,  coefi = yi - y(i-1)
#'     \item code_diff_forward:        intercept = ymean,  coefi = yi - y(i+)
#'     \item code_deviation:           intercept = ymean,  coefi = yi - ymean (drop last)
#'     \item code_deviation_first:     intercept = ymean,  coefi = yi - ymean (drop first)
#'     \item code_helmert:             intercept = ymean,  coefi = yi - mean(y0:(yi-1))
#'     \item code_helmert_forward:     intercept = ymean,  coefi = yi - mean(y(i+1):yp)
#' }
#' @param verbose TRUE or FALSE
#' @param n character vector
#' @param ... used for s3 dispatch
#' @return (explicitly coded) factor vector
#' @details
#' A General Linear Model contains:                                                                   \cr
#'   * An Intercept Coefficient: expressing some form of sample average                               \cr
#'   * For each numeric variable: a slope coefficient                                                 \cr
#'   * For each k-leveled factor: (k-1) Contrast Coefficients.                                        \cr
#'        The interpretation of (intercept and contrast) coefficients depends on the contrast coding function used.
#'        Several contrast coding functions are available in 'stats' and 'codingMatrices'
#'        But their (function and coefficient) namings are a bit confusing and unsystematic.
#'        Instead, the functions below offer an intuitive interface (to the otherwise powerful stats/codingMatrices packages).
#'        The names of these functions reflect the contrast coding used (treatment, backward, sum, or helmert contrasts).
#'        They also reflect the intercept interpretation (either first factor's first level or grand mean).
#'        They all produce intuitive coefficient names (e.g. 't1-t0' rather than just 't1').
#'        They all have unit scaling (a coefficient of 1 means a backward of 1).
#' @examples
#' # Coding functions
#'     x <- factor(paste0('t', 0:3))
#'     xlevels <- levels(x)
#'     contr.treatment(         xlevels)
#'     contr.treatment.explicit(xlevels)
#'     contr.diff(              xlevels)
#'     code_control(            xlevels)
#'     code_diff(               xlevels)
#'     code_diff_forward(       xlevels)
#'     code_deviation(          xlevels)
#'     code_deviation_first(    xlevels)
#'     code_helmert(            xlevels)
#'     code_helmert_forward(    xlevels)
#' 
#' # Code
#'     x %<>% code('contr.treatment')
#'     x %<>% code('contr.treatment.explicit')
#'     x %<>% code('contr.diff')
#'     x %<>% code('code_control')
#'     x %<>% code('code_diff')
#'     x %<>% code('code_diff_forward')
#'     x %<>% code('code_deviation')
#'     x %<>% code('code_deviation_first')
#'     x %<>% code('code_helmert')
#'     x %<>% code('code_helmert_forward')
#'
#' # Model
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#'     object %<>% linmod_limma(coding = 'contr.treatment') # default
#'     object %<>% linmod_limma(coding = 'contr.treatment.explicit')
#'     object %<>% linmod_limma(coding = 'contr.diff')
#'     object %<>% linmod_limma(coding = 'code_control')
#'     object %<>% linmod_limma(coding = 'code_diff')
#'     object %<>% linmod_limma(coding = 'code_diff_forward')
#'     object %<>% linmod_limma(coding = 'code_deviation')
#'     object %<>% linmod_limma(coding = 'code_deviation_first')
#'     object %<>% linmod_limma(coding = 'code_helmert')
#'     object %<>% linmod_limma(coding = 'code_helmert_forward')
#' @export
code <- function(object, ...)  UseMethod('code')


#' @rdname code
#' @export
code.factor <- function(object, coding, verbose = TRUE, ...){
# Assert
    if (is.null(coding))  return(object)
    assert_is_function(get(coding))
# Code
    k <- length(levels(object))
    contrasts(object) <- get(coding)(levels(object))
    if (verbose){
        contrastmat <- codingMatrices::mean_contrasts(contrasts(object))
        colnames(contrastmat) <- levels(object)
        rownames(contrastmat)[1] <- 'Intercept'
        names(dimnames(contrastmat)) <- c('coefficient', 'level')
        message_df('                    %s', contrastmat)
    }
# Return
    object
}


#' @rdname code
#' @export
code.character <- function(object, coding, verbose = TRUE, ...){
    code.factor(factor(object), coding = coding, verbose = verbose, ...)
}


#' @rdname code
#' @export
code.logical <- function(object, coding, verbose = TRUE, ...){
    code.factor(factor(object, coding = coding, verbose = verbose, ...))
}


#' @rdname code
#' @export
code.numeric <- function(object, coding, verbose = TRUE, ...){
    object
}
    

#' @rdname code
#' @export
code.data.table <- function(object, coding, vars = names(object), verbose = TRUE, ...){
# Assert
    if ( length(vars)==0)   return(object)      # when formula = ~1 
    if (is.null(coding)) return(object)
# Code
    for (var in vars){
        if (verbose)  cmessage('              Code `%s`', var)  # varname only at this level !
        object[[var]] %<>% code(coding, verbose = verbose)
    }
# Return
    object
}


#' @rdname code
#' @export
contr.treatment.explicit <- function(n){
    y <- contr.treatment(n)
    colnames(y) %<>% paste0('-', n[1])
    y
}


#' @rdname code
#' @export
code_control <- function(n){
    if (!installed('codingMatrices'))  return(n) 
    codingMatrices::code_control(n, abbreviate = FALSE)
}


#' @rdname code
#' @export
contr.diff <- function(n){
    if (!installed('codingMatrices'))   return(n) 
    codingMatrices::contr.diff(n, abbreviate = FALSE)
}

#' @rdname code
#' @export
code_diff <- function(n){
    if (!installed('codingMatrices'))  return(n) 
    codingMatrices::code_diff(n, abbreviate = FALSE)
}

#' @rdname code
#' @export
code_diff_forward <- function(n){
    if (!installed('codingMatrices'))  return(n) 
    codingMatrices::code_diff_forward(n, abbreviate = FALSE)
}

#' @rdname code
#' @export
code_deviation <- function(n){
    if (!installed('codingMatrices'))   return(n) 
    k <- length(n)
    contrastnames <- paste0(n, collapse = '+')
    contrastnames <- paste0('(', contrastnames, ')')
    contrastnames <- paste0(contrastnames, '/', length(n))
    contrastnames <- paste0(n[-k], '-', contrastnames) 
    y <- codingMatrices::code_deviation(n)
    colnames(y) <- contrastnames
    y
}

#' @rdname code
#' @export
code_deviation_first <- function(n){
    if (!installed('codingMatrices'))  return(n) 
    k <- length(n)
    contrastnames <- paste0(n, collapse = '+')
    contrastnames <- paste0('(', contrastnames, ')')
    contrastnames <- paste0(contrastnames, '/', length(n))
    contrastnames <- paste0(n[-1], '-', contrastnames) 
    y <- codingMatrices::code_deviation_first(n)
    colnames(y) <- contrastnames
    y
}

#' @rdname code
#' @export
code_helmert <- function(n){
    if (!installed('codingMatrices'))  return(n) 
    y <- codingMatrices::code_helmert(n) # properly scaled version of stats::contr.helmert
    for (i in seq(2, ncol(y)+1)){
        curlevel <- n[i]
        prevlevels <- n[seq(1,i-1)]
        helmertmean <- paste0(prevlevels, collapse = '+')
        if (i>2)  helmertmean <- paste0('(', helmertmean, ')/', i-1)
        colnames(y)[i-1] <- paste0(curlevel, '-', helmertmean)
    }
    y
}

#' @rdname code
#' @export
code_helmert_forward <- function(n){
    if (!installed('codingMatrices'))  return(n) 
    y <- codingMatrices::code_helmert_forward(n) # properly scaled version of stats::contr.helmert
    k <- length(n)
    for (i in seq(1, k-1)){
        curlevel <- n[i]
        nextlevels <- n[seq(i+1,k)]
        fwdmean <- nextlevels
        if (length(nextlevels)>1){
            fwdmean %<>% paste0(collapse = '+')
            fwdmean %<>% paste0('(', ., ')')
            fwdmean %<>% paste0('/', length(nextlevels))
        }
        colnames(y)[i] <- sprintf('%s-%s',  curlevel, fwdmean)
    }
    y
}



#==============================================================================
#
#                            contrast_coefs
#                            model_coefs
#
#==============================================================================


#' @rdname model_coefs
#' @export
contrast_coefs <- function(
       object, 
      formula = default_formula(object), 
         drop = varlevels_dont_clash(object, all.vars(formula)), 
       coding = 'code_control', 
       design = create_design(object, formula = formula, drop = drop, coding = coding, verbose = FALSE)
){
    
    if (ncol(design)==1)  colnames(design) else setdiff(colnames(design), 'Intercept')
}




#' Get model coefs
#' @param object   SummarizedExperiment
#' @param formula  formula
#' @param drop     TRUE or FALSE
#' @param coding   string: codingfunname
#' @param design   design matrix
#' @return SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% linmod_limma()
#'    model_coefs(object)
#' contrast_coefs(object)
#' @export
model_coefs <- function(
    object, 
    formula = default_formula(object), 
       drop = varlevels_dont_clash(object, all.vars(formula)), 
     coding = 'code_control', 
     design = create_design(object, formula = formula, drop = drop, coding = coding, verbose = FALSE)
){
    colnames(design)
}


#==============================================================================
#
#                            limma
#
#==============================================================================


# Old approach - interesting 
#     1. shows the metadata storage approach
#     2. shows how to pull out F, F.p, se
#
#
# .limmacontrast <- function(object, fit, formula){
#     # compute contrasts
#     design <- create_design(object, formula=formula, verbose = FALSE)
#     contrastmat <- makeContrasts(
#         contrasts = vectorize_contrastdefs(contrastdefs(object)),
#         levels    = design)
#     fit %<>% contrasts.fit(contrasts = contrastmat)
#     limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
#     } else { c('effect','rank','t','se','p','fdr','bonf') }
#     limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
#                             dimnames = list(feature  = rownames(fit),
#                                             contrast = colnames(fit),
#                                             quantity = limma_quantities))
#     limma(object)[,,'effect'] <- fit$coefficients
#     limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
#     #names(dimnames(limma(object)))[2] <- formula2str(formula)
#     # perform moderated t test
#     if (!all(fit$df.residual==0)){
#         fit %<>% eBayes()
#         pp <- fit$p.value
#         limma(object)[,,'t' ] <- fit$t
#         limma(object)[,,'se'] <- sqrt(fit$s2.post) * fit$stdev.unscaled
#         limma(object)[,,'p' ] <- pp
#         limma(object)[,,'rank'] <- apply(pp, 2, rank)
#         limma(object)[,,'fdr' ] <- apply(pp, 2, p.adjust, 'fdr')
#         limma(object)[,,'bonf'] <- apply(pp, 2, p.adjust, 'bonferroni')
#         fdata(object)$F.limma   <- fit$F
#         fdata(object)$F.p.limma <- fit$F.p
#     }
# }


#' Reset fit
#' @param object  SummarizedExperiment
#' @param fit     character vector
#' @param verbose TRUE or FALSE
#' @examples 
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %>% fdt()
#' object %>% linmod_limma() %>% fdt()
#' object %>% linmod_limma() %>% reset_fit() %>% fdt()
#' object %>% linmod_limma() %>% linmod_lm() %>% reset_fit('limma') %>% fdt()
#' object %>% linmod_limma() %>% linmod_lm() %>% reset_fit() %>% fdt()
#' @export
reset_fit <- function( object, fit = fits(object), verbose = TRUE ){
# Assert
    . <- NULL
    assert_is_valid_sumexp(object)
    if (is.null(fits(object)))  return(object)
    assert_is_a_bool(verbose)
# Reset fdt
    pattern <- sprintf('~(%s)$', paste0(fit, collapse = '|'))
    cols <- grep(pattern, fvars(object), value = TRUE)
    for (col in cols)  fdt(object)[[col]] <- NULL
    if (length(cols)>0)  if (verbose)  cmessage('%sRm %s', spaces(22), pattern)
# Reset metadata
    metadata(object)$survival <- NULL
# Return
    object
}


# object: SumExp
# fitres: data.table(p.contr1, p.contr2, effect.contr1, effect.contr2)
# stat:  'p', 'effect', 'fdr', 't'
merge_fit <- function(object, fitres, statistic = NULL){
    . <- NULL
    fitresdt <- data.table::copy(fitres)   # dont change in original
    firstcols <- intersect(c('feature_id', 'Intercept'), names(fitresdt))
    fitresdt %<>% extract(,c(firstcols, setdiff(names(.), firstcols)), with = FALSE)
    if (!is.null(statistic)) names(fitresdt)[-1] %<>% paste0(statistic, '~', .)
    object %<>% merge_fdt(fitresdt)
    object
}

mat2fdt <- function(mat)  mat2dt(mat, 'feature_id')

mat2sdt <- function(mat)  mat2dt(mat, 'sample_id')


#' General Linear Model
#'
#' @param object   SummarizedExperiment
#' @param formula  model formula
#' @param engine  'limma', 'lm', 'lme', 'lmer', or 'wilcoxon'
#' @param drop     TRUE or FALSE
#' @param coding   string: codingfunname
#' \itemize{
#'     \item 'contr.treatment':          intercept = y0,     coefi = yi - y0
#'     \item 'contr.treatment.explicit': intercept = y0,     coefi = yi - y0
#'     \item 'code_control':             intercept = ymean,  coefi = yi - y0
#'     \item 'contr.diff':               intercept = y0,     coefi = yi - y(i-1)
#'     \item 'code_diff':                intercept = ymean,  coefi = yi - y(i-1)
#'     \item 'code_diff_forward':        intercept = ymean,  coefi = yi - y(i+)
#'     \item 'code_deviation':           intercept = ymean,  coefi = yi - ymean (drop last)
#'     \item 'code_deviation_first':     intercept = ymean,  coefi = yi - ymean (drop first)
#'     \item 'code_helmert':             intercept = ymean,  coefi = yi - mean(y0:(yi-1))
#'     \item 'code_helmert_forward':     intercept = ymean,  coefi = yi - mean(y(i+1):yp)
#' }
#' @param design    design matrix
#' @param block     block svar (or NULL)
#' @param coefs     NULL or character vector: model coefs to record
#' @param contrasts NULL or character vector: posthoc contrasts to record
#' @param weightvar NULL or name of weight matrix in assays(object)
#' @param reset     TRUE/FALSE whether to wipe earlier modeling results
#' @param suffix    string: pvar suffix ("limma" in "p~t2~limma")
#' @param verbose   whether to msg
#' @param outdir    NULL or dir
#' @param writefun  'write_xl' or 'write_ods'
#' @param plotvolcano  TRUE or FALSE
#' @param plotexprs    TRUE or FALSE
#' @param argsvolcano  list: volcano args
#' @param argsexprs    list:  expr   args
#' @param opt          lme options
#' @return Updated SummarizedExperiment
#' @examples
#' # Standard usage
#'   file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'   object <- read_metabolon(file)
#'   LINMOD(object)                                # Default
#'   LINMOD(object, ~subgroup )                    # Custom formula
#'   LINMOD(object, ~subgroup, block = 'Subject')  # Block effect
#'   
#' # Alternative engines: argument 'engine' or dedicated function
#'   linmod_limma(   object, ~subgroup, block = 'Subject')  # Default engine
#'   linmod_lm(      object, ~subgroup, block = 'Subject')  # Traditional
#'   linmod_lme(     object, ~subgroup, block = 'Subject')  # Powerful random effects
#'   linmod_lmer(    object, ~subgroup, block = 'Subject')  # Yet more powerful random effects
#'   linmod_wilcoxon(object, ~subgroup, block = 'Subject')  # Non-parametric
#'     
#' # Alternative coding: backward diffs instead of baseline
#'   linmod_limma(object, ~ subgroup, block = 'Subject', coding = 'code_diff')
#'   linmod_lme(  object, ~ subgroup, block = 'Subject', coding = 'code_diff')
#'   linmod_lmer( object, ~ subgroup, block = 'Subject', coding = 'code_diff')
#'     
#' # Posthoc contrasts: limma-only, flexible, but sometimes approximate
#'   linmod_limma(object,     ~ subgroup, block = 'Subject', coding = 'code_control')
#'   linmod_limma(object, ~ 0 + subgroup, block = 'Subject', contrasts = 't1-t0')
#'       # flexible, but only approximate
#'       # stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html
#'         
#' # Top-level function also plots and writes
#'   LINMOD(object, block = 'Subject', coefs = 't1-t0')
#'   LINMOD(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE)
#'   LINMOD(object, block = 'Subject', coefs = 't1-t0',   plotexprs = TRUE)
#'   LINMOD(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE, plotexprs = TRUE)
#'   LINMOD(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE, plotexprs = TRUE, outdir = tempdir())
#' @export
LINMOD <- function(
       object, 
      formula = as.formula('~ subgroup'),
       engine = 'limma', 
         drop = varlevels_dont_clash(object, all.vars(formula)),
       coding = 'code_control', # if (engine == 'wilcoxon')  contr.treatment.explicit  else  contr.treatment , 
       design = create_design(object, formula = formula, drop = drop, coding = coding, verbose = FALSE),
        block = NULL,
        coefs = contrast_coefs(object, design = design),
    contrasts = NULL,
    weightvar = if ('weights' %in% assayNames(object)) 'weights'    else NULL,
       suffix = paste0('~', engine),
      verbose = TRUE, 
       outdir = NULL,
     writefun = 'write_xl',
  plotvolcano = FALSE, 
    plotexprs = FALSE, 
  argsvolcano = list(),
    argsexprs = list()
){
# Assert
    assert_scalar_subset(engine, c('limma', 'lme', 'lmer', 'wilcoxon', 'lm'))
    if (!is.null(outdir))  assert_all_are_dirs(outdir)
    assert_scalar_subset(writefun, c('write_xl', 'write_ods'))
    assert_is_a_bool(plotvolcano)
    assert_is_a_bool(plotexprs)
    if (plotvolcano)  assert_is_not_null(coefs)
    if (plotexprs)    assert_is_not_null(coefs)
    assert_is_list(argsvolcano)
    assert_is_list(argsexprs)
# Fit
    if (verbose)  cmessage('%sLinMod', spaces(4)) # unwanted when called during survival
    fitfun <- paste0('linmod_', engine)
    object %<>%  get(fitfun)( formula = formula,
                                 drop = drop,
                               coding = coding, 
                               design = design,
                                block = block, 
                                coefs = coefs,
                            contrasts = contrasts,
                            weightvar = weightvar,
                               suffix = suffix, 
                              verbose = verbose )
# Write tables
    if (!is.null(outdir)){
        outdir <- sprintf('%s/%s', outdir, formula2str(formula))
        dir.create(outdir, showWarnings = FALSE)
    }
    tableext <- switch(writefun, write_xl = 'xlsx', write_ods = 'ods')
    tablefile <- if (is.null(outdir)) NULL else sprintf('%s/tables.%s',    outdir, tableext)
    if (!is.null(outdir))  get(writefun)(object, tablefile) 
# Volcanoes
    if (plotvolcano){
    for (coef in coefs){
        file <- if (is.null(outdir)) NULL else sprintf('%s/%s.volcano.pdf', outdir, coef)
        title <- sprintf('%s', formula2str(formula))
        args <- list( object = object, fit = engine, coefs = coef, title = title, file = file )
        args %<>% c( argsvolcano )
        p <- do.call(plot_volcano, args)
        if (is.null(outdir))  print(p)
    }}
# Exprs
    if (plotexprs){
    for (coef in coefs){
        file <- if (is.null(outdir)) NULL else sprintf('%s/%s.exprs.pdf',   outdir, coef)
        title <- sprintf('%s', formula2str(formula))
        args <- list( object = object,  fit = engine, coefs = coef, title = title,  file = file, block = block )
        args %<>% c( argsexprs )
        p <- do.call(plot_exprs, args)
        if (is.null(outdir))  print(p)
    }}
# Return
    object
}



#' Get all variables from formulas or formula strings
#'
#' An extended version of \code{base::all.vars()} that also accepts character
#' strings (representing formulas), including vectors of such strings.
#' 
#' @param x   formula: scalar/vector with formula/string objects
#' @param ... additional arguments
#' @return    character vector
#' @examples
#' all.vars(   ~Time + Diabetes )
#' all.vars(  '~Time + Diabetes')
#' all.vars(c('~Time + Diabetes', '~Time'))
#' @export
all.vars <- function(x, ...)   UseMethod('all.vars')


# The @method tag is needed ! 
# To ensure dispatch of the function all.vars to class formula
# Rather than dispatching a function all to class vars.formula
#' @rdname all.vars
#' @method all.vars formula
#' @export
all.vars.formula <- function(x, ...)   base::all.vars(x, ...)


# The @method tag is needed ! 
# To ensure dispatch of the function all.vars to class character
# Rather than dispatching a function all to class vars.character
#' @rdname all.vars
#' @method all.vars character
#' @export
all.vars.character <- function(x, ...){
    formulas <- lapply(x, as.formula)
    vars <- lapply(formulas, base::all.vars)
    unique(unlist(vars))
}


#' Are varlevels unique
#' 
#' @param object SummarizedExperiment or data.table
#' @param vars character vector
#' @param ... required for s3 dispatch
#' @return TRUE or FALSE
#' @examples 
#' require(data.table)
#' object1 <- data.table(expand.grid(genome = c('WT', 'MUT'), treat = c('control', 'drug')))
#' object2 <- data.table(expand.grid(mutant = c('YES', 'NO'), treated = c('YES', 'NO')))
#' varlevels_dont_clash(object1)
#' varlevels_dont_clash(object2)
#' @export
varlevels_dont_clash <- function(object, ...)  UseMethod('varlevels_dont_clash')

#' @rdname varlevels_dont_clash
#' @export
varlevels_dont_clash.data.table <- function(
    object, vars = names(object), ...
){
    object                         %>% 
    extract(, vars, with = FALSE)  %>%
    lapply(factor)                 %>% 
    lapply(levels)                 %>% 
    unlist()                       %>% 
    duplicated()                   %>% 
    any()                          %>%
    magrittr::not()
}

#' @rdname varlevels_dont_clash
#' @export
varlevels_dont_clash.SummarizedExperiment <- function(
    object, vars = svars(object), ...
){
    varlevels_dont_clash.data.table(sdt(object), vars)
}


#' @rdname LINMOD
#' @export
linmod_limma <- function(
       object, 
      formula = as.formula('~ subgroup'),
         drop = varlevels_dont_clash(object, all.vars(formula)),
       coding = 'code_control',
       design = create_design(object, formula = formula, drop = drop, coding = coding, verbose = FALSE),
    contrasts = NULL,
        coefs = if (is.null(contrasts))  contrast_coefs(design = design) else NULL,
        block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
        reset = TRUE,
       suffix = '~limma',
      verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_valid_formula(formula, object)
    assert_is_a_bool(drop)
    assert_is_matrix(design)
    if (!is.null(block))      assert_is_subset(block, svars(object))
    if (!is.null(weightvar))  assert_scalar_subset(weightvar, assayNames(object))
    if (length(contrasts)==0 & length(coefs)==0)  return(object)  # awblinmod relies on this
    if (reset)  object %<>% reset_fit(fit = 'limma', verbose = verbose)
# Design/contrasts/block/weights
    . <- NULL
    blockvar <- NULL
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$dupcor)){
            if (verbose)  cmessage('%sDupcor `%s`', spaces(14), blockvar)
            metadata(object)$dupcor <- duplicateCorrelation(values(object), design = design, block = block)$consensus.correlation }
    }
    exprmat <-  values(object)[, rownames(design)]
    weightmat <- if (is.null(weightvar)){ NULL 
            } else {assert_is_a_string(weightvar)
                    assert_is_subset(weightvar, assayNames(object))
                    assays(object)[[weightvar]][, rownames(design)] }
# Fit
    if (verbose)  cmessage("%slinmod_limma( %s%s%s, coding = '%s' )", 
                      spaces(14),
                      formula2str(formula),
                      if(is.null(blockvar))  '' else paste0(' | ',blockvar),
                      if(is.null(weightvar)) '' else paste0(', weights = assays(object)$', weightvar),
                      coding )
    limmafit <- suppressWarnings(lmFit( object = exprmat, design = design, 
                    block = block, correlation = metadata(object)$dupcor, weights = weightmat))
    if (is.null(contrasts)){  limmafit %<>% contrasts.fit(coefficients = model_coefs(design = design))
    } else {                  limmafit %<>% contrasts.fit(contrasts = makeContrasts(contrasts = contrasts, levels = design)) }
    estimable <- !all(limmafit$df.residual==0)
    if (estimable)   limmafit %<>% eBayes()
    
# p/t/fdr
    fitdt <- data.table(feature_id = rownames(limmafit))
    dt0 <- data.table(limmafit$coefficients);         names(dt0) %<>% paste0('effect~', ., suffix); fitdt %<>% cbind(dt0)
    dt0 <- data.table(limmafit$t);                    names(dt0) %<>% paste0(     't~', ., suffix); fitdt %<>% cbind(dt0)
    dt0 <- data.table(limmafit$p.value);              names(dt0) %<>% paste0(     'p~', ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(total = limmafit$df.total);     names(dt0) %<>% paste0(    'df~', ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(prior = limmafit$df.prior);     names(dt0) %<>% paste0(    'df~', ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(resid = limmafit$df.residual);  names(dt0) %<>% paste0(    'df~', ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(sqrt(limmafit$s2.post) * limmafit$stdev.unscaled); names(dt0) %<>% paste0('se~', ., suffix); fitdt %<>% cbind(dt0)
# F statistics                                        # Suprising shorthand for intercept-free fstats !
    cols <- setdiff(colnames(limmafit), 'Intercept')  # https://support.bioconductor.org/p/65253/#65268
    fitdt[, (sprintf('PF~global%s', suffix)) := limmafit[, cols]$F.p.value ]
    fitdt[, (sprintf( 'F~global%s', suffix)) := limmafit[, cols]$F         ]
# Select
    if (is.null(coefs))  coefs <- contrasts
    fitdt %<>% extract(, c(1, which(split_extract_fixed(names(fitdt), '~', 2) %in% coefs)), with = FALSE)
    sumdt <- summarize_fit(fitdt, fit = 'limma', coefs = coefs)
    #if (verbose)  message('')
    if (verbose)  message_df('                          %s', sumdt)
# Return    
  # fdt(object)$F.limma   <- limmafit$F
  # fdt(object)$F.p.limma <- limmafit$F.p
    object %<>% merge_fdt(fitdt)
    object
}


#' @rdname LINMOD
#' @export
fit_limma <- function(...){ .Deprecated('linmod_limma'); linmod_limma(...)}


pull_level <- function(x, lev){
    assert_is_factor(x)
    if (lev %in% levels(x))  x %<>% 
        factor(levels = c(lev, setdiff(levels(x), lev)))
    x
}


#' Summarize fit
#' @param object  SummarizedExperiment or data.table
#' @param fit  'limma', 'lme', 'lm', 'lme', 'wilcoxon' or NULL
#' @param coefs string vector
#' @param ... S3 dispatch
#' @return data.table(contrast, nup, ndown)
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% linmod_limma()
#' object %<>% linmod_lm()
#' summarize_fit(object, coefs = c('t1-t0', 't2-t0', 't3-t0'))
#' @export
summarize_fit <- function(object, ...)  UseMethod('summarize_fit')


#' @rdname summarize_fit
#' @export
summarize_fit.data.table <- function(
    object, 
       fit = fits(object),
     coefs = autonomics::coefs(object, fit = fit), 
       ...
){
# Assert
    object %<>% copy()
    if (is.null(coefs))  return(NULL)
    statistic <- coefficient <- variable <- NULL
    effect <- p <- fdr <- NULL
# Summarize
    cols <- names(object) %>% extract(stri_detect_fixed(., '~'))
    object %<>% extract(, c('feature_id', cols), with = FALSE)
    object %<>% add_adjusted_pvalues(method = 'fdr', fit = fit, coefs = coefs)
    assert_has_no_duplicates(names(object))
        # Good to make sure!
        # Because if there are duplicate cols then the dcasting further down is no longer unique
        # And dcasting then resorts to meaningless length aggregation
    longdt <- object %>% melt.data.table(id.vars = 'feature_id')
    longdt[, statistic    := split_extract_fixed(variable, '~', 1) %>% factor(unique(.))]
    longdt[,  coefficient := split_extract_fixed(variable, '~', 2) %>% factor(unique(.))]
    longdt[,       fit    := split_extract_fixed(variable, '~', 3) %>% factor(unique(.))]
    longdt[, variable := NULL]
    
    sumdt <- dcast.data.table(longdt, feature_id + coefficient + fit ~ statistic, value.var = 'value')
    sumdt <- sumdt[, .(
        downfdr = sum(t < 0  & fdr < 0.05, na.rm = TRUE), 
        upfdr   = sum(t > 0  & fdr < 0.05, na.rm = TRUE),
        downp   = sum(t < 0  &   p < 0.05, na.rm = TRUE), 
        upp     = sum(t > 0  &   p < 0.05, na.rm = TRUE)), by = c('coefficient', 'fit') ]
    if (!is.null(fit)){
        idx <- sumdt$fit %in% fit
        sumdt %<>% extract(idx)
    }
    if (!is.null(coefs)){
        sumdt <- sumdt[coefficient %in% coefs]
    }
    sumdt
}


#' @rdname summarize_fit
#' @export
summarize_fit.SummarizedExperiment <- function(
    object, fit = fits(object), coefs = autonomics::coefs(object, fit = fit), ...
){
    summarize_fit.data.table(fdt(object), fit = fit, coefs = coefs)
}


#' Plot fit summary
#' @param sumdt data.table
#' @param nrow number
#' @param ncol number
#' @param order TRUE or FALSE
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% linmod_lm()
#' object %<>% linmod_limma(block = 'Subject')
#' sumdt <- summarize_fit(object, coefs = c('t1-t0', 't2-t0', 't3-t0'))
#' plot_fit_summary(sumdt)
#' @export
plot_fit_summary <- function(sumdt, nrow = NULL, ncol = NULL, order = FALSE){
    coefficient <- downfdr <- downp <- fit <- upfdr <- upp <- NULL
    if (order){
        sumdt <- sumdt[order(downfdr+upfdr, downp+upp)]
        sumdt[, coefficient := factor(coefficient, unique(coefficient))]
    }
    ggplot(sumdt) + facet_wrap(vars(fit), nrow = nrow, ncol = ncol) + 
    geom_col(aes(y = coefficient, x = -downp),   fill = 'firebrick',   alpha = 0.3) +
    geom_col(aes(y = coefficient, x =    upp),   fill = 'forestgreen', alpha = 0.3) + 
    geom_col(aes(y = coefficient, x = -downfdr), fill = 'firebrick',   alpha = 1) +
    geom_col(aes(y = coefficient, x =    upfdr), fill = 'forestgreen', alpha = 1) + 
    geom_text(data = sumdt[  downp>0], aes(y = coefficient, x = -max(downp), label = paste0(downp, ' | ', downfdr) ), hjust = +1) + 
    geom_text(data = sumdt[    upp>0], aes(y = coefficient, x =    max(upp), label = paste0(upfdr, ' | ', upp) ), hjust = 0) + 
    xlab('count') + 
    ylab(NULL) + 
    xlim(c(-max(sumdt$downp)-100, max(sumdt$upp)+100)) + 
    #scale_x_continuous(n.breaks = 20) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


setna <- function(dt, value){
   for (j in seq_len(ncol(dt)))  set(dt, which(is.na(dt[[j]])), j, value)
    dt
}


#' formula to string
#' @param formula formula
#' @return string
#' @examples 
#' formula2str(~0+subgroup)
#' @export
formula2str <- function(formula)  Reduce(paste, deparse(formula))


