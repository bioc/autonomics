
#==============================================================================
#
#                       create_design
#                           single_subgroup
#                           are_factor
#                           singlelevel
#                           multilevel
#                               nlevels
#
#==============================================================================

single_subgroup <- function(object){
    assert_is_subset('subgroup', svars(object))
    length(unique(object$subgroup))==1
}


are_factor <- function(df) vapply(df, is.factor, logical(1))

nlevels <- function(object, svar){
    if (!svar %in% svars(object))  return(0)
    length(unique(object[[svar]]))
}

singlelevel <- function(object, svar)   nlevels(object, svar) ==1
multilevel  <- function(object, svar)   nlevels(object, svar) > 1

#' @rdname default_formula
#' @export
default_subgroupvar <- function(object){
    if ('subgroup' %in% svars(object))  'subgroup' else NULL
}

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
#' @param object       SummarizedExperiment or data.frame
#' @param formula      formula with svars
#' @param drop         whether to drop predictor names
#' @param codingfun  factor coding function
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
#' unique(create_design(object, ~ Time, codingfun = code_control))
#' unique(create_design(object, ~ Time, codingfun = code_diff))
#' unique(create_design(object, ~ Time + Diabetes))
#' unique(create_design(object, ~ Time / Diabetes))
#' unique(create_design(object, ~ Time * Diabetes))
#' @export
create_design <- function(object, ...) UseMethod('create_design')


#' @rdname create_design
#' @export
create_design.SummarizedExperiment <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)), 
    codingfun = code_control,
    verbose   = TRUE, 
    ...
){
    create_design.data.table(sdt(object), 
                            formula   = formula,
                            codingfun = codingfun,
                            drop      = drop,
                            verbose   = verbose)
}

#' @rdname create_design
#' @export
create_design.data.table <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)), 
    codingfun = code_control,
    verbose   = TRUE, 
    ...
){
# Assert
    assert_is_subset(all.vars(formula), names(object))
    . <- NULL
# Contrast Code Factors
    object %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
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
            # For other cases it works wonderfully, so I keep it for now.
            # If it gives too many issues, roll back to doing the dropping only
            # for "subgroup" levels:
            #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
    }
# Return
    return(myDesign)
}



#' Model based prediction
#'
#' @param object     SummarizedExperiment or data.frame
#' @param fit        'limma', 'lm', 'lme', 'wilcoxon'
#' @param formula    formula
#' @param drop       TRUE or FALSE
#' @param codingfun  function
#' @return beta matrix (nlevel x nfeature)
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% limma(block = 'Subject') # intercept required!
#' beta(object)                    #    betas : nlevel x nfeature
#'    X(object)                    #   design : nlevel x nlevel
#'    X(object) %*% beta(object)   # response : nlevel x nfeature
#' @export
X <- function(
    object, 
    formula   = default_formula(object),
    drop      = varlevels_dont_clash(object, all.vars(formula)), 
    codingfun = code_control
){
    design <- create_design(object, formula = formula, drop = drop, codingfun = codingfun, verbose = FALSE)
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
    sep <- guess_fitsep(object)
    colnames(betas) %<>% split_extract_fixed(sep, 1)
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
#' @param vars svars
#' @param codingfun  factor coding function
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
#'     x %<>% code(contr.treatment)
#'     x %<>% code(contr.treatment.explicit)
#'     x %<>% code(contr.diff)
#'     x %<>% code(code_control)
#'     x %<>% code(code_diff)
#'     x %<>% code(code_diff_forward)
#'     x %<>% code(code_deviation)
#'     x %<>% code(code_deviation_first)
#'     x %<>% code(code_helmert)
#'     x %<>% code(code_helmert_forward)
#'
#' # Model
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#'     object %<>% limma(codingfun = contr.treatment) # default
#'     object %<>% limma(codingfun = contr.treatment.explicit)
#'     object %<>% limma(codingfun = contr.diff)
#'     object %<>% limma(codingfun = code_control)
#'     object %<>% limma(codingfun = code_diff)
#'     object %<>% limma(codingfun = code_diff_forward)
#'     object %<>% limma(codingfun = code_deviation)
#'     object %<>% limma(codingfun = code_deviation_first)
#'     object %<>% limma(codingfun = code_helmert)
#'     object %<>% limma(codingfun = code_helmert_forward)
#' @export
code <- function(object, ...)  UseMethod('code')


#' @rdname code
#' @export
code.factor <- function(object, codingfun, verbose = TRUE, ...){
# Assert
    if (is.null(codingfun))  return(object)
    assert_is_function(codingfun)
# Code
    k <- length(levels(object))
    contrasts(object) <- codingfun(levels(object))
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
code.character <- function(object, codingfun, verbose = TRUE, ...){
    code.factor(factor(object), codingfun = codingfun, verbose = verbose, ...)
}


#' @rdname code
#' @export
code.logical <- function(object, codingfun, verbose = TRUE, ...){
    code.factor(factor(object, codingfun = codingfun, verbose = verbose, ...))
}


#' @rdname code
#' @export
code.numeric <- function(object, codingfun, verbose = TRUE, ...){
    object
}
    

#' @rdname code
#' @export
code.data.table <- function(object, codingfun, vars = names(object), verbose = TRUE, ...){
# Assert
    if ( length(vars)==0)   return(object)      # when formula = ~1 
    if (is.null(codingfun)) return(object)
# Code
    for (var in vars){
        if (verbose)  cmessage('              Code `%s`', var)  # varname only at this level !
        object[[var]] %<>% code(codingfun, verbose = verbose)
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

#=============================================================================
#
#               contrast_coefs
#                   contrast_subgroup_cols
#                   contrast_subgroup_rows
#
#==============================================================================


#' Row/Col contrasts
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @return  matrix
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object$subgroup <- paste0(object$Diabetes, '.', object$Time)
#' subgroup_matrix(object, subgroupvar = 'subgroup')
#' contrast_subgroup_cols(object, subgroupvar = 'subgroup')
#' contrast_subgroup_rows(object, subgroupvar = 'subgroup')
#' @export
contrast_subgroup_cols <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s-%s',    # no space: as lm(contr.sdiff)
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s-%s',   # no space: as lm(contr.sdiff)
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s-%s',  # no space: as lm(contr.sdiff)
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s-%s', # no space: as lm(contr.sdiff)
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}


# contrast_coefs <- function(object, formula){
#     subgroupvar <- all.vars(formula)[1]
#     design <- create_design(object, formula = formula)
#     if (ncol(design)==1){
#         list(matrix(colnames(design), nrow=1, ncol=1), 
#             matrix(nrow=0, ncol=0))
#     } else if (all(design[, 1]==1)){
#         list(colnames(design)[-1][seq_len(nlevels(object, subgroupvar)-1)], 
#             matrix(nrow=0, ncol=0))
#     } else {
#         list(contrast_subgroup_cols(object, subgroupvar),
#             contrast_subgroup_rows( object, subgroupvar)) }
# }


#' @rdname model_coefs
#' @export
contrast_coefs <- function(
       object, 
      formula = default_formula(object), 
         drop = varlevels_dont_clash(object, all.vars(formula)), 
    codingfun = code_control, 
       design = create_design(object, formula = formula, drop = drop, codingfun = codingfun, verbose = FALSE)
){
    
    if (ncol(design)==1)  colnames(design) else setdiff(colnames(design), 'Intercept')
}


#' Get model coefs
#' @param object     SummarizedExperiment
#' @param formula    formula
#' @param drop       TRUE or FALSE
#' @param codingfun  coding function (e.g. contr.treatment)
#' @param design     design matrix
#' @return SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %<>% limma()
#' model_coefs(object)
#' contrast_coefs(object)
#' @export
model_coefs <- function(
    object, 
    formula = default_formula(object), 
       drop = varlevels_dont_clash(object, all.vars(formula)), 
  codingfun = code_control, 
     design = create_design(object, formula = formula, drop = drop, codingfun = codingfun, verbose = FALSE)
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

contrvec2mat  <- function(contrasts)  matrix(
                    contrasts, nrow=1, dimnames=list("", contrasts))

contrmat2list <- function(contrasts)  list(colcontrasts = contrasts)

vectorize_contrasts <- function(contrasts){
    unname(unlist(lapply(contrasts, function(x) na.exclude(c(t(x))))))
}


#' Reset fit
#' @param object  SummarizedExperiment
#' @param fit     character vector
#' @param verbose TRUE or FALSE
#' @examples 
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object %>% fdt()
#' object %>% limma() %>% fdt()
#' object %>% limma() %>% reset_fit() %>% fdt()
#' object %>% limma() %>% lm.() %>% reset_fit('limma') %>% fdt()
#' object %>% limma() %>% lm.() %>% reset_fit() %>% fdt()
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


#' Fit results separator
#' @examples
#' FITSEP
#' @export
FITSEP <- '~'


#' @rdname FITSEP
#' @export
PPATTERN <- paste0('p', FITSEP)


# object: SumExp
# fitres: data.table(p.contr1, p.contr2, effect.contr1, effect.contr2)
# stat:  'p', 'effect', 'fdr', 't'
merge_fit <- function(object, fitres, statistic = NULL){
    . <- NULL
    fitresdt <- data.table::copy(fitres)   # dont change in original
    firstcols <- intersect(c('feature_id', 'Intercept'), names(fitresdt))
    fitresdt %<>% extract(,c(firstcols, setdiff(names(.), firstcols)), with = FALSE)
    if (!is.null(statistic)) names(fitresdt)[-1] %<>% paste0(statistic,FITSEP,.)
    object %<>% merge_fdt(fitresdt)
    object
}

mat2fdt <- function(mat)  mat2dt(mat, 'feature_id')

mat2sdt <- function(mat)  mat2dt(mat, 'sample_id')


#' Formulate 
#'
#' Formulate model
#' @param modelvars vector
#' @param across    TRUE or FALSE
#' @param within    TRUE or FALSE
#' @param between   TRUE or FALSE
#' @return string
#' @examples
#' formulate(   'subgroup' )
#' formulate( c('Time', 'Diabetes'), across  = TRUE)
#' formulate( c('Time', 'Diabetes'), within  = TRUE)
#' formulate( c('Time', 'Diabetes'), between = TRUE)
#' formulate( c('Time', 'Diabetes'), across = TRUE, within  = TRUE, between = TRUE)
#' @export
formulate <- function(modelvars, across = FALSE, within = FALSE, between = FALSE){
    
    assert_is_character(modelvars)
    assert_is_a_bool(across)
    assert_is_a_bool(within)
    assert_is_a_bool(between)
    if (length(modelvars)>1)  assert_any_are_true(c(across, within, between))
    
    formula <- character(0)
    if (length(modelvars) == 1          ){  formula %<>% c(sprintf('~ %s', modelvars)                             ); names(formula)[length(formula)] <- 'default'                                     }
    if (length(modelvars) == 2 & across ){  formula %<>% c(sprintf('~ %s', paste0(    modelvars,  collapse = '+'))); names(formula)[length(formula)] <- 'across'                                      }
    if (length(modelvars) == 2 & within ){  formula %<>% c(sprintf('~ %s', paste0(    modelvars,  collapse = '/'))); names(formula)[length(formula)] <- paste0(rev(modelvars), collapse = '.within.') }
    if (length(modelvars) == 2 & within ){  formula %<>% c(sprintf('~ %s', paste0(rev(modelvars), collapse = '/'))); names(formula)[length(formula)] <- paste0(    modelvars,  collapse = '.within.') }
    if (length(modelvars) == 2 & between){  formula %<>% c(sprintf('~ %s', paste0(    modelvars,  collapse = '*'))); names(formula)[length(formula)] <- 'between'                                     }
    if (length(modelvars) >  2          ){  message('`modelvars` limited to two variables - use `formula` instead')}
    return(formula)
}


#' General Linear Model
#'
#' @param object    SummarizedExperiment
#' @param formula   model formula
#' @param engine    'limma', 'lm', 'lme', 'lmer', or 'wilcoxon'
#' @param drop      TRUE or FALSE
#' @param codingfun  factor coding function
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
#' @param design    design matrix
#' @param block     block svar (or NULL)
#' @param coefs     NULL or character vector: model coefs to record
#' @param contrasts NULL or character vector: posthoc contrasts to record
#' @param weightvar NULL or name of weight matrix in assays(object)
#' @param sep       string: pvar separator  ("~" in "p~t2~limma")
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
#'   object %<>% linmod()                                # Default
#'   object %<>% linmod(  ~subgroup )                    # Custom formula
#'   object %<>% linmod(  ~subgroup, block = 'Subject')  # Block effect
#'   summarize_fit(object)
#'   
#' # Alternative engines: argument 'engine' or dedicated function
#'   fdt(object) %<>% extract(, 'feature_id')
#'   object %<>% limma(   ~subgroup, block = 'Subject')  # Default engine
#'   object %<>% lm.(     ~subgroup, block = 'Subject')  # Traditional
#'   object %<>% lme(     ~subgroup, block = 'Subject')  # Powerful random effects
#'   object %<>% lmer(    ~subgroup, block = 'Subject')  # Yet more powerful random effects
#'   object %<>% wilcoxon(~subgroup, block = 'Subject')  # Non-parametric
#'   summarize_fit(object)
#'     
#' # Alternative coding: backward diffs instead of baseline
#'   fdt(object) %<>% extract(, 'feature_id')
#'   object %<>% limma(     ~ subgroup, block = 'Subject', codingfun = code_diff)
#'   object %<>% lme(       ~ subgroup, block = 'Subject', codingfun = code_diff)
#'   object %<>% lmer(      ~ subgroup, block = 'Subject', codingfun = code_diff)
#'   summarize_fit(object)
#'     
#' # Posthoc contrasts: limma-only, flexible, but sometimes approximate
#'   fdt(object) %<>% extract(, 'feature_id')
#'   object %<>% limma( ~ subgroup, block = 'Subject', codingfun = code_control)
#'   object %<>% limma( ~ 0 + subgroup, block = 'Subject', contrasts = 't1-t0')
#'       # flexible, but only approximate
#'       # stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html
#'         
#' # Custom separator
#'   fdt(object) %<>% extract(, 'feature_id')
#'   fdt( limma(object, sep = '.'))
#'   fdt( limma(object, block = 'Subject', sep = '.') )
#'
#' # Top-level function also plots and writes
#'   linmod(object, block = 'Subject', coefs = 't1-t0')
#'   linmod(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE)
#'   linmod(object, block = 'Subject', coefs = 't1-t0',   plotexprs = TRUE)
#'   linmod(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE, plotexprs = TRUE)
#'   linmod(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE, plotexprs = TRUE, outdir = tempdir())
#'   linmod(object, block = 'Subject', coefs = 't1-t0', plotvolcano = TRUE, plotexprs = TRUE, outdir = tempdir())
#' @export
linmod <- function(
       object, 
      formula = as.formula('~ subgroup'),
       engine = 'limma', 
         drop = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = code_control, # if (engine == 'wilcoxon')  contr.treatment.explicit  else  contr.treatment , 
       design = create_design(object, formula = formula, drop = drop, codingfun = codingfun, verbose = FALSE),
        block = NULL,
        coefs = contrast_coefs(object, design = design),
    contrasts = NULL,
    weightvar = if ('weights' %in% assayNames(object)) 'weights'    else NULL,
          sep = FITSEP,
       suffix = paste0(sep, engine),
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
    fitfun <- paste0('fit_', engine)
    object %<>%  get(fitfun)( formula = formula,
                                 drop = drop,
                            codingfun = codingfun, 
                               design = design,
                                block = block, 
                                coefs = coefs,
                            contrasts = contrasts,
                            weightvar = weightvar,
                                  sep = sep,
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


#' @rdname linmod
#' @export
fit_limma <- function(...){ .Deprecated('limma'); limma(...)}


#' @rdname linmod
#' @export
fit_lm <- function(...){ .Deprecated('lm.'); lm.(...)}


#' @rdname linmod
#' @export
fit_lme <- function(...){ .Deprecated('lme'); lme(...)}


#' @rdname linmod
#' @export
fit_lmer <- function(...){ .Deprecated('lmer'); lmer(...)}


#' @rdname linmod
#' @export
fit_wilcoxon <- function(...){ .Deprecated('wilcoxon'); wilcoxon(...)}



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


#' General Linear Model (awb interface)
#' @param object    SummarizedExperiment
#' @param engine 'limma', 'lm', 'lme', or 'lmer'
#' @param modelvars svars
#' @param across    TRUE/FALSE: fit across  model (additive)    ?
#' @param within    TRUE/FALSE: fit within  model (nested)      ?
#' @param between   TRUE/FALSE: fit between model (interaction) ?
#' @examples
#' object <- survobj()
#' svars(object)
#' awb(object, engine = 'limma', modelvars = c('age', 'sex'))
#' awb_limma(object, modelvars = c('age', 'sex'), block = 'replicate')
#' awb_lm(   object, modelvars = c('age', 'sex'), block = 'replicate')
#' awb_lme(  object, modelvars = c('age', 'sex'), block = 'replicate')
#' awb_lmer( object, modelvars = c('age', 'sex'), block = 'replicate')
awb <- function(
    object, 
    engine,
 modelvars,
    across = TRUE,
    within = TRUE,
   between = TRUE,
 codingfun = code_control,
      drop = TRUE, 
          ...
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(modelvars, svars(object))

# Model
    modelfun <- switch(engine, limma = gen_limma, lm = gen_lm, lme = gen_lme, lmer = gen_lmer)
    if (across){
        formula  <- paste0(modelvars, collapse = '+')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs  <- colnames(create_design(object,  formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% setdiff('Intercept')
        object %<>% modelfun(formula,  codingfun = codingfun, drop = drop, coefs =  coefs, ...)
    }
    if (within){
        formula <- paste0(modelvars, collapse = '/')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs <- colnames(create_design(object, formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% extract(stri_detect_fixed(., ':'))
        object %<>% modelfun(formula, codingfun = codingfun, drop = drop, coefs = coefs, ...)
        fvars(object) %<>% stri_replace_first_fixed(':', '/')
    }
    if (within){
        formula <- paste0(rev(modelvars), collapse = '/')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs <- colnames(create_design(object, formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% extract(stri_detect_fixed(., ':'))
        object %<>% modelfun(formula, codingfun = codingfun, drop = drop, coefs = coefs, ...)
        fvars(object) %<>% stri_replace_first_fixed(':', '/')
    }
    if (between){
        formula <- paste0(modelvars, collapse = '*')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs <- colnames(create_design(object, formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% extract(stri_detect_fixed(., ':'))
        object %<>% modelfun(formula, codingfun = codingfun, drop = drop, coefs = coefs, ...)
        fvars(object) %<>% stri_replace_first_fixed(':', '*')
    }
# Return
    object
}


#' @rdname awb
#' @export
awb_limma <- function(object, ...)  limoawb(object, engine = 'limma', ...)


#' @rdname awb
#' @export
awb_lm    <- function(object, ...)  limoawb(object, engine = 'lm', ...)


#' @rdname awb
#' @export
awb_lme   <- function(object, ...)  limoawb(object, engine = 'lme', ...)


#' @rdname awb
#' @export
awb_lmer  <- function(object, ...)  limoawb(object, engine = 'lmer', ...)


  
#' @rdname linmod
#' @export
gen_limma <- function(
       object, 
      formula = as.formula('~ subgroup'),
         drop = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = code_control,
       design = create_design(object, formula = formula, drop = drop, codingfun = codingfun),
    contrasts = NULL,
        coefs = if (is.null(contrasts))  contrast_coefs(design = design) else NULL,
        block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
          sep = FITSEP,
       suffix = paste0(sep, 'limma'),
      verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_valid_formula(formula, object)
    assert_is_a_bool(drop)
    assert_is_matrix(design)
    if (!is.null(block))      assert_is_subset(block, svars(object))
    if (!is.null(weightvar))  assert_scalar_subset(weightvar, assayNames(object))
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
    if (verbose)  cmessage('%slmFit(%s%s%s)', 
                    spaces(14),
                    formula2str(formula),
                    if(is.null(blockvar))  '' else paste0(' | ',blockvar),
                    if(is.null(weightvar)) '' else paste0(', weights = assays(object)$', weightvar))
    limmafit <- suppressWarnings(lmFit( object = exprmat, design = design, 
                    block = block, correlation = metadata(object)$dupcor, weights = weightmat))
    if (is.null(contrasts)){  limmafit %<>% contrasts.fit(coefficients = model_coefs(design = design))
    } else {                  limmafit %<>% contrasts.fit(contrasts = makeContrasts(contrasts = contrasts, levels = design)) }
    estimable <- !all(limmafit$df.residual==0)
    if (estimable)   limmafit %<>% eBayes()
    
# p/t/fdr
    fitdt <- data.table(feature_id = rownames(limmafit))
    dt0 <- data.table(limmafit$coefficients);                            names(dt0) %<>% paste0('effect',  sep, ., suffix); fitdt %<>% cbind(dt0)
    dt0 <- data.table(limmafit$t);                                       names(dt0) %<>% paste0('t',       sep, ., suffix); fitdt %<>% cbind(dt0)
    dt0 <- data.table(limmafit$p.value);                                 names(dt0) %<>% paste0('p',       sep, ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(total = limmafit$df.total);                        names(dt0) %<>% paste0('df',      sep, ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(prior = limmafit$df.prior);                        names(dt0) %<>% paste0('df',      sep, ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(resid = limmafit$df.residual);                     names(dt0) %<>% paste0('df',      sep, ., suffix); fitdt %<>% cbind(dt0)
   #dt0 <- data.table(sqrt(limmafit$s2.post) * limmafit$stdev.unscaled); names(dt0) %<>% paste0('se',      sep, ., suffix); fitdt %<>% cbind(dt0)
# F statistics                                        # Suprising shorthand for intercept-free fstats !
    cols <- setdiff(colnames(limmafit), 'Intercept')  # https://support.bioconductor.org/p/65253/#65268
    fitdt[, (sprintf('PF%sglobal%s', sep, suffix)) := limmafit[, cols]$F.p.value ]
    fitdt[, (sprintf( 'F%sglobal%s', sep, suffix)) := limmafit[, cols]$F         ]
# Select
    if (is.null(coefs))  coefs <- contrasts
    fitdt %<>% extract(, c(1, which(split_extract_fixed(names(fitdt), '~', 2) %in% coefs)), with = FALSE)
    sumdt <- summarize_fit(fitdt, fit = 'limma')
    if (verbose)  message_df('                  %s', sumdt)
# Return    
  # fdt(object)$F.limma   <- fitdt$F
  # fdt(object)$F.p.limma <- fitdt$F.p
    object %<>% merge_fdt(fitdt)
    object
}


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
#' object %<>% limma()
#' object %<>% lm.()
#' summarize_fit(object, coefs = c('t1-t0', 't2-t0', 't3-t0'))
#' @export
summarize_fit <- function(object, ...)  UseMethod('summarize_fit')


#' @rdname summarize_fit
#' @export
summarize_fit.data.table <- function(
    object, 
          fit = fits(object),
        coefs = autonomics::coefs(object, fit = fit), ...
){
# Assert
    object %<>% copy()
    if (is.null(coefs))  return(NULL)
    statistic <- coefficient <- variable <- NULL
    effect <- p <- fdr <- NULL
# Summarize
     sep <- guess_fitsep(object)
    cols <- names(object) %>% extract(stri_detect_fixed(., sep))
    object %<>% extract(, c('feature_id', cols), with = FALSE)
    object %<>% add_adjusted_pvalues(method = 'fdr', fit = fit, coefs = coefs)
    assert_has_no_duplicates(names(object))
        # Good to make sure!
        # Because if there are duplicate cols then the dcasting further down is no longer unique
        # And dcasting then resorts to meaningless length aggregation
    longdt <- object %>% melt.data.table(id.vars = 'feature_id')
    longdt[, statistic    := split_extract_fixed(variable, sep, 1) %>% factor(unique(.))]
    longdt[,  coefficient := split_extract_fixed(variable, sep, 2) %>% factor(unique(.))]
    longdt[,       fit    := split_extract_fixed(variable, sep, 3) %>% factor(unique(.))]
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
#' object %<>% lm.()
#' object %<>% limma(block = 'Subject')
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


