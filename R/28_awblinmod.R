#' General Linear Modeling (across-within-between interface)
#' @param object    SummarizedExperiment
#' @param engine 'limma', 'lm', 'lme', or 'lmer'
#' @param modelvars svars
#' @param across    TRUE/FALSE: fit across  model (additive)    ?
#' @param within    TRUE/FALSE: fit within  model (nested)      ?
#' @param between   TRUE/FALSE: fit between model (interaction) ?
#' @examples
#' object <- survobj()
#' svars(object)
#' awblinmod(object, 'limma', modelvars = 'age')
#' 
#' awblinmod(object, engine = 'limma', modelvars = c('age', 'sex'))
#' awblinmod_limma(object, modelvars = c('age', 'sex'), block = 'replicate')
#' awblinmod_lm(   object, modelvars = c('age', 'sex'), block = 'replicate')
#' awblinmod_lme(  object, modelvars = c('age', 'sex'), block = 'replicate')
#' awblinmod_lmer( object, modelvars = c('age', 'sex'), block = 'replicate')
awblinmod <- function(
    object, 
    engine,
 modelvars,
    across = TRUE,
    within = if (length(modelvars)==1) FALSE else TRUE,
   between = if (length(modelvars)==1) FALSE else TRUE,
 codingfun = code_control,
      drop = TRUE, 
          ...
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(modelvars, svars(object))

# Model
    modelfun <- switch(engine, limma = linmod_limma, lm = linmod_lm, lme = linmod_lme, lmer = linmod_lmer)
    if (across){
        formula  <- paste0(modelvars, collapse = '+')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs  <- colnames(create_design(object,  formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% setdiff('Intercept')
        object %<>% modelfun(formula,  codingfun = codingfun, drop = drop, coefs =  coefs, verbose = FALSE, ...)
    }
    if (within){
        formula <- paste0(modelvars, collapse = '/')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs <- colnames(create_design(object, formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% extract(stri_detect_fixed(., ':'))
        object %<>% modelfun(formula, codingfun = codingfun, drop = drop, coefs = coefs, verbose = FALSE, ...)
        fvars(object) %<>% stri_replace_first_fixed(':', '/')
    }
    if (within){
        formula <- paste0(rev(modelvars), collapse = '/')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs <- colnames(create_design(object, formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% extract(stri_detect_fixed(., ':'))
        object %<>% modelfun(formula, codingfun = codingfun, drop = drop, coefs = coefs, verbose = FALSE, ...)
        fvars(object) %<>% stri_replace_first_fixed(':', '/')
    }
    if (between){
        formula <- paste0(modelvars, collapse = '*')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        coefs <- colnames(create_design(object, formula, codingfun = codingfun, drop = drop, verbose = FALSE))
        coefs %<>% extract(stri_detect_fixed(., ':'))
        object %<>% modelfun(formula, codingfun = codingfun, drop = drop, coefs = coefs, verbose = FALSE, ...)
        fvars(object) %<>% stri_replace_first_fixed(':', '*')
    }
# Return
    object
}


#' @rdname awblinmod
#' @export
awblinmod_limma <- function(object, ...)  awblinmod(object, engine = 'limma', ...)


#' @rdname awblinmod
#' @export
awblinmod_lm    <- function(object, ...)  awblinmod(object, engine = 'lm', ...)


#' @rdname awblinmod
#' @export
awblinmod_lme   <- function(object, ...)  awblinmod(object, engine = 'lme', ...)


#' @rdname awblinmod
#' @export
awblinmod_lmer  <- function(object, ...)  awblinmod(object, engine = 'lmer', ...)


  
