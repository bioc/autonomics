#' Across/Within/Between model
#' @param object    SummarizedExperiment
#' @param engine 'limma', 'lm', 'lme', or 'lmer'
#' @param modelvars svars
#' @param across    TRUE/FALSE: fit across  model (additive)    ?
#' @param within    TRUE/FALSE: fit within  model (nested)      ?
#' @param between   TRUE/FALSE: fit between model (interaction) ?
#' @examples
#' object <- survobj()
#' svars(object)
#' awbmod(object, engine = 'limma', modelvars = c('age', 'sex'))
#' awbmod_limma(object, modelvars = c('age', 'sex'), block = 'replicate')
#' awbmod_lm(   object, modelvars = c('age', 'sex'), block = 'replicate')
#' awbmod_lme(  object, modelvars = c('age', 'sex'), block = 'replicate')
#' awbmod_lmer( object, modelvars = c('age', 'sex'), block = 'replicate')
awbmod <- function(
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
    modelfun <- switch(engine, limma = linmod_limma, lm = linmod_lm, lme = linmod_lme, lmer = linmod_lmer)
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


#' @rdname awbmod
#' @export
awbmod_limma <- function(object, ...)  awbmod(object, engine = 'limma', ...)


#' @rdname awbmod
#' @export
awbmod_lm    <- function(object, ...)  awbmod(object, engine = 'lm', ...)


#' @rdname awbmod
#' @export
awbmod_lme   <- function(object, ...)  awbmod(object, engine = 'lme', ...)


#' @rdname awbmod
#' @export
awbmod_lmer  <- function(object, ...)  awbmod(object, engine = 'lmer', ...)


  
