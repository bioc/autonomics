#' General Linear Modeling (across-within-between interface)
#' @param object    SummarizedExperiment
#' @param engine 'limma', 'lm', 'lme', or 'lmer'
#' @param modelvars svars
#' @param across    TRUE/FALSE: fit across  model (additive)    ?
#' @param within    TRUE/FALSE: fit within  model (nested)      ?
#' @param between   TRUE/FALSE: fit between model (interaction) ?
#' @param coding    character: codingfunname
#' @param drop      TRUE or FALSE
#' @param ...       passed to linmod
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' svars(object)
#' awblinmod_limma(object, modelvars = c('Diabetes', 'Time'), block = 'Subject')
#' awblinmod_lme(  object, modelvars = c('Diabetes', 'Time'), block = 'Subject')
#' awblinmod_lmer( object, modelvars = c('Diabetes', 'Time'), block = 'Subject')
#' awblinmod_lm(   object, modelvars = c('Diabetes', 'Time'))
#' awblinmod(object, engine = 'limma', modelvars = 'Time')
#' awblinmod(object, engine = 'limma', modelvars = c('Diabetes', 'Time'))
awblinmod <- function(
    object, 
    engine,
 modelvars,
    across = TRUE,
    within = if (length(modelvars)==1) FALSE else TRUE,
   between = if (length(modelvars)==1) FALSE else TRUE,
    coding = c('code_control', 'code_diff'),
      drop = TRUE, 
          ...
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(modelvars, svars(object))

# Model
    modelfun <- get(sprintf('linmod_%s', engine))
    if (across){
        formula  <- paste0(modelvars, collapse = '+')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        for (codi in coding){
            coefs  <- contrast_coefs(object,  formula, coding = codi, drop = drop)
            coefs %<>% setdiff(autonomics::coefs(object, fit = engine))
            object %<>% modelfun(formula,  coding = codi, drop = drop, coefs =  coefs, verbose = TRUE, reset = FALSE, ...)
        }
    }
    if (within){
        formula <- paste0(modelvars, collapse = '/')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        for (codi in coding){
            coefs  <- contrast_coefs(object, formula, coding = codi, drop = drop)
            coefs %<>% extract(stri_detect_fixed(., ':'))
            coefs %<>% setdiff(autonomics::coefs(object, fit = engine))
            object %<>% modelfun(formula, coding = codi, drop = drop, coefs = coefs, verbose = TRUE, reset = FALSE, ...)
        }
        fvars(object) %<>% stri_replace_first_fixed(':', '/') # needs to be out of the loop !
    }
    if (within){
        formula <- paste0(rev(modelvars), collapse = '/')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        for (codi in coding){
            coefs  <- contrast_coefs(object, formula, coding = codi, drop = drop)
            coefs %<>% extract(stri_detect_fixed(., ':'))
            coefs %<>% setdiff(autonomics::coefs(object, fit = engine))
            object %<>% modelfun(formula, coding = codi, drop = drop, coefs = coefs, verbose = TRUE, reset = FALSE, ...)
        }
        fvars(object) %<>% stri_replace_first_fixed(':', '/')
    }
    if (between){
        formula <- paste0(modelvars, collapse = '*')
        formula %<>% paste0('~', .)
        formula %<>% as.formula()
        for (codi in coding){
            coefs  <- contrast_coefs(object, formula, coding = codi, drop = drop)
            coefs %<>% extract(stri_detect_fixed(., ':'))
            object %<>% modelfun(formula, coding = codi, drop = drop, coefs = coefs, verbose = TRUE, reset = FALSE, ...)
        }
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


  
