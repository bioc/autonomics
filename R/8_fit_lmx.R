#==============================================================================
#
#                       fit_lmx: lm / lme / lmer
#
#==============================================================================

.cbindstats <- function(
    fitres, fname = 'F value', f.pname = 'Pr(>F)', pname = 'Pr(>|t|)', 
    tname = 't value', effectname = 'Estimate', sename = 'Std. Error'
){
    . <- NULL
    fval <- stats::anova(fitres)['subgroup', fname]
    f.p <- stats::anova(fitres)['subgroup', f.pname]
    fitres %<>% summary()
    fitres %<>% stats::coefficients()
    fitres
    
    pvalues <- as.data.table(as.list(fitres[, pname,      drop = FALSE]))
    tvalues <- as.data.table(as.list(fitres[, tname,      drop = FALSE]))
    effects <- as.data.table(as.list(fitres[, effectname, drop = FALSE]))
    stderrs <- as.data.table(as.list(fitres[, sename,     drop = FALSE]))
    names(pvalues) %<>% paste0('p.', .)
    names(tvalues) %<>% paste0('t.', .)
    names(effects) %<>% paste0('effect.', .)
    names(stderrs) %<>% paste0('se.', .)
    cbind(pvalues, tvalues, effects, stderrs, F = fval, F.p = f.p )
}

.lm <- function(sd, formula, block, weights, statvars, sep, optim = NULL){
    # Initialize
        value <- NULL
        formula <- as.formula(formula)
        environment(formula) <- environment()
    # Run mock lm on zero-imputed data to get all potential coefnames
    # before any get dropped because of all values in a block being NA (for this feature)
        sd0 <- copy(sd)
        sd0[is.na(value), value := 0]  
        fitres0 <- lm( formula = formula, data = sd0, weights = weights, na.action = stats::na.omit )
        fitres0 %<>% summary()
        fitres0 %<>% stats::coefficients()
        fitres0[] <- NA_real_              # Set coefvalues to NA, this is mock data only
    # Run actual lm on actual data
    # Rbind missing coefficients from mock lm
        fitres <- lm( formula = formula, data = sd,  weights = weights, na.action = stats::na.omit )
        Fp <- suppressWarnings(stats::anova(fitres))  # ANOVA F-tests on an essentially perfect fit are unreliable
        Fp <- Fp %>% extract(-nrow(.), , drop = FALSE)
        Fp <- Fp[, 'Pr(>F)'] %>% set_names(rownames(Fp))
        names(Fp) %<>% paste0('p~', .)
        fitres %<>% summary()                      # weights: stackoverflow.com/questions/51142338
        fitres %<>% stats::coefficients()
        rows <- setdiff(rownames(fitres0), rownames(fitres))
        if (!is.null(rows)){  fitres %<>% rbind( fitres0[ rows , , drop = FALSE]  ) 
                              fitres %<>% extract(rownames(fitres0), )  }
    # Reformat and Return
        colnames(fitres) %<>% stri_replace_first_fixed('Estimate', 'effect')
        colnames(fitres) %<>% stri_replace_first_fixed('Std. Error', 'se')
        colnames(fitres) %<>% stri_replace_first_fixed('t value',  't')
        colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|t|)', 'p')
        fitres %<>% extract(, statvars, drop = FALSE)
        fitmat <- matrix(fitres, nrow = 1)
        colnames(fitmat) <- paste(rep(colnames(fitres), each = nrow(fitres)), 
                                  rep(rownames(fitres), times = ncol(fitres)), sep = sep)
        data.table(cbind(fitmat , t(Fp)))
}

.lme <- function(sd, formula, block, weights, statvars, sep, opt = 'optim'){
    ctrl <- nlme::lmeControl(opt = opt)  # https://stats.stackexchange.com/a/40664
    fitres <- nlme::lme( fixed = formula, 
                        random = block, 
                          data = sd,
                     na.action = stats::na.omit, 
                       control = ctrl )
    Fp <- suppressWarnings(stats::anova(fitres)[-1, , drop = FALSE])
    Fp <- Fp[, 'p-value'] %>% set_names(rownames(Fp))
    names(Fp) %<>% paste0('p~', .)
    suppressWarnings(fitres %<>% summary())  # only 2 replicates in a group -> df = 0 -> p = NaN -> warning
    fitres %<>% stats::coefficients()
    colnames(fitres) %<>% stri_replace_first_fixed('Value', 'effect')
    colnames(fitres) %<>% stri_replace_first_fixed('Std.Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t-value', 't')
    colnames(fitres) %<>% stri_replace_first_fixed('p-value', 'p')
    fitres %<>% extract(, statvars, drop = FALSE)
    fitmat <- matrix(fitres, nrow = 1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each = nrow(fitres)), 
                        rep(rownames(fitres), times = ncol(fitres)), sep = sep )
    data.table(cbind(fitmat, t(Fp)))
}

.lmer <- function(sd, formula, block = NULL, weights, statvars, sep, optim = NULL){
    fitres <- lme4::lmer(  formula = formula,
                              data = sd,
                           weights = weights,
                         na.action = stats::na.omit, # https://stackoverflow.com/a/55367171
                           control = lme4::lmerControl(    check.conv.grad = lme4::.makeCC(action = 'ignore', tol=2e-3 ),
                                                       check.conv.singular = lme4::.makeCC(action = "ignore", tol=1e-4 ),
                                                           check.conv.hess = lme4::.makeCC(action = 'ignore', tol=1e-6 )))
    fitres %<>% lmerTest::as_lmerModLmerTest()
    Fp <- suppressWarnings(stats::anova(fitres))
    Fp <- Fp[, 'Pr(>F)'] %>% set_names(rownames(Fp))
    names(Fp) %<>% paste0('p~', .)
    fitres %<>% summary() %>% stats::coefficients()
    colnames(fitres) %<>% stri_replace_first_fixed('Estimate', 'effect')
    colnames(fitres) %<>% stri_replace_first_fixed('Std. Error', 'se')
    colnames(fitres) %<>% stri_replace_first_fixed('t value',  't')
    colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|t|)', 'p')
    fitres %<>% extract(, statvars, drop = FALSE)
    fitmat <- matrix(fitres, nrow=1)
    colnames(fitmat) <- paste(rep(colnames(fitres), each = nrow(fitres)), 
                        rep(rownames(fitres), times = ncol(fitres)), sep = sep )
    data.table(cbind(fitmat, t(Fp)))
}

.extractstat <- function(fitres, quantity){
    idx <- stri_startswith_fixed(names(fitres), paste0(quantity, FITSEP))
    mat <- as.matrix(fitres[, idx, with=FALSE])
    rownames(mat) <- fitres$feature_id
    colnames(mat) %<>% stri_replace_first_fixed(paste0(quantity, FITSEP), '')
    mat
}

#' Statistical models supported in autonomics
#' @examples
#' TESTS
#' @export
TESTS <- c('limma','lm','lme','lmer', 'wilcoxon')

addlhs  <- function(formula)  as.formula(paste0('value ', formula2str(formula)))
droplhs <- function(formula)  as.formula(stri_replace_first_regex(
    formula2str(formula), '^value[ ]*', ''))



#' Put block in lme-compatible format
#' @param formula  formula
#' @param block    block: charactervector or formula
#' @param verbose  TRUE or FALSE
#' @param ...      required for s3 dispatch
#' @examples
#' # lme: ensure lme-compatiblae block specification
#'     block2lme( block = list(subject = ~1, batch = ~1))
#'     block2lme( block =   ~1|subject)
#'     block2lme( block =   c('subject',    'batch'))
#' 
#' # lm: integrate block into formula as random effect
#'     formula2lm(   formula = ~ subgroup,  block = c('subject', 'batch') )
#' 
#' # lmer: integrate block into formula as fixed effect
#'     formula2lmer( formula = ~ subgroup,  block = c('subject',    'batch') )
#'     formula2lmer( formula = ~ subgroup         + (1|subject) + (1|batch ) )
#' @export
block2lme <- function(block, ...)  UseMethod('block2lme')

#' @rdname block2lme
#' @export
block2lme.list <- function(block, verbose = TRUE, ...)  return(block)

#' @rdname block2lme
#' @export
block2lme.formula <- function(block, verbose = TRUE, ...){
            block0 <- block
            block <- formula2str(block0)
            block %<>% substr(2,nchar(.)) %>% trimws()  # rm ~
            block %<>% stri_split_fixed('+') %>% unlist() %>% trimws()
       blocknames <- trimws(split_extract_fixed( block, '|', 2 ))
            block <- trimws(split_extract_fixed( block, '|', 1 ))
            block %<>% paste0('~', .)
            block %<>% lapply(as.formula)
            block %<>% set_names(blocknames)
     return(block)
}


#' @rdname block2lme
#' @export
block2lme.character <- function(block, verbose = TRUE, ...){
    block0 <- block
    block <- rep('~1', length(block0))
    names(block) <- block0
    block %<>% lapply(as.formula)
    block
}


#' @rdname block2lme
#' @export
formula2lmer <- function(formula, block){
    if (stri_detect_fixed(formula2str(formula), '|'))  return(formula)
    formula %<>% formula2str()
    block %<>% paste0('(1|', ., ')')
    block %<>% paste0(collapse = ' + ')
    formula %<>% paste0(' + ', block)
    formula %<>% as.formula()
    formula
}


#' @rdname block2lme
#' @export
formula2lm <- function(formula, block){
    if (is.null(block))  return(formula)
    formula %<>% formula2str()
    formula %<>% substr(2,nchar(.))
    formula %<>% trimws()
    formula %<>% c(block, .) %>% paste0(collapse = '+')
    formula %<>% paste0('~', .)
    formula %<>% as.formula()
    formula
}


#' @rdname block2lme
#' @export
block_vars <- function(formula){
    formula %<>% formula2str()
    formula %<>% stri_replace_first_regex('^[~][ ]*', '')
    formula %<>% stri_split_regex('[ ]*[+][ ]*')
    formula %<>% extract2(1)
    formula %<>% extract(stri_detect_fixed(., '|'))
    formula %<>% stri_replace_first_regex('[ ]*[(][ ]*', '') # lmer
    formula %<>% stri_replace_first_regex('[ ]*[)][ ]*', '') # lmer
    formula %<>% split_extract_regex('[ ]*[|][ ]*', 2)
    formula
}


fit_lmx <- function(
       object, 
          fit, 
      formula = as.formula('~ subgroup'),
         drop = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment.explicit,
        coefs = contrast_coefs(object, formula = formula, drop = drop, codingfun = codingfun),
        block = NULL, 
          opt = 'optim',
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
     statvars = c('effect', 'p', 'se', 't')[1:2],
        ftest = if (is.null(coefs))  TRUE else FALSE,
          sep = FITSEP,
       suffix = paste0(sep, fit),
      verbose = TRUE, 
         plot = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(fit, c('lm', 'lme', 'lmer'))
    assert_valid_formula(formula, object)
    assert_is_a_bool(drop)
    if (!is.null(weightvar)){   assert_is_character(weightvar)
                                assert_is_subset(weightvar, assayNames(object)) 
                                message('\t\t\tweights = assays(object)$', weightvar)  }
    assert_is_subset(statvars, c('effect', 'se', 't', 'p'))
    N <- value <- V1 <- NULL
# Filter / Customize
    obj <- object
    if (verbose)  cmessage('%sFilter', spaces(14))
    obj %<>% keep_replicated_features( formula, verbose = verbose)
    obj %<>% keep_connected_blocks(    block,   verbose = verbose)  # keep samples from fully connected blocks (in sdt, feature-specific NA values not considered)
    obj %<>% keep_connected_features(  block,   verbose = verbose)  # keep features with 2+ connected blocks
# Prepare
    object %<>% reset_fit(fit)
    if ( fit == 'lme'  ){     block %<>% block2lme(); mdlvars <-  unique(c(all.vars(formula), names(block)))   }
    if ( fit == 'lmer' ){   formula %<>% formula2lmer(block); mdlvars <- all.vars(formula)                     }
    if ( fit == 'lm'   ){   formula %<>% formula2lm(  block); mdlvars <- all.vars(formula)                     }
    fstr <- formula2str(formula)
    if (verbose)                  cmessage('%s%s',           spaces(14), R.utils::capitalize(fit))
    if (verbose & fit == 'lme' )  cmessage('%s%s, random = %s', spaces(22), fstr, capture.output(dput(block)))
    if (verbose & fit == 'lmer')  cmessage('%s%s',           spaces(21), fstr)
    if (verbose & fit == 'lm'  )  cmessage('%s%s',           spaces(21), fstr)
    fitmethod <- get(paste0('.', fit))
    if (is.null(weightvar)){ weightvar <- 'weights'; weights <- NULL }
    assays <- assayNames(object) %>% intersect(c(.[1], 'weights'))
    dt <- sumexp_to_longdt(obj, svars = mdlvars, assay = assays)
    lhsformula <- addlhs(formula)
# Fit
    fitres <- dt[, fitmethod( .SD,   formula = lhsformula, 
                                       block = block, 
                                     weights = get(weightvar),
                                    statvars = statvars, 
                                         sep = sep,
                                         opt = opt ),            by = 'feature_id' ]
    names(fitres) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    vars <- all.vars(formula)
    if (drop)   for (var in vars){     # t.p: p~subgroupt1 -> p~t1
                    pat <- sprintf('%s(.+)', var)   # f.p: p~subgroup
                    names(fitres) %<>% stri_replace_first_regex(pat, '$1') }
# Extract
                             pattern <- 'feature_id'
    if (!is.null( coefs))    pattern %<>% paste0('|', paste0(coefs, collapse = '|'))
    if (ftest)               pattern %<>% paste0('|', paste0( vars, collapse = '|'))
                             pattern %<>% paste0('(', ., ')$')
    fitres <- fitres[, .SD, .SDcols = patterns(pattern) ]
    names(fitres)[-1] %<>% paste0(suffix)
    if (verbose)  message_df('                      %s', summarize_fit(fitres, fit = fit, coefs = coefs))
# Merge back
    object %<>% merge_fit(fitres)
    formula %<>% droplhs() %<>% formula2str()
    
    if (!is.null(weights))  formula %<>% paste0(', weights = assays(object)$', weightvar)
    if (length(coefs) > 1)  coefs %<>% setdiff('Intercept')
    if (plot)  print(plot_volcano(object, fit = fit, coefs = coefs))
    object 
}


#' @rdname fit
#' @export
fit_lm <- function(
       object,
      formula = as.formula('~ subgroup'),
         drop = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment.explicit,
       design = NULL,  # only to make fit(.) work!
        block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
     statvars = c('effect', 'p', 'se', 't')[1:2],
          sep = FITSEP,
       suffix = paste0(sep, 'lm'),
        coefs = contrast_coefs(object, formula = formula, drop = drop, codingfun = codingfun), 
    contrasts = NULL,
        ftest = if (is.null(coefs)) TRUE else FALSE,
      verbose = TRUE
){
    
    sdt(object) %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
    fit_lmx(    object,
                   fit = 'lm', 
               formula = formula,
                  drop = drop,
             codingfun = codingfun,
                 block = block,
             weightvar = weightvar,
              statvars = statvars,
                 ftest = ftest,
                   sep = sep,
                suffix = suffix,
                 coefs = coefs,
               verbose = verbose )
}


#' @rdname fit
#' @export
fit_lme <- function(
       object, 
      formula = as.formula('~ subgroup'),
         drop = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment.explicit,
       design = NULL,  # only to make fit(.) work!
        block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
          opt = 'optim',
     statvars = c('effect', 'p', 'se', 't')[1:2],
          sep = FITSEP,
       suffix = paste0(sep, 'lme'),
        coefs = contrast_coefs(object, formula = formula, drop = drop, codingfun = codingfun), 
    contrasts = NULL,
        ftest = if (is.null(coefs))  TRUE else FALSE,
      verbose = TRUE
){
# Assert
    . <- NULL
    if (!requireNamespace('nlme', quietly = TRUE)){
        message("BiocManager::install('nlme'). Then re-run.")
        return(object)   }
# Fit
    sdt(object) %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
    fit_lmx(    object,
                   fit = 'lme', 
               formula = formula,
                  drop = drop,
             codingfun = codingfun,
                 block = block, 
             weightvar = weightvar,
              statvars = statvars,
                 ftest = ftest,
                   sep = sep,
                suffix = suffix,
                   opt = opt,
                 coefs = coefs, 
               verbose = verbose )
}


#' @rdname fit
#' @export
fit_lmer <- function(
       object, 
      formula = as.formula('~ subgroup'),
         drop = varlevels_dont_clash(object, all.vars(formula)),
    codingfun = contr.treatment.explicit,
       design = NULL,  # only to make fit(.) work!
        block = NULL, 
    weightvar = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
     statvars = c('effect', 'p', 'se', 't')[1:2],
          sep = FITSEP,
       suffix = paste0(sep, 'lmer'),
        coefs = contrast_coefs(object, formula = formula, drop = drop, codingfun = codingfun), 
    contrasts = NULL,
        ftest = if (is.null(coefs)) TRUE else FALSE,
      verbose = TRUE
){
# Assert
    . <- NULL
    if (!requireNamespace('lme4', quietly = TRUE)){
        message("`BiocManager::install('lme4')`. Then re-run.")
        return(object) }
    if (!requireNamespace('lmerTest', quietly = TRUE)){
        message("`BiocManager::install('lmerTest')`. Then re-run.")
        return(object) }
# Fit
    sdt(object) %<>% code(codingfun = codingfun, vars = all.vars(formula), verbose = verbose)
    fit_lmx(    object,
                   fit = 'lmer', 
               formula = formula,
                  drop = drop,
             codingfun = codingfun,
                 block = block, 
             weightvar = weightvar,
              statvars = statvars,
                 ftest = ftest,
                   sep = sep,
                suffix = suffix,
                 coefs = coefs, 
               verbose = verbose )
}
