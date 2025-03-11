

#' Survival analysis example
#' @return SummarizedExperiment
#' @examples
#' #  Background    BOC   CD4   PLG   XCL1
#' #                ===   ===   ===   ====
#' #                Flo   Fhi   Flo   Fhi    _______________
#' #                Fhi   Flo   Fhi   Flo    ____________  |
#' #                Mlo   Mhi   Mhi   Mlo    _____      |  |_____
#' #                Mhi   Mlo   Mlo   Mhi    __   |__   |__     |__
#' #                                           |     |     |
#' #                                           |__   |__   |_______
#' #                                              |     |
#' #                                              |__    |__
#' #                                                 |      |
#' survex()
#' @export 
survex <- function(){
    
    set.seed(1)
    expr1.m <- rnorm(10,3)
    expr1.f <- rnorm(10,3)
    expr2.m <- rnorm(10,5)
    expr2.f <- rnorm(10,5)# m.surv1  m.surv2  f.surv3  f.surv4
    mat <- rbind( geneA = c(expr1.m, expr2.m, expr1.f, expr2.f), 
                  geneB = c(expr2.m, expr1.m, expr2.f, expr1.f),
                  geneC = c(expr1.m, expr2.m, expr2.f, expr1.f), 
                  geneD = c(expr2.m, expr1.m, expr1.f, expr2.f))
    object <- SummarizedExperiment::SummarizedExperiment(list(exprs = mat))
    fdt(object)$feature_id <- fnames(object)
    object$sample_id <- snames(object)  <- c(sprintf('surv1.m.%d', 0:9),  sprintf('surv2.m.%d', 0:9), sprintf('surv3.f.%d', 0:9), sprintf('surv4.f.%d', 0:9))
    object$survgroup <- object$sample_id %>% split_extract_fixed('.', 1)
    object$sex       <- object$sample_id %>% split_extract_fixed('.', 2)
    object$replicate <- object$sample_id %>% split_extract_fixed('.', 3)

    surv1.m.time <- c( rep(1,4), rep(2,4), rep(3,2) ); surv1.m.event <- rep(1,10)
    surv2.m.time <- c( rep(2,2), rep(3,4), rep(4,4) ); surv2.m.event <- rep(1,10)
    surv3.f.time <- c( rep(4,4), rep(5,4), rep(7,2) ); surv3.f.event <- c(rep(1,8),rep(0,2))
    surv4.f.time <- c( rep(5,2), rep(6,2), rep(7,6) ); surv4.f.event <- c(rep(1,4),rep(0,6))
    object$timetoevent <- c(surv1.m.time,  surv2.m.time,  surv3.f.time,  surv4.f.time)
    object$event       <- c(surv1.m.event, surv2.m.event, surv3.f.event, surv4.f.event)
    object
    
}



#' @rdname dot-fit_survival
#' @export
survival_example <- function(){
    sampledt <- rbind(  data.table( subgroup = 'Control',  sample_id = 'C01', timetoevent = 3, event = 1), 
                        data.table( subgroup = 'Control',  sample_id = 'C02', timetoevent = 3, event = 1), 
                        data.table( subgroup = 'Control',  sample_id = 'C03', timetoevent = 4, event = 1),
                        data.table( subgroup = 'Control',  sample_id = 'C04', timetoevent = 4, event = 1),
                        data.table( subgroup = 'Control',  sample_id = 'C05', timetoevent = 4, event = 0),
                        data.table( subgroup = 'Control',  sample_id = 'C06', timetoevent = 4, event = 0),
                        data.table( subgroup = 'Control',  sample_id = 'C07', timetoevent = 4, event = 0), 
                        data.table( subgroup = 'Control',  sample_id = 'C08', timetoevent = 4, event = 0), 
                        data.table( subgroup = 'Control',  sample_id = 'C09', timetoevent = 4, event = 0),
                        data.table( subgroup = 'Control',  sample_id = 'C10', timetoevent = 4, event = 0),
                        data.table( subgroup = 'Diseased', sample_id = 'D01', timetoevent = 1, event = 1), 
                        data.table( subgroup = 'Diseased', sample_id = 'D02', timetoevent = 1, event = 1), 
                        data.table( subgroup = 'Diseased', sample_id = 'D03', timetoevent = 2, event = 1),
                        data.table( subgroup = 'Diseased', sample_id = 'D04', timetoevent = 2, event = 1),
                        data.table( subgroup = 'Diseased', sample_id = 'D05', timetoevent = 2, event = 0), # lets include right censoring examples too !
                        data.table( subgroup = 'Diseased', sample_id = 'D06', timetoevent = 2, event = 0),
                        data.table( subgroup = 'Diseased', sample_id = 'D07', timetoevent = 3, event = 1), 
                        data.table( subgroup = 'Diseased', sample_id = 'D08', timetoevent = 3, event = 1), 
                        data.table( subgroup = 'Diseased', sample_id = 'D09', timetoevent = 4, event = 1), 
                        data.table( subgroup = 'Diseased', sample_id = 'D10', timetoevent = 4, event = 1) )
    n <- nrow(sampledt)/2
    object <- rbind( ASGR1 = c(rnorm(n, mean =  2), rnorm(n, mean =  5)),
                       BOC = c(rnorm(n, mean =  2), rnorm(n, mean =  6)),
                       CD4 = c(rnorm(n, mean =  2), rnorm(n, mean =  7)),
                      LY86 = c(rnorm(n, mean =  2), rnorm(n, mean =  8)),
                       CFI = c(rnorm(n, mean =  6), rnorm(n, mean =  3)),
                       PLG = c(rnorm(n, mean =  7), rnorm(n, mean =  3)),
                      PROC = c(rnorm(n, mean =  8), rnorm(n, mean =  3)),
                      XCL1 = c(rnorm(n, mean =  9), rnorm(n, mean =  3)) )
    colnames(object) <- sampledt$sample_id
    object <- SummarizedExperiment(list(exprs = object))
    fdt(object) <- data.table(feature_id = fnames(object))
    sdt(object) <- sampledt
    object
}


#' Fit onefeature survival 
#' @param timetoevent  numeric (time to event)
#' @param event        numeric (1=event, 0=not)
#' @param xvalues      numeric (.coxph) or twolevel-factor (.survdiff, .logrank_test)
#' @param xname        string: used in fvar
#' @param drop         TRUE or FALSE : drop xname in output fvar ?
#' @examples
#' # Prepare
#'      sd <- survdt()
#'      sd[ , quantile := factor(dplyr::ntile(BOC, 2)) ]
#' # Survival
#'        .coxph(sd, Surv(timetoevent, event) ~ quantile)
#'     .survdiff(sd, Surv(timetoevent, event) ~ quantile)
#'      .logrank(sd, Surv(timetoevent, event) ~ quantile)
#' @rdname dot-coxph
#' @export
.coxph <- function(sd, formula){
    fitres <- survival::coxph(formula = formula, data = sd)
    #Fres <- suppressWarnings(stats::anova(fitres))
    #Fres <- Fres %>% extract(-1, , drop = FALSE)
    #pF <- Fres[, 'Pr(>|Chi|)' ] %>% set_names(paste0('PF~', rownames(Fres)))
    #tF <- Fres[, 'Chisq'      ] %>% set_names(paste0('F~', rownames(Fres)))
    fitres %<>% summary()
    fitres %<>% stats::coefficients()
    colnames(fitres) %<>% stri_replace_first_fixed('se(coef)', 'se') # dont reverse order of these two lines
    colnames(fitres) %<>% stri_replace_first_fixed('coef', 'effect')
    colnames(fitres) %<>% stri_replace_first_fixed('Pr(>|z|)', 'p')  # dont reverse order of these two lines
    colnames(fitres) %<>% stri_replace_first_fixed('z', 't')
    fitres %<>% extract(, c('effect', 't', 'p'), drop = FALSE)
    fitmat <- matrix(fitres, nrow = 1)
    colnames(fitmat) <- paste(rep(colnames(fitres),  each = nrow(fitres)), 
                              rep(rownames(fitres), times = ncol(fitres)), sep = '~')
    data.table(fitmat)
    #data.table(cbind(fitmat, t(tF), t(pF)))
}


#' @rdname dot-coxph
#' @export
.survdiff <- function(sd, formula){
    xvar <- labels(terms(formula))
    xvalues <- sd[[xvar]]
    xlevel1 <-     levels(xvalues)[1]
    xleveln <- rev(levels(xvalues))[1]
    
    sd %<>% extract(get(xvar) %in% c(xlevel1, xleveln))
    survout <- suppressWarnings(survival::survdiff(formula = formula, data = sd))
   meandiff <- sd[ , mean(timetoevent[get(xvar)==xleveln]) -
                     mean(timetoevent[get(xvar)==xlevel1]) ]
    outdt <- data.table(  effect = -meandiff,
                               t = -sign(meandiff) * survout$chisq,
                               p =  1 - pchisq(survout$chisq, 1)  )
    newnames <- sprintf('%s~%s%s-%s~survdiff', names(outdt), xvar, xleveln, xlevel1)
    setnames(outdt, names(outdt), newnames)
    outdt[]
}


#' @rdname dot-coxph
#' @export
.logrank <- function(sd, formula){
    xvar <- labels(terms(formula))
    xvalues <- sd[[xvar]]                                #   NOTE  The coin statistic is signed for twogroup comparisons 
    xlevel1 <-     levels(xvalues )[1]                        #   But unsigned for multigroup comparisons (which are anova like)
    xleveln <- rev(levels(xvalues))[1]                        #   But unsigned for multigroup comparisons (which are anova like)
            

    sd %<>% extract(get(xvar) %in% c(xlevel1, xleveln))         
    survout <- suppressWarnings(coin::logrank_test(formula = formula, data = sd))
  meandiff <- sd[ , mean(timetoevent[get(xvar)==xleveln]) -
                    mean(timetoevent[get(xvar)==xlevel1]) ]
    outdt <- data.table( effect = -meandiff,
                              t = -sign(meandiff) * abs(coin::statistic(survout)),
                              p = coin::pvalue(survout) )
    newnames <- sprintf('%s~%s%s-%s~logrank', names(outdt), xvar, xleveln, xlevel1)
    setnames(outdt, names(outdt), newnames )
    outdt[]
}



#' Bin/Factorize assay
#' @param object  SummarizedExperiment
#' @param assay   string
#' @param k       number of bins/levels
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' object <- survex()
#'       bin_assay(object, k = 4)
#' factorize_assay(object, k = 4)
#' @export
bin_assay <- function(object, assay = assayNames(object)[1], k = 3, verbose = TRUE){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, assayNames(object))
    assert_is_a_number(k)
    assert_is_a_bool(verbose)
# Bin
    mat <- assays(object)[[assay]]
    mat %<>% apply(1, dplyr::ntile, n = k) %>% t()
    colnames(mat) <- colnames(object)
# Add
    newassayname <- sprintf('%s%dbins', assay, k)
    if (verbose)   cmessage('%sAdd `%s`', spaces(8), newassayname)
    assays(object)[[newassayname]] <- mat
    object
}


#' @rdname bin_assay
#' @export
factorize_assay <- function(object, assay = assayNames(object)[1], k = 3, verbose = TRUE){
# Bin (assertions done during binning)
    object %<>% bin_assay(assay = assay, k = k, verbose = verbose)
# Factorize
    binnedassay <- sprintf('%s%dbins', assay, k)
    mat <- assays(object)[[binnedassay]]
    mode(mat) <- 'character'
    mat %<>% paste0('xpr', .)
    dim(mat) <- dim(object)
    dimnames(mat) <- dimnames(object)
# Add
    newassayname <- sprintf('%s%dlevels', assay, k)
    if (verbose)   cmessage('%sAdd `%s`', spaces(8), newassayname)
    assays(object)[[newassayname]] <- mat
    object
}


#' Fit survival
#' 
#' Compute survival effect of svars, exprs, and their interactions
#' 
#' @param object    SummarizedExperiment
#' @param formula   Formula
#' @param bins      Number of value bins. Zero means unbinned.
#' @param bintype  'factor' or 'numeric'
#' @param engine   'coxph', 'survdiff', or 'logrank'
#' @param drop      Whether to drop factor varname in coefnames
#' @param codingfun (factor) coding function
#' @param verbose   TRUE or FALSE
#' @examples
#' # Load/Transform
#'    object <- survex()
#'    object %<>% bin_assay(k = 2)
#'    object %<>% factorize_assay(k = 2)
#' # coxph{survival}
#'   .fit_survival(object)
#'   .fit_survival(object, ~ exprs)                         #      expr effect
#'   .fit_survival(object, ~ exprs2bins)                    #   exprbin effect
#'   .fit_survival(object, ~ exprs2levels)                  # exprlevel effect
#'   .fit_survival(object, formula = ~ sex)                 #       sex effect
#'   .fit_survival(object, formula = ~ sex + exprs2levels)  #       sex effect ACROSS exprlevels,  exprlevel effect ACROSS sexes.
#'   .fit_survival(object, formula = ~ sex / exprs2levels)  # exprlevel effect WITHIN sex,               sex effect ACROSS exprlevels.
#'   .fit_survival(object, formula = ~ exprs2levels / sex)  #       sex effect WITHIN exprlevel,   exprlevel effect ACROSS sexes
#'   .fit_survival(object, formula = ~ exprs2levels * sex)  #       sex effect differences BETWEEN exprlevels
#' # survdiff{survival}
#'   .fit_survival(object, formula = ~ exprs2levels, engine = 'survdiff')
#'   .fit_survival(object, formula = ~ exprs2levels, engine = 'logrank')
.fit_survival <- function( 
       object,
       formula = as.formula(sprintf('~%s', assayNames(object)[1])),
        engine = c('coxph', 'survdiff', 'logrank')[1],
          drop = TRUE,
     codingfun = code_control,
       verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(c('timetoevent', 'event'), svars(object))
    assert_is_formula(formula)
    assert_is_subset(all.vars(formula), c(assayNames(object), svars(object)))
    assert_scalar_subset(engine, c('coxph', 'survdiff', 'logrank'))
    assert_is_a_bool(drop)
    assert_is_function(codingfun)
    assert_is_a_bool(verbose)
    if (engine == 'logrank')  if (!installed('coin'))  return(NULL)
    object %<>% filter_samples(!is.na(event) & !is.na(timetoevent))
# Code
    survivalvars <- c('timetoevent', 'event')
      samplevars <- intersect(all.vars(formula),      svars(object))
        assayvar <- intersect(all.vars(formula), assayNames(object))
    if (is_empty(assayvar)){   
        dt <- sdt(object)[, c('sample_id', samplevars, survivalvars), with = FALSE]
        dt[, feature_id := formula2str(formula)] 
    } else {
        assert_is_scalar(assayvar)   # currently written for scalar assayvar
        if (engine == 'survdiff')   assert_character_matrix(assays(object)[[assayvar]], .xname = sprintf('assays(object)$%s', assayvar))
        dt <- sumexp_to_longdt(object, svars = c(samplevars, survivalvars), assay = assayvar, value.name = assayvar) 
        charactercols <- vapply(dt, is.character, logical(1))           # lower-level functions expect factors
        charactercols %<>% extract(. == TRUE)
        charactercols %<>% names()
        for (col in charactercols)   dt[ , (col) := factor(get(col)) ]  # this ensures level order
    }
    dt %<>% code(codingfun = codingfun, vars = c(assayvar, samplevars), verbose = verbose)
# Fit
    twosideformula <- formula
    twosideformula %<>% formula2str()
    twosideformula %<>% paste0('Surv(timetoevent, event)', .)
    if (verbose)  cmessage('%s%s(%s)', spaces(8), engine, twosideformula)
    twosideformula %<>% as.formula()
    if (engine == 'coxph')     fitres <- dt[,    .coxph(.SD, twosideformula), by = 'feature_id']
    if (engine == 'survdiff')  fitres <- dt[, .survdiff(.SD, twosideformula), by = 'feature_id']
    if (engine == 'logrank')   fitres <- dt[,  .logrank(.SD, twosideformula), by = 'feature_id']
    
    if (drop)   for (var in c(assayvar, samplevars)){  
                    pat <- sprintf('%s(.+)', var)
                    names(fitres) %<>% stri_replace_first_regex(pat, '$1')  }
    names(fitres)[-1] %<>% paste0('~coxph')
# Merge    
    if (verbose)  message_df('                      %s', summarize_fit(fitres))
    #if ('expr' %in% all.vars(formula)){  object %<>% merge_fit(fitres)
    #} else {                              metadata(object)$survival <- fitres[, -1] }
    #object
    fitres
}



#' Fit survival 
#' 
#' Compute association between survival and expression (or svar)
#' 
#' Compute association between survival and expression (or svar)
#' \verb{    } Continuous for \code{coxph}.                                         \cr
#' \verb{    } Categorical for \code{survdiff} or \code{logrank}                    \cr
#' \verb{        } Samples are split into \code{ntile} expression groups.           \cr
#' \verb{        } Survival is compared between highest and lowest expressors.      \cr 
#' 
#' Three statistics recorded per engine                                             \cr
#' \verb{        } \code{p}                                                         \cr
#' \verb{   } \code{effect: } coef (\code{coxph})                                   \cr
#' \verb{           } mean survival difference (\code{survdiff, logrank})           \cr
#' \verb{        } \code{t: } \eqn{z}  (\code{coxph})                               \cr
#' \verb{           }         \eqn{ \chi^2} (\code{survdiff}, \code{logrank})       \cr
#' \verb{                } sign reflects whether expression                         \cr
#' \verb{                } increases (positive) or decreases (negative) survival
#' @param object      SummarizedExperiment
#' @param splitvar    svar or assayName
#' @param engine     'coxph' 'survdiff' or 'logrank'
#' @param drop        TRUE or FALSE
#' @param ntile       number
#' @param splitvar       string
#' @param sep         fvar string separator : e.g. '~' gives p~surv~LR50 
#' @param verbose     TRUE or FALSE
#' @param plot        TRUE or FALSE
#' @param width       number
#' @param height      number
#' @param n           number of features to plot
#' @param nrow        number of rows
#' @param ncol        number of cols
#' @param outdir      dir
#' @param writefunname 'write_xl' or 'write_ods'
#' @return SummarizedExperiment
#' @examples                                                 # Innerfun
#'  object <- survival_example()                             #     returns data.table
#'                                                           #     accepts scalar args
#' .fit_survival(object)                                     #         effect of exprquantile (2-1) on survival
#' .fit_survival(object, ntile = 3)                          #         effect of exprquantile (3-1) on survival
#' .fit_survival(object, splitvar = 'subgroup')              #         effect of subgroup (Disease-Control) on survival
#' .fit_survival(object, engine = 'coxph')                   #         effect of exprvalue (continous) on survival
#'                                                           # Outerfun 
#'  fit_survival(object)                                     #     returns SummarizedExperiment
#'                                                           #     accepts vector args
#'  fit_survival(object, ntile = c(2,3))                     #             ntile: multiple contrasts: exprquantiles contrasts: (2-1) and (3-1)
#'  fit_survival(object, splitvar = c('subgroup', 'exprs'))  #          splitvar: multiple exprquantile  (2-1) and subgroups
#'  fit_survival(object, engine = c('survdiff', 'coxph'))    #            engine: survdiff (categorical) and coxph (continous)
#'                                                           #     Writes
#'                                                           #     Plots
#' @export
fit_survival <- function(
        object, 
       formula = as.formula(sprintf('~%s', assayNames(object)[1])),
        engine = c('coxph', 'survdiff', 'logrank')[1],
          drop = TRUE,
     codingfun = code_control,
       verbose = TRUE,
        outdir = NULL,
          plot = if (is.null(outdir)) FALSE else TRUE,
         width = 7,
        height = 7,
             n = min(nrow(object), 9),
          ncol = 3,
          nrow = 3,
  writefunname = 'write_xl'
){
    if (verbose)  cmessage('%sSurvival', spaces(8))
# Compute
    for (eng in engine){
        outdt <- .fit_survival(  object = object, 
                                formula = formula,
                                 engine = engine,
                                   drop = drop,
                              codingfun = codingfun,
                                verbose = verbose )
                 object %<>% merge_fdt(outdt[feature_id != formula2str(formula)])
        metadata(object)$survival  <-  outdt[feature_id == formula2str(formula)]
    }
# Write
    if (!is.null(outdir)){
        outdir <- sprintf('%s/survival', outdir)
        dir.create(outdir, showWarnings = FALSE)
        tableext <- switch(writefunname, write_xl = 'xlsx', write_ods = 'ods')
        tablefile <- if (is.null(outdir)) NULL else sprintf('%s/survival.%s',    outdir, tableext)
        get(writefunname)(object, tablefile)
    }
# Plot
    if (plot){
        file <- if (is.null(outdir)) NULL else file.path(outdir, 'survival.pdf')
        print( plot_survival(
                    object = object, 
                     assay = splitvar, 
                    engine = engine, 
                     ntile = ntile,
                      file = file, 
                     width = width, 
                    height = height,
                         n = n, 
                      nrow = nrow, 
                      ncol = ncol
        ) )
    }
# Return
    object
}


svar_formula <- function(formula, object)  all(all.vars(formula) %in% svars(object))


#' Is package installed?
#' @param pkg package (string)
#' @return TRUE or FALSE
#' @export
installed <- function(pkg){
    txt <- sprintf("        `BiocManager::install('%s')`. Then rerun.", pkg)
    if (requireNamespace(pkg, quietly = TRUE)){  return(TRUE )
    } else {                       message(txt); return(FALSE)
    }
}


#' Plot survival
#' 
#' @param object     SummarizedExperiment
#' @param assay      value in assayNames(object)
#' @param engine    'coxph', 'survdiff' or 'logrank'
#' @param ntile  number of quantiles
#' @param title      string
#' @param subtitle   string
#' @param file       filepath
#' @param width      number
#' @param height     number
#' @param n          number of features to plot
#' @param ncol       number of columns
#' @param nrow       number of rows
#' @return ggplot
#' @examples
#' # ~ survgroup
#'     object <- survex()
#'     object %<>% fit_survival(~ survgroup)
#'     plot_survival(object, formula = ~ survgroup)
#'     plot_survival(object, formula = ~ survgroup)
#' 
#' # ~ exprs2levels
#'     object <- survex()
#'     object %<>% factorize_assay(k = 2)
#'     object %<>% fit_survival(~ exprs2levels)
#'     plot_survival(object, formula = ~ exprs2levels)
#'
#' # ~ sex
#'     object <- survex()
#'     object %<>% fit_survival(~ sex)
#'     plot_survival(object, formula = ~ sex)
#'
#' #' ~ expr2levels + subgroup
#'      object <- survex()
#'      object %<>% factorize_assay(k = 2)
#'      object %<>% fit_survival(~ sex + exprs2levels)
#'      plot_survival(object, formula = ~ sex + exprs2levels, nrow = 2, ncol = 2) + scale_x_continuous(breaks = 0:8)
#' 
#' # Engines
#'     object <- survival_example()
#'     object %<>% fit_survival(engine = c('coxph', 'survdiff', 'logrank'))
#'     plot_survival(object)
#' # Pdf
#'     # plot_survival(object, file = file.path('testdir', 'survival', 'survival.pdf'))
#' @export
plot_survival <- function(
      object, 
      formula = as.formula(sprintf('~%s', assayNames(object)[1])),
      engine = c('coxph', 'survdiff', 'logrank') %>% intersect(fits(object)) %>% extract(1),
        coef = coefs(object, fit = engine)[1],
       title = if (svar_formula(formula, object)) NULL else formula2str(formula) , # svar_formula becomes facethdr
    subtitle = sprintf('%s', paste0(engine, collapse = '      ')),
        file = NULL,
       width = 7,
      height = 7,
           n = if (svar_formula(formula, object)) length(all.vars(formula))  else min(nrow(object),9),
        ncol = if (svar_formula(formula, object)) length(all.vars(formula))  else 3,
        nrow = if (svar_formula(formula, object)) length(all.vars(formula))  else 3
){
# Assert
    if (!installed('ggtext'))   return(NULL) 
    if (!installed('ggstance')) return(NULL)
    assert_is_valid_sumexp(object)
    assert_is_subset(all.vars(formula), c(svars(object), assayNames(object)))
    assert_scalar_subset(engine, fits(object))
    assert_scalar_subset(coef, coefs(object, fit = engine))
    event <- timetoevent <- NULL      # svar
    curOut <- facet <- label <- nalive <- nout <- totDead <- totObs <- survival <- y <- NULL
# Prepare
    object %<>% extract_coef_features(fit = engine, coefs = coef, n = n)
    assayvar <- all.vars(formula) %>% intersect(assayNames(object))
  samplevars <- all.vars(formula) %>% intersect(svars(object))
    if (length(assayvar)==0){
        plotdt <- sdt(object)[, c('sample_id', samplevars, 'timetoevent', 'event'), with = FALSE]
        plotdt[, feature_id := formula2str(formula)]
    } else {
        assert_is_a_string(assayvar)
        assert_character_matrix(assays(object)[[assayvar]], .xname = sprintf('assays(object)$%s', assayvar))
        plotdt <- sumexp_to_longdt(object, assay = assayvar, svars = c(samplevars, 'timetoevent', 'event'), value.name = assayvar)
    }
    plotdt[, alive := 1-event]
    setorderv(plotdt, c('feature_id', all.vars(formula), 'timetoevent', 'alive'))
    
    plotdt[, survivalgroup := do.call(paste, c(.SD, sep = '.')), .SDcols = all.vars(formula) ]
    plotdt[ , totObs   := .N - cumsum(1-event),     by = c('feature_id', 'survivalgroup')   ]
    plotdt[ , totDead := cumsum(event),             by = c('feature_id', 'survivalgroup')   ]
    plotdt <- plotdt[ , .(totObs  = max(totObs), 
                          totDead = max(totDead), 
                          curOut  = sum(event==0)), by = c('feature_id', 'survivalgroup', 'timetoevent')]
    plotdt[, survival := 100*(totObs-totDead)/totObs]
    setorderv(plotdt, c('feature_id', 'survivalgroup', 'timetoevent'))
    plotdt0 <- plotdt[ , .SD[ 1] , by = c('feature_id', 'survivalgroup')][, timetoevent := 0 ][, totDead := 0 ][, survival := 100 ][, curOut := 0]
    plotdtn <- plotdt[ , .SD[.N] , by = c('feature_id', 'survivalgroup')][, timetoevent := max(timetoevent)+1][, curOut := 0]
    plotdtn <- plotdtn[totDead!=totObs]  # vertically end survival curve when all dead
    plotdt <- rbind(plotdt0, plotdt, plotdtn)
# Statistics
    pcols <- pvar(object, fit = engine, coef = coef) # `copy` is very important !
    tcol  <- tvar(object, fit = engine, coef = coef) # without modifies are performed in the object!
    statdt <- if (length(assayvar)==0){  copy(metadata(object)$survival)
              } else {                   fdt(object)[, c('feature_id', pcols, tcol), with = FALSE] }
    statdt[, (pcols) := lapply(.SD, formatC, format = 'g', digits = 2), .SDcols = pcols]
    statdt[, facet := paste0(.SD, collapse = '      '), .SDcols = pcols, by = 'feature_id']
    statdt[, (pcols) := NULL]
    #statdt[, facet := sprintf('%s\n%s', paste0(engine, collapse = spaces(8)), facet)]
    statdt[, facet := sprintf('%s\n%s', feature_id, facet)]
    plotdt %<>% merge(statdt, by = 'feature_id')
    plotdt %<>% extract(order(get(tcol)))
    plotdt[, facet := factor(facet, unique(facet))]
# Plot
    maxtime <- max(plotdt$timetoevent)     # stringi::stri_escape_unicode("°")   # \u00b0
    maxsurvival <- max(plotdt$survival)    # stringi::stri_escape_unicode("†")   # \u2020
    maxtotal <- max(plotdt$totObs)         # stringi::stri_escape_unicode("•")   # \u2022
    maxdigits <- ceiling(log10(maxtotal))
    ndt <- plotdt[, .(totObs  = totObs[1], 
                      totDead = totDead[.N], 
                       nout   = totObs[1] - totObs[.N]), by = c('facet', 'survivalgroup')]
    ndt[ , nalive := totObs-totDead-nout ]
    ndt[, label := sprintf('%d<sup>\u00b0</sup> %d<sup>\u2020</sup> %d<sup>\u2022</sup>', nalive, totDead, nout)]
    survivalgroups <- unique(ndt$survivalgroup)
    colordt <- data.table(survivalgroup = survivalgroups, color = make_colors(survivalgroups))
    ndt %<>% merge(colordt, by = 'survivalgroup')
    ndt[ , label := sprintf("<span style='color:%s'>%s</span>", color, label) ]
    ndt <- ndt[, .(label = paste0(label, collapse = '<br>')), by = 'facet' ]
    nfacets <- nrow(ndt)
    npages <- if (is.null(nrow) | is.null(ncol)) 1 else ceiling(nfacets / nrow / ncol)
    if (!is.null(file))  pdf(file, width = width, height = height)
    for (i in seq_len(npages)){
        p <- ggplot(plotdt) + 
             theme_bw() + 
             facet_wrap_paginate(vars(facet), nrow = nrow, ncol = ncol, page = i) + 
             ggtitle(title, subtitle = subtitle) + 
             theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5),
                  panel.grid  = element_blank()) + 
             #ggtext::geom_richtext(data = ndt, aes(x = maxtime, y = maxsurvival, label = label), 
             #                      hjust = 1, vjust = 1, show.legend = FALSE, label.color = 'NA') +
                # Place text before lines to give the latter more prominence
             geom_step(aes(x = timetoevent, y = survival, group = survivalgroup, color = survivalgroup)) + 
             scale_color_manual(values = colordt$color %>% set_names(colordt$survivalgroup)) #+ 
             #geom_point(data = plotdt[curOut>0], aes(x = timetoevent, y = survival, color = survivalgroup), size = 1, show.legend = FALSE) + 
                # Note that here the dropout is placed after the stepdown.
                # This is because each dropout changes the denominator.
                # So changes the survival percentage.
                # But this approach seems to deviate from convention.
                # survminer flags the dropout before the stepdown.
                # It is possible that a future implementation will switch to that behaviour.
        if (!is.null(file))  print(p)
    }
    if (is.null(file))  return(p) else dev.off()
}

