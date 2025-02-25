
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
#'      object <- survival_example()
#'      sd <- sumexp_to_longdt(object, svars = c('timetoevent', 'event'))
#'      sd %<>% extract(feature_id %in% feature_id[1])
#'      sd[ , quantile := factor(dplyr::ntile(value, 2)) ]
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


#' Survival engines
#' @export
#' @examples
#' SURVIVALENGINES
SURVIVALENGINES <- c('coxph', 'survdiff', 'logrank')


#' Bin/Factorize assay
#' @param object  SummarizedExperiment
#' @param assay   string
#' @param k       number of bins/levels
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' object <- survival_example()
#'       bin_assay(object)
#' factorize_assay(object)
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
    newassayname <- paste0(assay, '3bins')
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
    newassayname <- paste0(assay, '3levels')
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
#'    object <- survival_example()
#'    object %<>% bin_assay()
#'    object %<>% factorize_assay()
#' # coxph{survival}
#'   .fit_survival(object)
#'   .fit_survival(object, ~ exprs)                             #      expr effect
#'   .fit_survival(object, ~ exprs3bins)                        #   exprbin effect
#'   .fit_survival(object, ~ exprs3levels)                      # exprlevel effect
#'   .fit_survival(object, formula = ~ subgroup)                #  subgroup effect
#'   .fit_survival(object, formula = ~ subgroup + expr3levels)  #  subgroup effect ACROSS exprlevels,  exprlevel effect ACROSS subgroups.
#'   .fit_survival(object, formula = ~ subgroup / expr3levels)  # exprlevel effect WITHIN subgroup,     subgroup effect ACROSS exprlevels.
#'   .fit_survival(object, formula = ~ expr3levels / subgroup)  #  subgroup effect WITHIN exprlevel,   exprlevel effect ACROSS subgroups.
#'   .fit_survival(object, formula = ~ expr3levels * subgroup)  #  subgroup effect differences BETWEEN exprlevels
#' # survdiff{survival}
#'   .fit_survival(object, formula = ~ exprs3levels, engine = 'survdiff')
#'   .fit_survival(object, formula = ~ exprs3levels, engine = 'logrank')
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
    if (engine == 'logrank'){
        if (!requireNamespace('coin', quietly = TRUE))  message("BiocManager::install('coin'). Then rerun")
    }
    object %<>% filter_samples(!is.na(event) & !is.na(timetoevent))
# Code
    survivalvars <- c('timetoevent', 'event')
      samplevars <- intersect(all.vars(formula),      svars(object))
        assayvar <- intersect(all.vars(formula), assayNames(object))
    if (is_empty(assayvar)){   
        dt <- sdt(object)[, c('sample_id', samplevars, survivalvars), with = FALSE]
        dt[, feature_id := 'dummy'] 
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
    if ('expr' %in% all.vars(formula)){  object %<>% merge_fit(fitres)
    } else {                              metadata(object)$survival <- fitres[, -1] }
    object
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
         ntile = 2,
        engine = c('survdiff', 'coxph', 'logrank')[1],
      splitvar = assayNames(object)[1],
           sep = FITSEP,
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
    for (sva in splitvar){
    for (eng in engine){
    for (nti in ntile){
        outdt <- .fit_survival(  object = object, 
                               splitvar = sva, 
                                  ntile = nti, 
                                 engine = eng, 
                                verbose = FALSE  )
                 object %<>% merge_fdt(outdt[feature_id != 'dummy'])
        metadata(object)$survival <-  outdt[feature_id == 'dummy']
    }}}
# Message
    if (verbose)  message_df('                                   %s', 
                             summarize_fit(  rbind(fdt(object), metadata(object)$survival),
                                             fit = engine ))
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
#' # Defaults
#'     object <- survival_example()
#'     object %<>% fit_survival()
#'     plot_survival(object)
#' # Engines
#'     object <- survival_example()
#'     object %<>% fit_survival(engine = c('coxph', 'survdiff', 'logrank'))
#'     plot_survival(object)
#' # Pdf
#'     # plot_survival(object, file = file.path('testdir', 'survival', 'survival.pdf'))
#' @export
plot_survival <- function(
      object, 
       assay = assayNames(object)[1],
      engine = intersect(fits(object), c('coxph', 'survdiff', 'logrank')),
       ntile = 2,
       title = sprintf('surv ~ expr'), 
    subtitle = sprintf('%s', paste0(engine, collapse = '      ')),
        file = NULL,
       width = 7,
      height = 7,
           n = min(nrow(object),9),
        ncol = 3,
        nrow = 3
){
# Prevent check notes
    if (!requireNamespace('ggtext', quietly = TRUE)){   
        message("BiocManager::install('ggtext'). Then rerun")
        return(NULL)
    }
    event <- timetoevent <- NULL      # svar
    value <- NULL                     # sumexp_to_longdt                                     # plotdt
    color <- curOut <- facet <- label <- nalive <- nout <- totDead <- totObs <- survival <- y <- NULL
# Prepare
    obj <- extract_coef_features(object, fit = engine[1], n = n)
    plotdt <- sumexp_to_longdt(obj, assay = assay, svars = c('timetoevent', 'event'))
    plotdt[, quantile := NA_character_]
    plotdt[, quantile := dplyr::ntile(value, ntile), by = 'feature_id']
    plotdt <- plotdt[quantile %in% c(1, ntile)]
    plotdt %<>% extract(order(feature_id, quantile, timetoevent))
    plotdt <- plotdt[order(feature_id, quantile, timetoevent, -event)]
    plotdt[ , totObs   := .N - cumsum(1-event),     by = c('feature_id', 'quantile')   ]
    plotdt[ , totDead := cumsum(event),             by = c('feature_id', 'quantile')   ]
    plotdt <- plotdt[ , .(totObs  = max(totObs), 
                          totDead = max(totDead), 
                          curOut  = sum(event==0)), by = c('feature_id', 'quantile', 'timetoevent')]
    plotdt[, survival := 100*(totObs-totDead)/totObs]
    plotdt %<>% extract(order(feature_id, quantile, timetoevent))
    plotdt0 <- plotdt[ , .SD[ 1] , by = c('feature_id', 'quantile')][, timetoevent := 0 ][, totDead := 0 ][, survival := 100 ][, curOut := 0]
    plotdtn <- plotdt[ , .SD[.N] , by = c('feature_id', 'quantile')][, timetoevent := max(timetoevent)+1][, curOut := 0]
    plotdt <- rbind(plotdt0, plotdt, plotdtn)
# Statistics
    pcols <- pvar(object, fit = engine)
    tcol  <- tvar(object, fit = engine[1])
    statdt <- fdt(object)[, c('feature_id', pcols, tcol), with = FALSE]
    statdt[, (pcols) := lapply(.SD, formatC, format = 'g', digits = 2), .SDcols = pcols]
    statdt[, facet := paste0(.SD, collapse = '      '), .SDcols = pcols, by = 'feature_id']
    statdt[, (pcols) := NULL]
    #statdt[, facet := sprintf('%s\n%s', paste0(engine, collapse = spaces(8)), facet)]
    statdt[, facet := sprintf('%s\n%s', feature_id, facet)]
    plotdt %<>% merge(statdt, by = 'feature_id')
    plotdt %<>% extract(order(get(tcol)))
    plotdt[, facet := factor(facet, unique(facet))]
    plotdt[, quantile := paste0('Q', quantile)]
# Plot
    maxtime <- max(plotdt$timetoevent)     # stringi::stri_escape_unicode("°")   # \u00b0
    maxsurvival <- max(plotdt$survival)    # stringi::stri_escape_unicode("†")   # \u2020
    maxtotal <- max(plotdt$totObs)         # stringi::stri_escape_unicode("•")   # \u2022
    maxdigits <- ceiling(log10(maxtotal))

    ndt <- plotdt[, .(totObs  = totObs[1], 
                      totDead = totDead[.N], 
                      nout   = totObs[1] - totObs[.N]), by = c('facet', 'quantile')]
    ndt[ , nalive := totObs-totDead-nout ]
    ndt[, label := sprintf('%d<sup>\u00b0</sup> %d<sup>\u2020</sup> %d<sup>\u2022</sup>', nalive, totDead, nout)]
    quantiles <- unique(ndt$quantile)
    colordt <- data.table(quantile = quantiles, color = make_colors(quantiles))
    ndt %<>% merge(colordt, by = 'quantile')
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
             ggtext::geom_richtext(data = ndt, aes(x = maxtime, y = maxsurvival, label = label), 
                                   hjust = 1, vjust = 1, show.legend = FALSE, label.color = 'NA') +
                # Place text before lines to give the latter more prominence
             geom_step(aes(x = timetoevent, y = survival, group = quantile, color = quantile)) + 
             geom_point(data = plotdt[curOut>0], aes(x = timetoevent, y = survival, color = quantile), size = 1, show.legend = FALSE)
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

