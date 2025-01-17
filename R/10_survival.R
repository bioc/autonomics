
#' @rdname fit_survival
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
#' @param timetoevent numeric (time to event)
#' @param event       numeric (1=event, 0=not)
#' @param expr        numeric (.coxph) or twolevel-factor (.survdiff, .logrank_test)
#' @examples
#' # Prepare
#'          object <- survival_example()
#'     timetoevent <- object$timetoevent
#'           event <- object$event
#'            expr <- values(object)[1,]
#'        quantile <- factor(dplyr::ntile(expr, 2))
#' # Survival
#'        .coxph(timetoevent, event, expr)
#'     .survdiff(timetoevent, event, quantile)
#'      .logrank(timetoevent, event, quantile)
#' # Sumexp
#'          fit_survival(object)
#' @rdname dot-coxph
#' @export
.coxph <- function(timetoevent, event, expr){
    survout <- suppressWarnings(stats::coef(summary(coxph(Surv(timetoevent, event)~expr))))
    data.table( `p~expr~coxph` =    survout[, 'Pr(>|z|)'], 
           `effect~expr~coxph` = -1*survout[, 'coef'    ], 
                `t~expr~coxph` = -1*survout[, 'z'       ])
}

#' @rdname dot-coxph
#' @export
.survdiff <- function(timetoevent, event, expr){
    nexpr <- length(levels(expr))
    survout <- suppressWarnings(survdiff(   Surv(timetoevent, event) ~ expr   ))
    meandiff <- mean(timetoevent[expr==rev(levels(expr))[1]]) - 
                mean(timetoevent[expr==   (levels(expr))[1]])
    dtout <- data.table(  `p~expr~survdiff` =  1 - pchisq(survout$chisq, 1), 
                     `effect~expr~survdiff` = meandiff,
                          `t~expr~survdiff` = survout$chisq * sign(meandiff) )
    oldnames <- newnames <- names(dtout)
    newnames %<>% stri_replace_first_fixed('~expr~', sprintf('~expr%s~', rev(levels(expr))[1]))
    setnames(dtout, oldnames, newnames)
    dtout[]
}


#' @rdname dot-coxph
#' @export
.logrank <- function(timetoevent, event, expr){
    nexpr <- length(levels(expr))
    survout <- suppressWarnings(coin::logrank_test(   Surv(timetoevent, event) ~ expr   ))
    meandiff <- mean(timetoevent[expr==rev(levels(expr))[1]]) - 
                mean(timetoevent[expr==   (levels(expr))[1]])
    dtout <- data.table(  `p~expr~logrank` = coin::pvalue(survout), 
                     `effect~expr~logrank` = meandiff,
                          `t~expr~logrank` = coin::statistic(survout) * sign(meandiff) )
    oldnames <- newnames <- names(dtout)
    newnames %<>% stri_replace_first_fixed('~expr~', sprintf('~expr%s~', rev(levels(expr))[1]))
    setnames(dtout, oldnames, newnames)
    dtout[]
}


#' Survival engines
#' @export
#' @examples
#' SURVIVALENGINES
SURVIVALENGINES <- c('coxph', 'survdiff', 'logrank')


#' Fit/Plot survival 
#' @param object      SummarizedExperiment
#' @param engine     'coxph' (survival), 'survdiff' (survival), 'logrank' (coin)
#' @param ntile       number
#' @param assay       string
#' @param sep         fvar string separator : e.g. '~' gives p~surv~LR50 
#' @param verbose     TRUE or FALSE
#' @param plot        TRUE or FALSE
#' @param outdir      dir
#' @param writefunname 'write_xl' or 'write_ods'
#' @return ggsurvplot
#' @examples
#' # Defaults
#'     object <- survival_example()
#'     fit_survival(object)
#' # Engines
#'     fit_survival(object, engine = c('coxph', 'survdiff'))
#'     fit_survival(object, engine = c('coxph', 'survdiff', 'logrank'))
#' # Quantiles
#'     fit_survival(object, engine = 'logrank')
#'     fit_survival(object, engine = 'logrank', ntile = 4)
#' # Plot
#'     fit_survival(object)
#'     fit_survival(object, plot = TRUE)
#'     fit_survival(object, engine = c('coxph', 'survdiff', 'logrank'), plot = TRUE)
#' @export
fit_survival <- function(
        object, 
         ntile = 2,
        engine = c('survdiff', 'coxph', 'logrank')[1:2],
         assay = assayNames(object)[1],
           sep = FITSEP,
       verbose = TRUE,
          plot = if (is.null(outdir)) FALSE else TRUE,
        outdir = NULL,
  writefunname = 'write_xl'
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(engine, c('coxph', 'survdiff', 'logrank'))
    assert_scalar_subset(assay, assayNames(object))
    event <- exprlevel <- timetoevent <- value <- NULL
    if ('logrank' %in% engine){
        if (!requireNamespace('coin'))  message("BiocManager::install('coin'). Then rerun")}
# Prepare
    if (verbose)  cmessage('%sSurvival', spaces(8))
    object %<>% filter_samples(!is.na(event) & !is.na(timetoevent))       # Filter
    dt <- sumexp_to_longdt(object, svars = c('timetoevent', 'event'))
# Analyze
    dt[, quantile := dplyr::ntile(value, ntile), by = 'feature_id']   # Quantile
    dt <- dt[quantile %in% c(1, ntile)]
   #dt <- dt[, .SD[sum(event==1 & !is.na(value))>=3], by = c('feature_id', 'quantile')]  #    3 events     per feature/exprlevel
    dt <- dt[, .SD[    length(unique(na.exclude(quantile)))==2], by = c('feature_id')]   #    2 exprlevels per feature
    dt[, quantile := factor(quantile)]
    txt <- '                                   %s'
    if ('survdiff' %in% engine){
        if (verbose)  cmessage('%ssurvdiff%d: surv ~ ntile(exprs,%d)', spaces(8+8), ntile, ntile)
        outdt <- dt[ , .survdiff(timetoevent, event, quantile), by = 'feature_id' ]
        object %<>% merge_fdt(outdt)
    }
    if ('logrank'  %in% engine){
        if (verbose)  cmessage('%slogrank%d: surv ~ ntile(exprs,%d)', spaces(8+8+1), ntile, ntile)
        outdt <- dt[ ,  .logrank(timetoevent, event, quantile), by = 'feature_id' ]
        object %<>% merge_fdt(outdt)
    }
    if ('coxph' %in% engine){
        if (verbose)  cmessage('%scoxph: surv ~ exprs', spaces(8+8+4))
        outdt <- dt[ , .coxph(timetoevent, event, value), by = 'feature_id' ]
        object %<>% merge_fdt(outdt)
    }
    if (verbose)  message_df(txt, summarize_fit(object, fit = engine))
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
        print( plot_survival( object = object, 
                               assay = assay, 
                                file = file, 
                              engine = engine, 
                               ntile = ntile ) )
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
    event <- timetoevent <- NULL                                  # svar
    value <- NULL                                                 # sumexp_to_longdt
    facet <- label <- totDead <- totObs <- survival <- y <- NULL    # plotdt 
# Prepare
    obj <- extract_coef_features(object, fit = engine[1], n = n)
    plotdt <- sumexp_to_longdt(obj, assay = assay, svars = c('timetoevent', 'event'))
    plotdt[, quantile := NA_character_]
    plotdt[, quantile := dplyr::ntile(value, ntile), by = 'feature_id']
    plotdt <- plotdt[quantile %in% c(1, ntile)]
    plotdt %<>% extract(order(feature_id, quantile, timetoevent))
    plotdt <- plotdt[order(feature_id, quantile, timetoevent, -event)]
    plotdt[ , totObs   := .N - cumsum(1-event),                                               by = c('feature_id', 'quantile')   ]
    plotdt[ , totDead := cumsum(event),                                                       by = c('feature_id', 'quantile')   ]
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
    maxtime <- max(plotdt$timetoevent)
    maxsurvival <- max(plotdt$survival)
    maxtotal <- max(plotdt$totObs)
    maxdigits <- ceiling(log10(maxtotal))

    ndt <- plotdt[, .(totObs  = totObs[1], 
                      totDead = totDead[.N], 
                      nout   = totObs[1] - totObs[.N]), by = c('facet', 'quantile')]
    ndt[ , nalive := totObs-totDead-nout ]
    ndt[, label := sprintf('%d<sup>°</sup> %d<sup>†</sup> %d<sup>•</sup>', nalive, totDead, nout)]
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

