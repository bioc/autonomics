
empty_survplot <- function(){
    # https://stackoverflow.com/questions/61907987/produce-empty-plot-with-ggsurvplot
    dt  <- data.table(exprlevel = c(rep("low", 10), rep("high", 10)), 
                      value   = c(rnorm(10,mean = 2), rnorm(10,mean = 3)))
    fit <- survival::survfit(survival::Surv(value) ~ exprlevel, data = dt)
    survminer::ggsurvplot(
        fit, 
        data             = dt, 
        surv.median.line = "none", 
        palette          = rep("white", 2), 
        legend           = "none")
}


.dichotomize_exprs <- function(subdt, percentile){
    value <- NULL
    if (all(is.na(subdt$value)))  return(cbind(subdt, exprlevel = 'no values available'))
    subdt %<>% extract(!is.na(value))
    subdt %<>% extract(order(value))
    n <- floor(0.01*percentile*nrow(subdt))
    lowervalue <- subdt$value[n]
    uppervalue <- rev(subdt$value)[n]
    if (length(lowervalue)==0 | lowervalue==uppervalue){
        subdt <- cbind(subdt[0], exprlevel = character(0))
    } else {
        lowergroup <- paste0('Lo', percentile) #paste0('<', signif(lowervalue,1))
        uppergroup <- paste0('Hi', percentile) #paste0('>', signif(uppervalue,1)) 
        subdt <- rbind(cbind(subdt[value<=lowervalue], exprlevel = lowergroup),
                       cbind(subdt[value>=uppervalue], exprlevel = uppergroup))
        #subdt$exprlevel %<>% factor(c(lowergroup, uppergroup))
    }
    subdt
}

dichotomize_exprs <- function(dt, percentile){
    dt %<>% extract(, .dichotomize_exprs(.SD, percentile = percentile), by = 'feature_id')
    dt[, exprlevel := factor(exprlevel, sprintf('%s%d', c('Lo', 'Hi'), percentile))]
    dt[]
}

.fit_survival <- function(subdt, sep, samples = FALSE){
    timetoevent <- event <- exprlevel <- NULL
    logrank <- suppressWarnings(survdiff(Surv(timetoevent, event) ~ exprlevel, data = subdt))
     cph <- suppressWarnings(coef(summary(coxph(Surv(subdt$timetoevent, subdt$event)~subdt$value))))
    exprlevels <- levels(subdt$exprlevel)
    # We want right tail logrank pvalue only
    # Left tail indicates more similarity than expected
    # Can be interesting to detect fraud (dropping bad data)
    # But not useful for our purpose
    # https://stats.stackexchange.com/questions/22347/is-chi-squared-always-a-one-sided-test
    dt <- data.table( `p~hi-lo~logrank` = 1 - pchisq(logrank$chisq, 1), # we want onesided value
                 `effect~hi-lo~logrank` = logrank$chisq * sign(cph[,'coef']),
                      `t~hi-lo~logrank` = logrank$chisq * sign(cph[,'coef']),
                      `p~expr~cph` = cph[,'Pr(>|z|)'],
                 `effect~expr~cph` = cph[,'coef'    ], 
                      `t~expr~cph` = cph[,'z'       ] )   # effect, p
    if (samples){
         low <- unique(subdt[exprlevel == exprlevels[1]])$sample_id
        high <- unique(subdt[exprlevel == exprlevels[2]])$sample_id
         low %<>% as.character() %>% commonify_strings()
        high %<>% as.character() %>% commonify_strings()
        dt[ , lo := low  ] #  lower survival samples
        dt[ , hi := high ] # higher survival samples
    }
    return(dt)
}

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
                        data.table( subgroup = 'Diseased', sample_id = 'D05', timetoevent = 2, event = 1),
                        data.table( subgroup = 'Diseased', sample_id = 'D06', timetoevent = 2, event = 1),
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

#' Fit/Plot survival 
#' @param object      SummarizedExperiment
#' @param assay       string
#' @param percentile  percentage (not greater than 50)
#' @param sep         fvar string separator : e.g. '~' gives p~surv~LR50 
#' @param samples     TRUE or FALSE : record which samples in which stratum ?
#' @param verbose     TRUE or FALSE
#' @param outdir      dir
#' @param writefunname 'write_xl' or 'write_ods'
#' @param plot        TRUE or FALSE
#' @param n           number
#' @param ncol        number
#' @param nrow        number
#' @param file        filepath
#' @param width       number
#' @param height      number
#' @param title       string
#' @param subtitle    string
#' @param palette     color vector
#' @param conf.int    TRUE or FALSE
#' @return ggsurvplot
#' @examples 
#' object <- survival_example()
#' fit_survival(object)
#' @export
fit_survival <- function(
        object, 
         assay = assayNames(object)[1],
    percentile = 25, 
           sep = FITSEP,
       samples = if (ncol(object) < 50) TRUE else FALSE,
       verbose = TRUE, 
        outdir = NULL,
  writefunname = 'write_xl',
          plot = if (is.null(outdir)) FALSE else TRUE,
             n = 4,
          ncol = 4,
          nrow = length(percentile),
         width = 7*ncol,
        height = 7*nrow
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, assayNames(object))
    assert_is_numeric(percentile)
    assert_all_are_in_left_open_range(percentile, 0, 50)
    event <- exprlevel <- timetoevent <- value <- NULL
# Fit
    if (verbose)  cmessage('%ssurvival ~ exprs  cphmodel', spaces(8))                         # Filter across
    object %<>% filter_samples(!is.na(event) & !is.na(timetoevent))
    for (pct in percentile){
        dt <- sumexp_to_longdt(object, assay = assay, svars = c('event', 'timetoevent'))       # Melt
        if (verbose)  cmessage("%s~ expr%d logranktest", spaces(17), pct)   # Dichotomize
        dt %<>% dichotomize_exprs(percentile = pct)                                            # Filter within 
       #dt <- dt[, .SD[sum(event==1 & !is.na(value))>=3], by = c('feature_id', 'exprlevel')]   #    3 events     per feature/exprlevel
        dt <- dt[, .SD[    length(unique(exprlevel))==2], by = c('feature_id')             ]   #    2 exprlevels per feature
        dt %<>% extract(, .fit_survival(.SD, sep = sep, samples = samples), by = 'feature_id') # Fit survival
        oldnames <- newnames <- names(dt)
        newnames %<>% stri_replace_all_fixed( 'hi-', sprintf('hi%d-', pct))
        newnames %<>% stri_replace_all_fixed( '-lo', sprintf('-lo%d', pct))
        newnames %<>% stri_replace_all_regex('^hi$', sprintf('hi%d', pct))
        newnames %<>% stri_replace_all_regex('^lo$', sprintf('lo%d', pct))
        setnames(dt, oldnames, newnames) 
        for (col in newnames)  object[[col]] <- NULL
        if (verbose)  message_df('                                   %s', 
                                 summarize_fit(dt, fit = c('logrank', 'cph')))
        object %<>% merge_fdt(dt)
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
        plot_survminer(object = object, assay = assay, n = n, 
                      ncol = ncol, nrow = nrow, width = width, height = height, file = file)
    }
# Return
    object
}


#' Plot survival
#' @param object   SummarizedExperiment
#' @param assay    value in assayNames(object)
#' @param coefs    autonomics::coefs(object) subset
#' @param splitvar split svar
#' @param title    string
#' @examples
#' object <- survival_example()
#' object %<>% fit_survival(percentile = 25)
#' plot_survival(object)
#' @export
plot_survival <- function(
      object, 
       assay = assayNames(object)[1],
  percentile = 25,
    splitvar = 'exprlevel',
       title = sprintf('survival ~ expr'),
    ordervar = tvar(object, fit = 'cph', coef = 'expr'),
        pvar = c('p~expr~cph', sprintf('p~hi%d-lo%d~logrank', percentile, percentile))
                 
){
# Prepare
    logrankvar <- sprintf('p~hi%d-lo%d~logrank', percentile, percentile)
    cphvar <- 'p~expr~cph'
    plotdt <- sumexp_to_longdt(object, assay = assay, svars = c('timetoevent', 'event'))
    if (!is.null(percentile))  plotdt %<>% dichotomize_exprs(percentile = percentile)
    plotdt %<>% extract(order(feature_id, get(splitvar), timetoevent))
    plotdt[ , ntotal := .N , by = c('feature_id', splitvar)]
    plotdt <- plotdt[ , .(ntotal = unique(ntotal),                    ndead = sum(event)) ,   by = c('feature_id', splitvar, 'timetoevent')]
    plotdt <- plotdt[ , .(ntotal = ntotal, timetoevent = timetoevent, ndead = cumsum(ndead)), by = c('feature_id', splitvar)]
    plotdt[, survival := 100*(ntotal-ndead)/ntotal]
    plotdt %<>% extract(order(feature_id, get(splitvar), timetoevent))
    plotdt0 <- plotdt[ , .SD[ 1] , by = c('feature_id', splitvar)][, timetoevent := 0 ][, ndead := 0 ][, survival := 100 ]
    plotdtn <- plotdt[ , .SD[.N] , by = c('feature_id', splitvar)][, timetoevent := max(timetoevent)+1]
    plotdt <- rbind(plotdt0, plotdt, plotdtn)
# Statistics
    plotdt %<>% merge(fdt(object)[, .SD, .SDcols = patterns('feature_id|logrank|cph')], by = 'feature_id', sort = FALSE)
    plotdt[, facet := feature_id]
    plotdt[, facet := sprintf('%s\ncph                lr%d', facet, percentile)]
    plotdt[, facet := sprintf('%s\n%s         %s', facet, formatC(get(cphvar),     format = 'e', digits = 1),
                                                          formatC(get(logrankvar), format = 'e', digits = 1)) , by = 'feature_id']
    plotdt %<>% extract(order(get(ordervar)))
    plotdt[, facet := factor(facet, unique(facet))]
# Plot
    # pdt <- plotdt[, .(label = sprintf('cph %s\nlr%d %s', formatC(get(cphvar    )[1], format = 'e', digits = 1), 
    #                                          percentile, formatC(get(logrankvar)[1], format = 'e', digits = 1)),
    #                       x = median(timetoevent),
    #                       y = median(survival)), by = 'facet']
     ndt <- plotdt[, .(x = min(timetoevent), 
                       y = max(survival)*(1+0.1-0.2*as.numeric(get(splitvar))), 
                   label = sprintf('%d', ntotal[1])), by = c('facet', splitvar)]
    ggplot(plotdt) + 
         theme_bw() + 
         facet_wrap(vars(facet)) + 
         ggtitle(title) + 
         theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
              panel.grid  = element_blank()) + 
         geom_step(aes(x = timetoevent, y = survival, group = !!sym(splitvar), color = !!sym(splitvar))) + 
         geom_text(data = ndt, aes(x = x, y = y, label = label, color = !!sym(splitvar)), hjust = -0.1, vjust = 1, show.legend = FALSE)
#         geom_text(data = pdt, aes(x = -Inf, y = -Inf, label = label), hjust = -0.1, vjust = -0.2)
}


.plot_survminer <- function(
        object,
         assay = assayNames(object)[1],
         coefs = autonomics::coefs(object, fit = 'logrank'),
         title = paste0(assay, ' ', coefs),
      subtitle = NULL,
       palette = c("#009999", "#ff5050"),
      conf.int = FALSE
){
# Assert
    if (!requireNamespace('survminer', quietly = TRUE)){
        message("BiocManager::install('survminer'). Then re-run.")
        return(object) 
    }
    assert_is_valid_sumexp(object)
    if (nrow(object)==0)  return(empty_survplot())
    assert_is_subset(c('event', 'timetoevent'), svars(object))
    assert_is_identical_to_true(nrow(object)==1)
    feature <- unique(fdata(object)$feature_id)
    title %<>% paste(feature, ., sep = ' : ')
    assert_is_scalar(feature)
    value <- exprlevel <- NULL
# Prepare
    subdt <- sumexp_to_longdt( object, assay = assay, svars = c('event', 'timetoevent') )
    subdt %<>% dichotomize_exprs( percentile = as.numeric(substr(coefs,3,4)) )
# Plot
    fit <- survfit(Surv(timetoevent, event) ~ exprlevel, data = subdt)
    # legend.labs is passed to ggtext
    # ggtext thinks `<` is a tag to be parsed, but is able to process it and throws an error.
    # Adding \u200B (a zero space unicode character) breaks the tag and fixes the error.
    # Source: https://stackoverflow.com/questions/67890410
    legend.labs <- sprintf('\u200B%s', unique(subdt$exprlevel))
    survminer::ggsurvplot(
        fit, data = subdt, conf.int = conf.int, palette = palette,
        risk.table = TRUE, risk.table.col = 'strata', risk.table.height = 0.25, 
        pval = TRUE, ggtheme = theme_bw(), title = title, subtitle = subtitle,
        legend.labs = legend.labs, legend.title = assay)
}


#' survival percentiles
#' @param object SummarizedExperiment
#' @return numeric vector
#' @export
percentiles <- function(object)  as.numeric(substr(coefs(object, fit = 'logrank'), 3,4))


plot_survminer <- function(
        object, 
         assay = assayNames(object)[1],
         coefs = autonomics::coefs(object, fit = 'logrank'),
         title = paste0(assay, ' ', coefs),
      subtitle = NULL,
       palette = c("#009999", "#ff5050"),
      conf.int = FALSE,
             n = 4,
          ncol = 4, 
          nrow = length(coefs), 
          file = NULL, 
         width = 7*ncol, 
        height = 7*nrow
    
){
# Extract
    object %<>% order_on_p(fit = 'logrank', coefs = coefs[1], verbose = FALSE)
    n %<>% min(nrow(object))
    object %<>% extract(1:n, )
    object %<>% order_on_t(fit = 'logrank', coefs = coefs[1], verbose = FALSE)
# Plot
    if (!is.null(file)){
        cmessage('%s%s', spaces(21), file)
        pdf(file, width = width, height = height)
    }
    npages <- ceiling(nrow(object)/ncol)
    for (i in 1:npages){
        cmessage('\t\t\tPage %02d/%02d', i, npages)
        idx1 <- (i-1)*ncol+1
        idxn <- min(i*ncol, nrow(object))
        idx <- idx1:idxn
        objlist <- object[idx, ]
        objlist %<>% split_features(by = 'feature_id')
        plots <- mapply(
            .plot_survminer, 
            object     = rep(objlist, each = length(coefs)), 
            coefs = rep(coefs, times = length(objlist)),
            MoreArgs = list(assay = assay, palette = palette, conf.int = conf.int), SIMPLIFY = FALSE)
        survminer::arrange_ggsurvplots(plots, nrow = nrow, ncol = ncol)
    }
    if (!is.null(file))  dev.off()
}


