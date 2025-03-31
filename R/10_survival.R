

#' Survival analysis example
#' @return SummarizedExperiment
#' @examples
#' survobj()
#' @export 
survobj <- function(){
    
    set.seed(1)
    mat <- rbind( GENA = c( rnorm(10,3),
                            rnorm(10,4),  
                            rnorm(10,5), 
                            rnorm(10,6)),   # age + sex increase expression
                  GENB = c( rnorm(10,6),    # age + sex decrease expression
                            rnorm(10,4),  
                            rnorm(10,3), 
                            rnorm(10,2)),
                  GENC = c( rnorm(10,3),
                            rnorm(10,3),  
                            rnorm(10,6), 
                            rnorm(10,6)),   # age increases expression
                  GEND = c( rnorm(10,6),    # age decreases expression
                            rnorm(10,6),
                            rnorm(10,3),
                            rnorm(10,3)),
                  GENE = c( rnorm(10,3),
                            rnorm(10,6),
                            rnorm(10,3),
                            rnorm(10,6)),   # female sex increases expression
                  GENF = c( rnorm(10,6),    # female sex decreases expression
                            rnorm(10,3),
                            rnorm(10,6),
                            rnorm(10,3)), 
                  GENG = c( rnorm(10,3),
                            rnorm(10,3),
                            rnorm(10,6),
                            rnorm(10,3)),   # m:age increases expression
                  GENH = c( rnorm(10,3),    # f:age increases expression
                            rnorm(10,3),
                            rnorm(10,3),
                            rnorm(10,6)),
                  GENI = c( rnorm(10,3),
                            rnorm(10,3),
                            rnorm(10,3),
                            rnorm(10,6)),   # junior:f increases expression
                  GENJ = c( rnorm(10,3),    # senior:f increases expression
                            rnorm(10,6),
                            rnorm(10,3),
                            rnorm(10,3)), 
                  GENK = c( rnorm(10,3),    # flat around three
                            rnorm(10,3),
                            rnorm(10,3),
                            rnorm(10,3)), 
                  GENL = c( rnorm(10,4),    # flat around four
                            rnorm(10,4),
                            rnorm(10,4),
                            rnorm(10,4)), 
                  GENM = c( rnorm(10,5),    # flat around five
                            rnorm(10,5),
                            rnorm(10,5),
                            rnorm(10,5)), 
                  GENN = c( rnorm(10,6),    # flat around six
                            rnorm(10,6),
                            rnorm(10,6),
                            rnorm(10,6))
    )
    object <- SummarizedExperiment(list(exprs = mat))
    fdt(object)$feature_id <- fnames(object)
    object$sample_id <- snames(object)  <- c(sprintf('senior.m.%d', 0:9),  
                                             sprintf('senior.f.%d', 0:9), 
                                             sprintf('junior.m.%d', 0:9), 
                                             sprintf('junior.f.%d', 0:9))
    object$age       <- object$sample_id %>% split_extract_fixed('.', 1)
    object$sex       <- object$sample_id %>% split_extract_fixed('.', 2)
    object$replicate <- object$sample_id %>% split_extract_fixed('.', 3)

    time.senior.m <- c( rep(1,4), rep(2,4), rep(3,2) )
    time.senior.f <- c( rep(2,2), rep(3,4), rep(4,4) )
    time.junior.m <- c( rep(4,4), rep(5,4), rep(7,2) )
    time.junior.f <- c( rep(5,2), rep(6,2), rep(7,6) )
    
    event.senior.m <- rep(1,10)
    event.senior.f <- rep(1,10)    
    event.junior.m <- c(rep(1,8),rep(0,2))
    event.junior.f <- c(rep(1,4),rep(0,6))    
    
    object$timetoevent <- c( time.senior.m,  time.senior.f,  time.junior.m,  time.junior.f )
    object$event       <- c(event.senior.m, event.senior.f, event.junior.m, event.junior.f )
    object %<>% factorize_assay(k = 2)
    object
}



#' Fit onefeature survival 
#' @param sd       data.table
#' @param formula  model formula
#' @examples
#' # Prepare
#'      sd <- sumexp_to_longdt(survobj()[1,], svars = c('timetoevent', 'event'), assay = 'exprs2levels')
#'      sd[ , value := factor(value)]
#' # Survival
#'        .coxph(sd, survival::Surv(timetoevent, event) ~ value)
#'     .survdiff(sd, survival::Surv(timetoevent, event) ~ value)
#'      .logrank(sd, survival::Surv(timetoevent, event) ~ value)
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
    colnames(fitmat) %<>% paste0('~coxph')
    data.table(fitmat)
    #data.table(cbind(fitmat, t(tF), t(pF)))
}


#' @rdname dot-coxph
#' @export
.survdiff <- function(sd, formula){
    
    timetoevent <- NULL
    xvar <- labels(stats::terms(formula))
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
    
    timetoevent <- NULL
    xvar <- labels(stats::terms(formula))
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
#' object <- survobj()
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
    mat %<>% paste0('bin', .)
    dim(mat) <- dim(object)
    dimnames(mat) <- dimnames(object)
# Add
    newassayname <- sprintf('%s%dlevels', assay, k)
    if (verbose)   cmessage('%sAdd `%s`', spaces(8), newassayname)
    assays(object)[[newassayname]] <- mat
    object
}


# Fit survival
# 
# Compute survival effect of svars, exprs, and their interactions
# 
# @param object    SummarizedExperiment
# @param formula   Formula
# @param bins      Number of value bins. Zero means unbinned.
# @param bintype  'factor' or 'numeric'
# @param engine   'coxph', 'survdiff', or 'logrank'
# @param drop      Whether to drop factor varname in coefnames
# @param codingfun (factor) coding function
# @param verbose   TRUE or FALSE
# @examples
# # Load/Transform
#       object <- survobj()
# # coxph{survival}
#      .fit_survival(object)
#      .fit_survival(object, ~ exprs)                         #      expr effect
#      .fit_survival(object, ~ exprs2bins)                    #   exprbin effect
#      .fit_survival(object, ~ exprs2levels)                  # exprlevel effect
#      .fit_survival(object, formula = ~ sex)                 #       sex effect
#      .fit_survival(object, formula = ~ sex + exprs2levels)  #       sex effect ACROSS exprlevels,  exprlevel effect ACROSS sexes.
#      .fit_survival(object, formula = ~ sex / exprs2levels)  # exprlevel effect WITHIN sex,               sex effect ACROSS exprlevels.
#      .fit_survival(object, formula = ~ exprs2levels / sex)  #       sex effect WITHIN exprlevel,   exprlevel effect ACROSS sexes
#      .fit_survival(object, formula = ~ exprs2levels * sex)  #       sex effect differences BETWEEN exprlevels
#      .fit_survival(object, formula = ~ exprs2levels, engine = 'survdiff')
#      .fit_survival(object, formula = ~ exprs2levels, engine = 'logrank')
# # survdiff
#       fit_survival(object, ~ exprs2levels)                  # exprlevel effect
#' @rdname fit_survival
#' @export
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
    event <- timetoevent <- NULL
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

    if (drop){ # drop varname from non-numeric vars
        anum <- assays(object)
        snum <- sdt(object)[, samplevars, with = FALSE]
        anum  %<>% vapply(is.not.numeric, logical(1))
        snum %<>% vapply(is.not.numeric, logical(1))
        anum <- names(anum)[anum]
        snum <- names(snum)[snum]
        anum %<>% intersect(assayvar)
        snum %<>% intersect(samplevars)
        for (var in c(anum, snum)){  
            pat <- sprintf('%s(.+)', var)
            names(fitres) %<>% stri_replace_first_regex(pat, '$1')  
    }}
# Merge    
    if (verbose)  message_df('                      %s', summarize_fit(fitres))
    #if ('expr' %in% all.vars(formula)){  object %<>% merge_fit(fitres)
    #} else {                              metadata(object)$survival <- fitres[, -1] }
    #object
    fitres
}


is.not.numeric <- function(x)  !is.numeric(x)
    

#' Fit/Plot survival
#' 
#' @param object        SummarizedExperiment
#' @param formula       model formula: contains svars/assayNames
#' @param engine       'coxph', 'survdiff' or 'logrank'
#' @param drop          TRUE or FALSE : whether to drop var in coefname
#' @param codingfun     coding function
#' @param verbose       TRUE or FALSE
#' @param outdir        output directory
#' @param plot          TRUE or FALSE
#' @param width         number
#' @param height        number
#' @param n             number of features to plot
#' @param n_col         number of columns
#' @param n_row         number of rows
#' @param writefunname  'write_xl' or 'write_ods'
#' @param order         coefs to order plots
#' @param stats         coefs to print stats for
#' @param title         string
#' @param file          filepath
#' @return SummarizedExperiment/ggplot
#' @examples
#' # survival ~ svars
#'   object <- survobj()
#'   object %>% fit_survival(~age)                          %>% plot_survival(~age)
#'   object %>% fit_survival(~age, engine = 'survdiff')     %>% plot_survival(~age, engine = 'survdiff')
#'   object %>% fit_survival(~sex)     %>% plot_survival(~sex)
#'   object %>% fit_survival(~age+sex) %>% plot_survival(~age+sex)
#'   object %>% fit_survival(~age/sex) %>% plot_survival(~age/sex)
#' 
#' # survival ~ assay
#'   object <- survobj()
#'   object %>% fit_survival(~exprs)
#'   object %>% fit_survival(~exprs2levels) %>% plot_survival(~exprs2levels)
#'
#' # survival ~ svar / assay
#'   object <- survobj()
#'   object %<>% fit_survival(~age/exprs2levels)
#'   object %>% plot_survival(~age/exprs2levels, stats = c('senior:bin2-bin1', 'junior:bin2-bin1'))
#'
#' # Plot
#'   fit_survival(object, ~exprs2levels, plot = TRUE)
#'  #fit_survival(object, ~exprs2levels, plot = TRUE, outdir = 'outdir', n = Inf)
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
             n = if (svar_formula(formula, object)) 1  else min(nrow(object),4), # Inf works
         n_col = n %>% min(nrow(object)) %>% sqrt() %>% floor()   %>% min(4),
         n_row = n %>% min(ncol(object)) %>% sqrt() %>% ceiling() %>% min(4),
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
        if (all(all.vars(formula) %in% svars(object))){  metadata(object)$survival <-  outdt
        } else {                                         object %<>% merge_fdt(outdt)  }
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
        print( plot_survival(   object = object, 
                               formula = formula,
                                engine = engine,
                                  file = file, 
                                 width = width, 
                                height = height,
                                     n = n, 
                                 n_row = n_row, 
                                 n_col = n_col      ) )
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


#' @rdname fit_survival
#' @export
prep_survival <- function(
      object, 
     formula = as.formula(sprintf('~%s', assayNames(object)[1])),
      engine = c('coxph', 'survdiff', 'logrank') %>% intersect(fits(object)) %>% extract(1),
       order = autonomics::coefs(object, fit = engine)[1],
       stats = autonomics::coefs(object, fit = engine),
           n = if (svar_formula(formula, object)) 1  else min(nrow(object),9)#,
#        title = if (svar_formula(formula, object)) NULL else formula2str(formula) , # svar_formula becomes facethdr
#     subtitle = sprintf('%s', paste0(engine, collapse = '      ')),
# dodge_height = 2,
#         file = NULL,
#        width = 7,
#       height = 7,
#         ncol = if (svar_formula(formula, object)) length(all.vars(formula))  else 3,
#         nrow = if (svar_formula(formula, object)) length(all.vars(formula))  else 3
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(all.vars(formula), c(svars(object), assayNames(object)))
    assert_scalar_subset(engine, fits(object))
    assert_is_subset(order, autonomics::coefs(object, fit = engine))
    assert_is_subset(stats, autonomics::coefs(object, fit = engine))
    event <- timetoevent <- NULL      # svar
    curOut <- facet <- label <- nalive <- nout <- totDead <- totObs <- survival <- y <- NULL
    alive <- coef <- p <- NULL
# Prepare
    object %<>% extract_coef_features(fit = engine, coefs = order, n = n)
    assayvar <- all.vars(formula) %>% intersect(assayNames(object))
  samplevars <- all.vars(formula) %>% intersect(svars(object))
    if (length(assayvar)==0){
        plotdt <- sdt(object)[, c('sample_id', samplevars, 'timetoevent', 'event'), with = FALSE]
        plotdt[, feature_id := formula2str(formula)] # allows for generic data.table code
        plotdt[, feature_id := factor(feature_id)]   # allows for levels(.) to work later
    } else {
        assert_is_a_string(assayvar)
        assert_character_matrix(assays(object)[[assayvar]], .xname = sprintf('assays(object)$%s', assayvar))
        plotdt <- sumexp_to_longdt(object, assay = assayvar, svars = c(samplevars, 'timetoevent', 'event'), value.name = assayvar)
    }
    plotdt[, alive := 1-event]
    setorderv(plotdt, c('feature_id', all.vars(formula), 'timetoevent', 'alive'))
    
    plotdt[ , totObs   := .N - cumsum(1-event),     by = c('feature_id', all.vars(formula))   ]
    plotdt[ , totDead := cumsum(event),             by = c('feature_id', all.vars(formula))   ]
    plotdt <- plotdt[ , .(totObs  = max(totObs), 
                          totDead = max(totDead), 
                          curOut  = sum(event==0)), by = c('feature_id', all.vars(formula), 'timetoevent')]
    plotdt[, survival := 100*(totObs-totDead)/totObs]
    setorderv(plotdt, c('feature_id', all.vars(formula), 'timetoevent'))
    plotdt0 <- plotdt[ , .SD[ 1] , by = c('feature_id', all.vars(formula))][, timetoevent := 0 ][, totDead := 0 ][, survival := 100 ][, curOut := 0]
    plotdtn <- plotdt[ , .SD[.N] , by = c('feature_id', all.vars(formula))][, timetoevent := max(timetoevent)+1][, curOut := 0]
    plotdtn <- plotdtn[totDead!=totObs]                  # Vertically end survival curve when all dead
    features <- plotdt[, levels(feature_id)]
    plotdt <- rbind(plotdt0, plotdt, plotdtn)            # Preserve tvar order
    plotdt[, feature_id := factor(feature_id, features)] # Note that for svar-formula feature_id is the formula!
    plotdt <- plotdt[order(feature_id)]
# Statistics
    plongdt <- pdt(object, fit = engine, coef = stats)
    tlongdt <- tdt(object, fit = engine, coef = stats)
    plongdt %<>% melt.data.table(id.vars = 'feature_id', variable.name = 'coef', value.name = 'p')
    tlongdt %<>% melt.data.table(id.vars = 'feature_id', variable.name = 'coef', value.name = 't')
    statdt <- merge(plongdt, tlongdt, by = c('feature_id', 'coef'))
    statdt[, coef := split_extract_fixed(coef, '~', 1)]
    statdt[, p := formatC(p, format = 'g', digits = 2)]
    statdt[sign(t)=='-1', p := sprintf('-%s', p)]
    statdt[, p := sprintf('p = %s', p ) ]
    statdt[sign(t)=='-1', p := sprintf('-%s', p)]
    statdt[, p := stri_pad_both(p, nchar(coef))]
    statdt[, coef := stri_pad_both(coef, nchar(p))]
    statdt <- statdt[, .(coef = paste0(coef, collapse = '        '), 
                            p = paste0(p,    collapse = '        ') ), by = 'feature_id']
    statdt[, facet := sprintf('%s\n%s', coef, p) , by = 'feature_id']
    if (any(all.vars(formula) %in% assayNames(object)))  statdt[, facet := sprintf('%s\n%s', feature_id, facet)]
    plotdt %<>% merge(statdt, by = 'feature_id', sort = FALSE)
    #setorderv(plotdt, tcol)
    plotdt[, facet := factor(facet, unique(facet))]
    plotdt[]
}

#' @rdname fit_survival
#' @export
plot_survival <- function(
      object,
     formula = as.formula(sprintf('~%s', assayNames(object)[1])), 
      engine = c('coxph', 'survdiff', 'logrank') %>% intersect(fits(object)) %>% extract(1),
       order = autonomics::coefs(object, fit = engine)[1],
       stats = autonomics::coefs(object, fit = engine),
       title = sprintf('%s ~ %s', engine, formula2str(formula) %>% substr(2,nchar(.))),
#dodge_height = 0,      # `color` and `linetype` are hardmapped from `all.vars(formula)`
        file = NULL,    #  softmapping them formula-agnostically doesnt work
       width = 7,       #  Only for formula group is sample property (e.g. sex) sharing guaranteed
      height = 7,
           n = if (svar_formula(formula, object)) 1  else min(nrow(object),4), # Inf works
       n_col = n %>% min(nrow(object)) %>% sqrt() %>% floor()   %>% min(4),
       n_row = n %>% min(ncol(object)) %>% sqrt() %>% ceiling() %>% min(4)
){
# Assert
    if (!installed('ggtext'))   return(NULL) 
    if (!installed('ggstance')) return(NULL)
    if (is.infinite(n))  n <- nrow(object)
    totObs <- totDead <- nalive <- nout <- label <- color <- facet <- NULL
    timetoevent <- survival <- NULL
# Plot
    plotdt <- prep_survival(object = object, formula = formula, engine = engine, order = order, stats = stats, n = n)
    maxtime <- max(plotdt$timetoevent)     # stringi::stri_escape_unicode("°")   # \u00b0
    maxsurvival <- max(plotdt$survival)    # stringi::stri_escape_unicode("†")   # \u2020
    maxtotal <- max(plotdt$totObs)         # stringi::stri_escape_unicode("•")   # \u2022
    maxdigits <- ceiling(log10(maxtotal))
    ndt <- plotdt[, .(totObs  = totObs[1], 
                      totDead = totDead[.N], 
                       nout   = totObs[1] - totObs[.N]), by = c('facet', all.vars(formula))]
    ndt[ , nalive := totObs-totDead-nout ]
    ndt[, label := sprintf('%d<sup>\u00b0</sup> %d<sup>\u2020</sup> %d<sup>\u2022</sup>', nalive, totDead, nout)]
    paste. <- function(...) paste(..., sep = '.')
    ndt[ , color := make_colors(do.call(paste., .SD)) , .SDcols = all.vars(formula) ]
    ndt[ , label := sprintf("<span style='color:%s'>%s</span>", color, label) ]
    ndt <- ndt[, .(label = paste0(label, collapse = '<br>')), by = 'facet' ]
    npages <- if (is.null(n_row) | is.null(n_col)) 1 else ceiling(n / n_row/ n_col)
    if (!is.null(file))  pdf(file, width = width, height = height)
    for (i in seq_len(npages)){
        subtitle <- if (svar_formula(formula, object)) NULL else paste0(order, collapse = '  ')
        p <- ggplot(plotdt) + 
             theme_bw() + 
             facet_wrap_paginate(vars(facet), nrow = n_row, ncol = n_col, page = i) + 
             ggtitle(title, subtitle = subtitle) + 
             theme( plot.title    = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5),
                      panel.grid  = element_blank())
             #ggtext::geom_richtext(data = ndt, aes(x = maxtime, y = maxsurvival, label = label), 
             #                      hjust = 1, vjust = 1, show.legend = FALSE, label.color = 'NA') +
                # Place text before lines to give the latter more prominence
        groupsyms <- syms(all.vars(formula))
         colorsym <- sym( all.vars(formula)[[1]])
      linetypesym <- if (length(all.vars(formula))<2) quo(NULL) else sym( all.vars(formula)[[2]])
        p <- p + geom_step( mapping = aes(  x = timetoevent, 
                                            y = survival,              
                                        group = interaction(!!!groupsyms),  # !!! for syms
                                        color = !!colorsym,                 #  !! for sym
                                     linetype = !!linetypesym )) #, 
                           #position = ggstance::position_dodgev(dodge_height))   # dodging slows code
                            #+ 
             #scale_color_manual(values = colordt$color %>% set_names(colordt$color)) #+ 
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

