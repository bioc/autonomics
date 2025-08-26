

#' Example objects for binding
#' @return SummarizedExperiment
#' @examples
#' object1()
#' object2()
#' @export
object1 <- function(){
    
    obj <- survobj(verbose = FALSE)                   # Differentiate
    sdt(obj)$sampleid  <- seq_len(ncol(obj))          #     svars
    fdt(obj)$featureid <- seq_len(nrow(obj))          #     fvars
    sdt(obj)$sex %<>% toupper()                       #     svalues
    fdt(obj)$gene <- paste0('GENE ', LETTERS[1:14])   #     fvalues
    obj %<>% extract(,-1)                             #     samples
    obj %<>% extract(-1,)                             #     features
    assays(obj)$exprs2levels <- NULL                  #     assaynames
    obj
}


#' @rdname object1
#' @export
object2 <- function(){
    
    obj <- survobj(verbose = FALSE)                   # Differentiate
    sdt(obj)$sid <- seq_len(ncol(obj))                #     svars
    fdt(obj)$fid <- seq_len(nrow(obj))                #     fvars
    sdt(obj)$sex %<>% tolower()                       #     svalues
    fdt(obj)$gene <- paste0('gene ', letters[1:14])   #     fvalues
    obj %<>% extract(,-ncol(.))                       #     samples
    obj %<>% extract(-nrow(.),)                       #     features
    assays(obj)$exprs2bins <- NULL                    #     assaysnames
    obj
    
    
}


#' Survival analysis example
#' @param verbose TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' survobj()
#' @export 
survobj <- function(verbose = TRUE){
    set.seed(1)
    mat <- rbind(   GENA = c( rnorm(10,3), rnorm(10,4), rnorm(10,5), rnorm(10,6)), # age + sex increase expression
                    GENB = c( rnorm(10,6), rnorm(10,4), rnorm(10,3), rnorm(10,2)), # age + sex decrease expression
                    GENC = c( rnorm(10,3), rnorm(10,3), rnorm(10,6), rnorm(10,6)), # age increases expression
                    GEND = c( rnorm(10,6), rnorm(10,6), rnorm(10,3), rnorm(10,3)), # age decreases expression
                    GENE = c( rnorm(10,3), rnorm(10,6), rnorm(10,3), rnorm(10,6)), # female sex increases expression
                    GENF = c( rnorm(10,6), rnorm(10,3), rnorm(10,6), rnorm(10,3)), # female sex decreases expression
                    GENG = c( rnorm(10,3), rnorm(10,3), rnorm(10,6), rnorm(10,3)), # m:age increases expression
                    GENH = c( rnorm(10,3), rnorm(10,3), rnorm(10,3), rnorm(10,6)), # f:age increases expression
                    GENI = c( rnorm(10,3), rnorm(10,3), rnorm(10,3), rnorm(10,6)), # junior:f increases expression
                    GENJ = c( rnorm(10,3), rnorm(10,6), rnorm(10,3), rnorm(10,3)), # senior:f increases expression
                    GENK = c( rnorm(10,3), rnorm(10,3), rnorm(10,3), rnorm(10,3)), # flat around three
                    GENL = c( rnorm(10,4), rnorm(10,4), rnorm(10,4), rnorm(10,4)), # flat around four
                    GENM = c( rnorm(10,5), rnorm(10,5), rnorm(10,5), rnorm(10,5)), # flat around five
                    GENN = c( rnorm(10,6), rnorm(10,6), rnorm(10,6), rnorm(10,6))) # flat around six
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
    object %<>% factorize(k = 2, verbose = verbose)
    object %<>%       bin(k = 2, verbose = verbose)
    object
}


#' @rdname factor.vars
#' @export
left.vars <- function(formula){
    assert_is_formula(formula)
    if (length(formula) == 2) character(0) else all.vars(formula[[2]])
}


#' @rdname factor.vars
#' @export
right.vars <- function(formula){
    assert_is_formula(formula)
    all.vars(formula[[length(formula)]])
}


#' Get factor variables
#' @param formula formula
#' @param object SummarizedExperiment or data.table
#' @return character vector
#' @examples
#' object <- survobj()
#' formula <- survival::Surv(timetoevent, event) ~ age/exprs2levels
#'    all.vars(formula)
#'   left.vars(formula)
#'  right.vars(formula)
#' factor.vars(formula, object)
#' @export
setGeneric('factor.vars',   function(formula, object)  standardGeneric('factor.vars'))

#' @rdname factor.vars
#' @exportMethod factor.vars
setMethod( 'factor.vars', signature(formula = 'formula', object = 'SummarizedExperiment'), 
                            function(formula, object){
                                samplevars <- intersect( right.vars(formula),      svars(object) ) 
                                assayvars  <- intersect( right.vars(formula), assayNames(object) )
                                modeldt <- sumexp_to_longdt(object, svars = samplevars, assay = assayvars, value.name = assayvars)
                                factor.vars(formula, modeldt)
                            }
)

#' @rdname factor.vars
#' @exportMethod factor.vars
setMethod( 'factor.vars', signature(formula = 'formula', object = 'data.table'), 
                            function(formula, object){
                                object %<>% extract(, right.vars(formula), with = FALSE)
                                factorvars <- vapply(object, function(x) is.character(x) | is.factor(x), logical(1))
                                factorvars %<>% extract(.)
                                factorvars %<>% names()
                                factorvars
                            }
)


#' Fit onefeature survival 
#' @param sd               data.table
#' @param formula          model formula
#' @examples
#' # Dataset
#'      sd <- survobj()
#'      sd %<>% sumexp_to_longdt( svars = c('timetoevent', 'event', 'age', 'sex'), assay = 'exprs2levels')
#'      sd[, value := code(factor(value), code_control)]
#'      sd[,   age := code(factor(age  ), code_control)]
#'      sd[,   sex := code(factor(sex  ), code_control)]
#'      
#' # Singlefactor - coxph, survdiff, logrank
#'     .survdiff(sd, survival::Surv(timetoevent, event) ~ value)
#'      .logrank(sd, survival::Surv(timetoevent, event) ~ value)
#'        .coxph(sd, survival::Surv(timetoevent, event) ~ value)
#'        .coxph(sd, survival::Surv(timetoevent, event) ~ age/value)
#' @rdname dot-coxph
#' @export
.coxph <- function(sd, formula){
# coxph
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
         coefs = NULL, 
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
    if (verbose)  cmessage('%sModel %s(%s)', spaces(14), engine, twosideformula) # Align with Code `exprs2levels`
    twosideformula %<>% as.formula()
    if (engine == 'coxph')     fitdt <- dt[,    .coxph(.SD, twosideformula), by = 'feature_id']
    if (engine == 'survdiff')  fitdt <- dt[, .survdiff(.SD, twosideformula), by = 'feature_id']
    if (engine == 'logrank')   fitdt <- dt[,  .logrank(.SD, twosideformula), by = 'feature_id']

    if (drop){ # drop varname from non-numeric vars
        anum <- assays(object)
        snum <- sdt(object)[, samplevars, with = FALSE]
        anum %<>% vapply(is_non_numeric, logical(1))
        snum %<>% vapply(is_non_numeric, logical(1))
        anum <- names(anum)[anum]
        snum <- names(snum)[snum]
        anum %<>% intersect(assayvar)
        snum %<>% intersect(samplevars)
        for (var in c(anum, snum)){  
            pat <- sprintf('%s(.+)', var)
            names(fitdt) %<>% stri_replace_first_regex(pat, '$1')  
    }}
# Select/Return
    if (!is.null(coefs)){  idx <- c(1, which(split_extract_fixed(names(fitdt), '~', 2) %in% coefs))
                           fitdt %<>% extract(, idx, with = FALSE)  }
    if (verbose)  message_df('                      %s', summarize_fit(fitdt))
    fitdt
}


#' @rdname all_non_numeric
#' @export
is_non_numeric <- function(x)  !is.numeric(x)


#' Are all variables non-numeric ?
#' @param object  SummarizedExperiment
#' @param formula formula
#' @param x vector
#' @return TRUE or FALSE
#' @examples
#' all_non_numeric(survobj(), ~ age)
#' all_non_numeric(survobj(), ~ exprs2levels)
#' all_non_numeric(survobj(), ~ age/exprs2levels)
#' all_non_numeric(survobj(), ~ age/exprs)
#' @export
all_non_numeric <- function(object, formula){
    samplevars <- intersect(all.vars(formula),      svars(object))
    assayvars  <- intersect(all.vars(formula), assayNames(object))
    snon <- sdt(object)[, samplevars, with = FALSE]
    anon <- assays(object)[assayvars]
    snon %<>% vapply(is_non_numeric, logical(1))
    anon %<>% vapply(is_non_numeric, logical(1))
    all(c(snon, anon))
}


#' Fit/Plot survival
#' 
#' @param object        SummarizedExperiment
#' @param formula       model formula: contains svars/assayNames
#' @param assaylevels   NULL or vector: assaylevels to be used (for plotting)
#' @param engine       'coxph', 'survdiff' or 'logrank'
#' @param drop          TRUE or FALSE : whether to drop var in coefname
#' @param codingfun     coding function
#' @param coefs         NULL or character (coefs to be stored in object)
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
#' @param dodge         number
#' @param file          filepath
#' @return SummarizedExperiment/ggplot
#' @examples
#' # Formula
#'     # Samplevar-based
#'           fit_survival(survobj(), ~age)           # age
#'           fit_survival(survobj(), ~sex)           # sex
#'           fit_survival(survobj(), ~age + sex)     # age across  sexlevels, sex across agelevels
#'           fit_survival(survobj(), ~age / sex)     # sex within  agelevel
#'           fit_survival(survobj(), ~age * sex)     # sex between agelevels (=age between sexlevels)
#'   
#'     # Assayvar-based
#'           fit_survival(survobj(), ~exprs)         #   numerical coding
#'           fit_survival(survobj(), ~exprs2bins)    #     integer coding
#'           fit_survival(survobj(), ~exprs2levels)  # categorical coding
#'
#'     # Samplevar/Assayvar-based
#'           fit_survival(survobj(), ~age+exprs2levels, order = 'senior-junior'          ) #  age effect across exprlevels
#'           fit_survival(survobj(), ~age+exprs2levels, order = '2-1'                    ) # expr effect across agelevels
#'           fit_survival(survobj(), ~age/exprs2levels, order = 'senior:2-1'             ) # expr effect within agelevel
#'           fit_survival(survobj(), ~age*exprs2levels, order = 'senior-junior:2-1'      ) # expr effect differences between agelevels (or vice versa)
#'   
#' # Other arguments
#'     # engine: 'coxph' -> 'survdiff'
#'           fit_survival(survobj(), ~ exprs2levels)                        # coxph
#'           fit_survival(survobj(), ~ exprs2levels, engine = 'survdiff')   # survdiff
#' 
#'     # drop: drop varname in coefnames -> dont
#'           fit_survival(survobj(), ~ exprs2levels)                # 2-1
#'           fit_survival(survobj(), ~ exprs2levels, drop = FALSE)  # exprs2levels2-1
#' 
#'     # codingfun: code_control -> contr.treatment
#'           fit_survival(survobj(), ~ exprs2levels)                              # code_control
#'           fit_survival(survobj(), ~ exprs2levels, codingfun = contr.treatment) # contr.treatment
#'
#'     # outdir: print to object/screen -> print to xlsx/pdf
#'           fit_survival(survobj(), ~ exprs2levels)                                                 # print to object/screen
#'           fit_survival(survobj(), ~ exprs2levels, outdir = tempdir())                             # print to   xlsx/pdf
#'           fit_survival(survobj(), ~ exprs2levels, outdir = tempdir(), writefunname = 'write_ods') # print to    ods/pdf
#' 
#'     # plot: plot -> dont
#'           fit_survival(survobj(), ~ exprs2levels)                # plot
#'           fit_survival(survobj(), ~ exprs2levels, plot = FALSE)  # dont
#' 
#'     # order: order on first coef -> order on custom coef
#'           fit_survival(survobj(), ~ age+exprs2levels)                  # order on 'senior-junior'
#'           fit_survival(survobj(), ~ age+exprs2levels, order = '2-1')   # order on '2-1'
#' 
#'     # stats: show stats for all coefs -> show stats for custom coefs
#'           fit_survival(survobj(), ~ age+exprs2levels)                          # show stats for 'senior-junior' and 'bin2-bin1'
#'           fit_survival(survobj(), ~ age+exprs2levels, stats = 'senior-junior') # show stats for 'senior-junior'
#' 
#'     # dodge: overlap curves -> dodge curves
#'           fit_survival(survobj(), ~ age+exprs2levels)            # overlap curves
#'           fit_survival(survobj(), ~ age+exprs2levels, dodge = 2) # dodge curves
#' 
#'     # n: (plot) top2 -> top4
#'           fit_survival(survobj(), ~ age+exprs2levels)         # top2
#'           fit_survival(survobj(), ~ age+exprs2levels, n = 4)  # top4
#' 
#'     # n_row n_col: 1 row 2 col -> 2 row 1 col
#'           fit_survival(survobj(), ~ age+exprs2levels)                       # 1 row 2 col
#'           fit_survival(survobj(), ~ age+exprs2levels, n_row = 2, n_col = 1) # 2 row 1 col
#' @export
fit_survival <- function(
        object, 
       formula = as.formula(sprintf('~%s', assayNames(object)[1])),
        engine = c('coxph', 'survdiff', 'logrank')[1],
          drop = TRUE,
     codingfun = code_control,
         coefs = NULL,
       verbose = TRUE,
        outdir = NULL,
          plot = FALSE,
         order = coefs(object, fit = engine)[1],
         stats = coefs(object, fit = engine),
         dodge = 0,
             n = if (svar_formula(formula, object)) 1  else min(nrow(object),2), # Inf works
         n_col = n %>% min(nrow(object)) %>% sqrt() %>% ceiling() %>% min(4),
         n_row = n %>% min(ncol(object)) %>% sqrt() %>% floor()   %>% min(4),
         width = 3*n_col,       #  Only for formula group is sample property (e.g. sex) sharing guaranteed
        height = 3*n_row,
  writefunname = 'write_xl'
){
    if (verbose)  cmessage('%sSurvival', spaces(4))
# Compute
    for (eng in engine){
        fitdt <- .fit_survival(  object = object, 
                                formula = formula,
                                  coefs = coefs,
                                 engine = engine,
                                   drop = drop,
                              codingfun = codingfun,
                                verbose = verbose )
        if (all(all.vars(formula) %in% svars(object))){  metadata(object)$survival <-  fitdt
        } else {                                         object %<>% merge_fdt(fitdt)  }
    }
# Write
    if (!is.null(outdir)){
        cmessage('%sPrint', spaces(14))
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
                                 order = order,
                                 stats = stats,
                                 dodge = dodge,
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
 assaylevels = NULL, 
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
    object %<>% extract_contrast_features(fit = engine, coefs = order, n = n)
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
        if (!is.null(assaylevels))  plotdt <- plotdt[ get(assayvar) %in% assaylevels ]
    }
    plotdt[, alive := 1-event]
    setorderv(plotdt, c('feature_id', all.vars(formula), 'timetoevent', 'alive'))
    
    plotdt[ , totObs   := .N - cumsum(1-event),     by = c('feature_id', all.vars(formula))   ]
    plotdt[ , totDead := cumsum(event),             by = c('feature_id', all.vars(formula))   ]
    
# Aggregate per timepoint
    plotdt <- plotdt[ , .(totObs  = max(totObs), 
                          totDead = max(totDead), 
                          curOut  = sum(event==0)), by = c('feature_id', all.vars(formula), 'timetoevent')]
    
# Survival = delayed(alive/(alive+dead))
    plotdt[, survival := 100*(totObs-totDead)/totObs]
    plotdt[curOut!=0, survival := NA_real_]
    plotdt[, survival := data.table::nafill(survival, type = 'locf'), by = c('feature_id', all.vars(formula)) ]
    
# Start with 100% survival
    setorderv(plotdt, c('feature_id', all.vars(formula), 'timetoevent'))
    plotdt0 <- plotdt[ , .SD[ 1] , by = c('feature_id', all.vars(formula))][, timetoevent := 0 ][, totDead := 0 ][, survival := 100][, curOut := 0]
    plotdt <- rbind(plotdt0, plotdt)
    setorderv(plotdt, c('feature_id', all.vars(formula), 'timetoevent'))
    plotdt
    
# Vertically end curve when all dead
    plotdtn <- plotdt[ , .SD[.N] , by = c('feature_id', all.vars(formula))][, timetoevent := max(timetoevent)+1][, curOut := 0]
    plotdtn <- plotdtn[totDead!=totObs]
    features <- plotdt[, levels(feature_id)]
    plotdt <- rbind(plotdt, plotdtn)            # Preserve tvar order
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


#' Compute step auc
#' @param x    numeric vector
#' @param y    numeric vector
#' @param plot TRUE or FALSE
#' @return number
#' @examples
#' x <- c(  0, 4,   8, 27)
#' y <- c(100, 67, 33,  0)
#' stepauc(x, y, plot = TRUE)
#' @export
stepauc <- function(x, y, color = 'group1', plot = FALSE){
        ord <- order(x)
          x <- x[ord]
          y <- y[ord]
         dx <- diff(x)
    heights <- tail(y, -1)
    auc <- sum(dx*heights)
    if (plot){   p <- ggplot(data = data.table(x = x, y = y, color = color))
                 p <- p + theme_bw() + theme(panel.grid = element_blank())
                 p <- p + geom_step(   aes(x = x, y = y, color = color))
                 p <- p + geom_rect(aes(xmin = x, xmax = dplyr::lead(x), ymin = 0, ymax = y,  fill = color), alpha = 0.2)
                 p <- p + geom_text(aes(x = min(x), y = min(y), label = auc, color = color), hjust = -3, vjust = -3)
                 print(p)  }
    auc
}



#' @rdname fit_survival
#' @export
plot_survival <- function(
      object,
     formula = as.formula(sprintf('~%s', assayNames(object)[1])), 
 assaylevels = NULL,
      engine = c('coxph', 'survdiff', 'logrank') %>% intersect(fits(object)) %>% extract(1),
       order = autonomics::coefs(object, fit = engine)[1],
       stats = autonomics::coefs(object, fit = engine),
       title = sprintf('%s ~ %s', engine, formula2str(formula) %>% substr(2,nchar(.))),
       dodge = 0,      # `color` and `linetype` are hardmapped from `all.vars(formula)`
        file = NULL,    #  softmapping them formula-agnostically doesnt work
           n = if (svar_formula(formula, object)) 1  else min(nrow(object),4), # Inf works
       n_col = n %>% min(nrow(object)) %>% sqrt() %>% ceiling() %>% min(4),
       n_row = n %>% min(ncol(object)) %>% sqrt() %>% floor()   %>% min(4),
       width = 3*n_col,       #  Only for formula group is sample property (e.g. sex) sharing guaranteed
      height = 3*n_row
){
# Assert
    if (!installed('ggtext'))   return(NULL) 
    if (!installed('ggstance')) return(NULL)
    if (is.infinite(n))  n <- nrow(object)
    totObs <- totDead <- nalive <- nout <- label <- color <- facet <- NULL
    timetoevent <- survival <- NULL
# Plot
    plotdt <- prep_survival(object = object, formula = formula, assaylevels = assaylevels, 
                            engine = engine, order = order, stats = stats, n = n)
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
    if (!is.null(file))  cmessage('%s%s', spaces(21), file)
    if (!is.null(file))  pdf(file, width = width, height = height)
    for (i in seq_len(npages)){
        subtitle <- if (svar_formula(formula, object)) NULL else paste0(order, collapse = '  ')
        # Steps
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
             colorvar <-      all.vars(formula)[[1]]
             colorsym <- sym( all.vars(formula)[[1]])
             alphavar <- if (length(all.vars(formula))<2)     NULL  else      all.vars(formula)[[2]]
             alphasym <- if (length(all.vars(formula))<2) quo(NULL) else sym( all.vars(formula)[[2]])
          alphalevels <- if (length(all.vars(formula))<2)     NULL  else seq( from = 0.3, to = 1, length.out = length(unique(plotdt[[all.vars(formula)[[2]]]])) )
            p <- p + geom_step( mapping = aes(  x = timetoevent, 
                                                y = survival,              
                                            group = interaction(!!!groupsyms),  # !!! for syms
                                            color = !!colorsym,                 #  !! for sym
                                            alpha = !!alphasym ) ,           # position_identity speedsup code 2.5 times
                                         position = if (dodge == 0) position_identity() else ggstance::position_dodgev(dodge))
            p <- p + scale_alpha_manual(values = alphalevels)
        
        # Pluses (censorings)
            p <- p + geom_point(data = plotdt[curOut==1], mapping = aes(x = timetoevent, y = survival, color = !!colorsym, alpha = !!alphasym), size = 1.5)
            p <- p + geom_point(data = plotdt[curOut >1], mapping = aes(x = timetoevent, y = survival), size = 2)
            
        # Labels
            labeldt <- plotdt[ , .( x = max(timetoevent), 
                                    y = max(survival), 
                                  auc = stepauc(timetoevent,survival),
                                label = sprintf('%s%s: %s*  %s†  %s•', 
                                                              get(colorvar), 
                                                              if (is.null(alphavar)) '' else  paste0('.', get(alphavar)), 
                                                              totObs[.N]-totDead[.N], 
                                                              totDead[.N],
                                                              totObs[1] - totObs[.N]
                                                )) , by = c('facet', colorvar, alphavar)]
            labeldt[, x := max(x)]
            labeldt <- labeldt[, .SD[rev(order(auc))], by = 'facet']
            labeldt[, i := seq(0,.N-1) , by = 'facet']
            labeldt[, y := (1-i*0.12)*y ]
            p <- p + geom_text(data = labeldt, mapping = aes(x = x, y = y, label = label, color = !!colorsym, alpha = !!alphasym), hjust = 1, vjust = 1, show.legend = FALSE)
           #p <- p + guides(color = 'none', alpha = 'none')
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

