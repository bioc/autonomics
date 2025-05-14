#==============================================================================
#
#                           rm_unmatched_samples
#                           rm_singleton_samples
#
#==============================================================================

#' rm unmatched/singleton samples
#'
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup variable (string)
#' @param subgroupctr  control subgroup (string)
#' @param block        block variable (string)
#' @param verbose      TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- system.file('extdata/atkin.somascan.adat', package = 'autonomics')
#' object <- read_somascan(file)
#' object %<>% filter_samples(subgroup %in% c('t1', 't2'), verbose = TRUE)
#' rm_singleton_samples(object, subgroupvar = 'Subject')
#' rm_unmatched_samples(object, subgroupvar = 'subgroup', block = 'Subject')
#' @export
rm_unmatched_samples <- function(
    object,
    subgroupvar = 'subgroup',
    subgroupctr = slevels(object, subgroupvar)[1],
    block,
    verbose = TRUE
){
    snames1 <- sdt(object)[,
                .SD[(sum(get(subgroupvar)==subgroupctr)==1) &
                    (sum(get(subgroupvar)!=subgroupctr) >0)],  by = block]$sample_id
    n <- length(snames1)
    if (verbose & n < ncol(object)){
        message('\t\tRetain ', n, '/', ncol(object), ' samples with matching ', subgroupctr) }
    object %<>% extract(, snames1)
    object
}

#' @rdname rm_unmatched_samples
#' @export
rm_singleton_samples <- function(object, subgroupvar = 'subgroup', verbose = TRUE){
    selectedsamples <- sdt(object)[, .SD[.N>1], by = subgroupvar][['sample_id']]
    n <- length(selectedsamples)
    if (verbose & n < ncol(object)){
        message('\t\tRetain ', length(selectedsamples), '/',
                ncol(object), ' samples with replicated ', subgroupvar) }
    object[, selectedsamples]
}


#==============================================================================
#
#                           log2transform()
#                           zscore()
#                           quantnorm()
#                           invnorm()
#
#==============================================================================

# mat <- cbind(s1=c(-1,-1), s2=c(-1,1), s3=c(1,-1), s4=c(0.1,0.1))
# which.medoid(mat)
which.medoid <- function(mat){
    idx <- matrixStats::rowAnyNAs(mat)
    assert_all_are_not_na(idx)
    if (any(idx))  cmessage('\t\t\t\tUse %d/%d non-NA rows to compute spatial median', 
                           sum(idx), length(idx))
    mat %<>% extract(!idx, )

    tryCatch(
        {   spatmed <- ICSNP::spatial.median(t(mat))
            which.min(sqrt(colSums((sweep(mat, 1, spatmed))^2)))
        }, 
        error = function(cond){
            message('    spatial median failed - using simple centroid instead')
            centroid <- colMeans(mat)
            centroid %<>% matrix(nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
            which.min(colSums((mat - centroid)^2))
        }
    )
}

.filter_medoid <- function(object, verbose = FALSE){
    if (ncol(object)==1)  return(object)
    medoid <- which.medoid(values(object))
    object %<>% extract(, medoid)
    if (verbose)  message('\t\t\t\t', object$sample_id)
    object
}

#' Filter medoid sample
#' @param object SummarizedExperiment
#' @param by svar
#' @param verbose  whether to message
#' @return SummarizedExperiment
#' @examples 
#' file <- system.file('extdata/billing19.rnacounts.txt', package = 'autonomics')
#' object <- read_rnaseq_counts(file, plot = FALSE)
#' object %<>% filter_medoid(by = 'subgroup', verbose=TRUE)
#' @export
filter_medoid <- function(object, by = NULL, verbose = FALSE){
    if (!installed('ICSNP'))    return(object)
    if (is.null(by))  return(.filter_medoid(object, verbose=verbose))
    object %<>% split_samples(by)
    if (verbose)  message('\t\t\tRetain medoid sample')
    object %<>% lapply(.filter_medoid, verbose=verbose)
    do.call(BiocGenerics::cbind, object)
}

.subtract_baseline <- function(
    object, subgroupvar, subgroupctr, assaynames
){
    idx <- object[[subgroupvar]] == subgroupctr
    controls <- object[, idx]
    perturbs <- object[, !idx]
    controls %<>% filter_medoid()
    for (ass in assaynames){
        assays(perturbs)[[ass]] %<>% sweep(1, assays(controls)[[ass]]) 
        # Confirm that sweep works as thought (it does)
        # controlmat <- assays(controls)[[ass]]
        # assert_is_identical_to_true(ncol(controlmat)==1)
        # controlmat %<>% extract(, 1)
        # controlmat %<>% matrix(byrow = FALSE, nrow = length(controlvec), ncol = ncol(perturbs))
        # assays(perturbs)[[ass]] %<>% subtract(controlmat)
    }
    perturbs
}

#' Subtract baseline
#' 
#' Subtract baseline level within block
#' 
#' \code{subtract_baseline} subtracts baseline levels within block, using the 
#' medoid baseline sample if multiple exist. \cr
#' 
#' \code{subtract_pairs} also subtracts baseline level within block. 
#' It cannot handle multiple baseline samples, but has instead been optimized
#' for many blocks \cr
#' 
#' \code{subtract_differences} subtracts differences between subsequent levels, 
#' again within block
#' @param  object       SummarizedExperiment
#' @param  subgroupvar  subgroup svar
#' @param  subgroupctr  control subgroup
#' @param  block        block svar (within which subtraction is performed)
#' @param  assaynames   which assays to subtract for
#' @param  verbose      TRUE/FALSE
#' @return SummarizedExperiment
#' @examples
#' # read 
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object0 <- read_metabolon(file)
#'     pca(object0, plot = TRUE, color = 'Time')
#' 
#' # subtract_baseline: takes medoid of baseline samples if multiple
#'     object <- subtract_baseline(object0, block = 'Subject', subgroupvar = 'Time')
#'     pca(object, plot = TRUE, color = 'Time')
#' 
#' # subtract_pairs: optimized for many blocks
#'     object <- subtract_pairs(object0, block = 'Subject', subgroupvar = 'Time')
#'     pca(object, plot = TRUE, color = 'Time')
#' 
#' # subtract_differences
#'     object <- subtract_differences(object0, block = 'Subject', subgroupvar = 'Time')
#'     values(object) %<>% na_to_zero()
#'     pca(object, plot = TRUE, color = 'Time')
#' @export 
subtract_baseline <- function(
    object, subgroupvar, subgroupctr = slevels(object, subgroupvar)[1], 
    block = NULL, assaynames = setdiff(assayNames(object), c('weights', 'pepcounts')), 
    verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(subgroupvar, svars(object))
    assert_scalar_subset(subgroupctr, unique(object[[subgroupvar]]))
    if (!is.null(block))  assert_scalar_subset(block, svars(object))
    assert_is_subset(assaynames, assayNames(object))
    assert_is_a_bool(verbose)
# Subtract
    if (verbose){ 
        message("\t\tSubtract controls")
        message("\t\t\tcontrols  : ", subgroupvar, "=", subgroupctr," (medoid)")
        if (!is.null(block))  message("\t\t\tin block  : ", block)
        message("\t\t\tfor assays: ", paste0(assaynames, collapse = ', '))
    }
    objlist <- if (is.null(block)) list(object) else split_samples(object, block)
    objlist %<>% lapply(.subtract_baseline, 
                        subgroupvar = subgroupvar, subgroupctr = subgroupctr, 
                        assaynames = assaynames)
# Return
    objlist %<>% do.call(BiocGenerics::cbind, .)
    objlist
}


#' @rdname subtract_baseline
#' @export
subtract_pairs <- function(
    object, 
    subgroupvar = 'subgroup', 
    subgroupctr = slevels(object, subgroupvar)[1], 
    block,
    assaynames = assayNames(object)[1], verbose = TRUE
){
# Report
    # PRO: optimized for many block levels
    # CON: works only with exactly one ref per block 
    . <- NULL
    if (verbose){ 
        txt <- paste0("\t\tSubtract ", subgroupvar, "==", subgroupctr, ' ')
        txt %<>% paste0(paste0(assaynames, collapse = '/'))
        if (!is.null(block))  txt %<>% paste0(" per ", block)
        message(txt)
    }
# Ensure single ref per block
    sdt1 <- sdt(object)[, c('sample_id', subgroupvar, block), with = FALSE]
    singlerefperblock <- sdt1[, sum(get(subgroupvar)==subgroupctr)==1, by=block]$V1
    assert_is_identical_to_true(all(singlerefperblock))
# Subtract ref
    splitobjects <- split_samples(object, subgroupvar)
    refobj <- splitobjects[[subgroupctr]]
    splitobjects %<>% extract(setdiff(names(splitobjects), subgroupctr))
    splitobjects %<>% lapply(function(obj){
        idx <- match(obj[[block]], refobj[[block]])
        refobj %<>% extract(, idx)
        assert_is_identical_to_true(all(obj[[block]] == refobj[[block]]))
        for (assayname in assaynames){
            assays(obj)[[assayname]] %<>% subtract(assays(refobj)[[assayname]])
            assayNames(obj) %<>% stri_replace_first_fixed(
                                    assayname, paste0(assayname, 'ratios'))
        }
        obj        
    })
    splitobjects %<>% do.call(S4Vectors::cbind, .)
    idx <- na.exclude(match(object$sample_id, splitobjects$sample_id))
    splitobjects[, idx]
}


#' @rdname subtract_baseline
#' @export
subtract_differences <- function(object, block, subgroupvar, verbose=TRUE){
    # PRO: robust (works with missing levels, or multiple replicates per level)
    # CON: performance not optimized for many block levels
    #      currently only performed on first assay (can off course be updated)
    sample_id <- NULL
    if (verbose){ 
        message("\t\tSubtract differences")
        if (!is.null(block))  message("\t\t\tin block  : ", block)
        message("\t\t\tfor assays: ", assayNames(object)[1])
    }
    fvars0 <- intersect(c('feature_id', 'feature_name'), fvars(object))
    dt <- sumexp_to_longdt(object, fvars=fvars0, svars=c(subgroupvar, block))
    subgroups <- slevels(object, subgroupvar)
    n <- length(subgroups)
    formula <- paste0(c(fvars0, block), collapse = ' + ')
    formula %<>% paste0(' ~ ', subgroupvar)
    formula %<>% as.formula()
    dt %<>% dcast(formula, value.var ='value')
    
    newdt  <- dt[, c(fvars0, block), with = FALSE]
    diffdt <- dt[, setdiff(names(dt), c(fvars0, block)), with=FALSE]
    diffdt  <-  diffdt[, subgroups[-1], with=FALSE] - 
                diffdt[, subgroups[-n], with=FALSE]
    names(diffdt) %<>% paste0('_', subgroups[-n])
    newdt %<>% cbind(diffdt)
    
    newdt %<>% melt.data.table(
                id.vars =  c(fvars0, block), variable.name = subgroupvar)
    data.table::setorderv(newdt, c('feature_id', block, subgroupvar))
    newdt[, sample_id := paste0(get(block), '.', get(subgroupvar))]
    newobject <- dt2sumexp(newdt)
    assayNames(newobject)[1] <- assayNames(object)[1]
    newobject
}


#' Transform values
#' 
#' @param  object   SummarizedExperiment
#' @param  mat      matrix
#' @param  assay    character vector : assays for which to perform transformation
#' @param  pseudo   number           : pseudo value to be added prior to transformation
#' @param  verbose  TRUE or FALSE    : whether to msg
#' @param  delog    TRUE or FALSE (vsn)
#' @param  relog    TRUE or FALSE (vsn)
#' @return Transformed sumexp
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#'
#' object                                %>% plot_sample_densities()
#' invnorm(object)                       %>% plot_sample_densities()
#'
#' object                                %>% plot_sample_densities()
#' quantnorm(object)                     %>% plot_sample_densities()
#'
#' object                                %>% plot_sample_densities()
#'#vsn(object)                           %>% plot_sample_densities()  # dataset too small
#'
#' object                                %>% plot_sample_densities()
#' zscore(object)                        %>% plot_sample_densities()
#'
#' object                                %>% plot_sample_densities()
#' exp2transform(object)                 %>% plot_sample_densities()
#' log2transform(exp2transform(object))  %>% plot_sample_densities()
#' @export
log2transform <- function(
    object, 
    assay   = assayNames(object)[1], 
    pseudo  = 0,
    verbose = FALSE
){
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_not_null(assayNames(object))
    assert_is_subset(assay, assayNames(object))
    . <- NULL
    for (ass in assay){
        i <- match(ass, assayNames(object))
        if (verbose)  if (pseudo != 0)  message('\t\tAdd pseudo value ', pseudo)
        assays(object)[[i]] %<>% add(pseudo) 
         
        if (verbose)  cmessage('%slog2 %s', spaces(14), ass)
        assays(object)[[i]] %<>% log2()
        assayNames(object)[i] %<>% paste0('log2', .)
    }
    object
}


#' @rdname log2transform
#' @export
exp2transform <- function(
    object, 
    assay   = assayNames(object)[1], 
    verbose = FALSE
){
  assert_is_all_of(object, 'SummarizedExperiment')
  assert_is_not_null(assayNames(object))
  assert_is_subset(assay, assayNames(object))
  . <- NULL
  for (ass in assay){
    i <- match(ass, assayNames(object))
    if (verbose)  cmessage('%sexp2 %s', spaces(14), ass)
    assays(object)[[i]] %<>% magrittr::raise_to_power(2, .)
    if (stri_startswith_fixed(assayNames(object)[i], "log2")) {
      assayNames(object)[i] %<>% stri_replace_first_regex('^log2\\s*', '')
    } else {
      assayNames(object)[i] %<>% paste0('log2', .)
    }
  }
  object
}

#' @rdname log2transform
#' @export
zscore <- function(object, verbose = FALSE){
    values(object) %<>% sscale(verbose = verbose)
    object
}

#' @rdname log2transform
#' @export
sscale <- function(mat, verbose = FALSE){
    assert_is_matrix(mat)
    if (verbose)  message('\t\tZscore samples')
    scale(mat)
}

#' @rdname log2transform
#' @export
fscale <- function(mat, verbose = FALSE){
    assert_is_matrix(mat)
    if (verbose)  message('\t\tZscore features')
    t(scale(t(mat)))
}

#' Center samples
#' @param object   SummarizedExperiment
#' @param selector logical vector (length = nrow(object))
#' @param fun      aggregation function (string)
#' @param verbose  TRUE/FALSE
#' @param ...      parameters handed through to center()
#' @return SummarizedExperiment
#' @examples
#' require(matrixStats)
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' fdt(object)$housekeeping <- FALSE
#' fdt(object)$housekeeping[order(rowVars(values(object)))[1:5]] <- TRUE
#' values(object)[, object$subgroup=='Adult'] %<>% magrittr::add(5)
#' plot_sample_densities(object)
#' plot_sample_densities(center(object))
#' plot_sample_densities(center(object, housekeeping))
#' @export
center <- function(
    object, selector = rep(TRUE, nrow(object))==TRUE, fun = 'median', verbose = TRUE
){
    selector <- enexpr(selector)
    selector <- rlang::eval_tidy(selector, data = fdata(object))
    if (verbose)  message('\t\t', fun, ' center samples on ', 
                            nrow(object[selector, ]), ' features')
    correction_factors <- apply(values(object[selector, ]), 2, fun, na.rm=TRUE)
    correction_factors[is.na(correction_factors)] <- 0
    values(object) %<>% sweep(2, correction_factors)
    object
}

#' @rdname center
#' @export
center_mean <- function(object, ...)
{
  object %<>% center(fun = 'mean', ...)
}

#' @rdname center
#' @export
center_median <- function(object, ...)
{
  object %<>% center(fun = 'median', ...)
}

#' @rdname log2transform
#' @export
quantnorm <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tQuantnorm')
    values(object) %<>% limma::normalizeBetweenArrays()
    object
}


#' @rdname log2transform
#' @export
invnorm <- function(object, verbose = FALSE){
    if (verbose)  message('\t\tInvnorm')
    values(object) %<>% apply(2, transform_to_fitting_normal)
    object
}

#' @rdname log2transform
#' @export
vsn <- function(object, delog = TRUE, relog = delog, verbose = FALSE)
{
    assert_is_valid_sumexp(object)
    assert_is_a_bool(delog)
    assert_is_a_bool(relog)
    assert_is_a_bool(verbose)
    if (verbose) message('\t\tVSN')
    if (delog) object %<>% exp2transform(verbose = FALSE)
    values(object) %<>% vsn::justvsn(verbose = FALSE)
    if (relog) object %<>% log2transform(verbose = FALSE)
    object
}

#' Transform vector to fitting normal distribution
#' @param x numeric vector
#' @return transformed vector
#' @noRd
transform_to_fitting_normal <- function(x){
    pars <- estimate_mean_sd(x)
    transform_to_normal(x, mean = pars[['mean']], sd = pars[['sd']])
}

estimate_mean_sd <- function(x){
    . <- NULL
    x %<>% extract(!is.na(.) & !is.infinite(.))
    MASS::fitdistr(x, 'normal')[['estimate']]
}



#' Transform vector to normal distribution
#' @param x numeric vector
#' @param mean  mean
#' @param sd    standard deviation
#' @return transformed vector
#' @noRd
transform_to_normal <- function(x, mean, sd){
    selector <- !is.na(x) & !is.nan(x) & !is.infinite(x)
    pvals <- rank(x[selector]) / (length(x[selector]) + 1)
    y <- x
    y[selector] <- qnorm(pvals, mean = mean, sd = sd)
    y
}


#' Transform vector to standard normal distribution
#' @param x numeric vector
#' @return transformed vector
#' @noRd
transform_to_standard_normal <- function(x){
    transform_to_normal(x, mean = 0, sd = 1)
}


#==============================================================================
#
#                        biplot_transforms
#                        biplot_transforms_assays
#                        plot_densities_transforms
#                        plot_violins_transforms
#
#==============================================================================

gglegend<-function(p){
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(vapply(
                    tmp$grobs, function(x) x$name, character(1))=="guide-box")
    if (length(leg)==0)  grid::nullGrob()  else  tmp$grobs[[leg]]
}

#' Visually evaluate transformation effects
#' 
#' @param  object       SummarizedExperiment
#' @param  assay        string           : assay name to operate on
#' @param  assays       character vector : assay names to operate on
#' @param  subgroupvar  svar
#' @param  transforms   character vector : transformations explored
#' @param  method       string           : dimension reduction technique
#' @param  dims         numbers          : biplot dimensions
#' @param  color        svar
#' @param  shape        svar
#' @param  size         svar
#' @param  alpha        svar
#' @param  group        svar
#' @param  label        svar
#' @param  ncol         integer          : columns for facet wraping
#' @param  nrow         integer          : rows for facet wraping
#' @param  ...                           : further plotting parameters
#' @param  fixed        list             : fixed  aesthetics
#' @param  verbose      TRUE/FALSE       : message?
#' @return ggplot2 object
#' @author Johannes Graumann
#' @rdname explore-transforms
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#'
#' # `vsn` implemented, but example data set to small
#' transformations <- c(
#'   'center_mean', 'center_median', 'invnorm', 'quantnorm', 'zscore')
#'
#' # object %>% plot_densities_transforms(transforms = transformations) # Requires package ggridges
#' object %>% plot_violins_transforms(transforms = transformations)
#' 
#' object %>% biplot_transforms(
#'   method  = 'pca', transforms = transformations, nrow = 2)
#' object %>% biplot_transforms(
#'   method  = 'pls', transforms = transformations, nrow = 2)
#'   
#' object[['replicate']] <- gsub('^.*\\.(.+)$', '\\1', object[['sample_id']])
#' object %>%
#'   biplot_transforms(
#'   transforms = transformations, label = 'replicate')
#' @author Johannes Graumann
#' @export
plot_densities_transforms <- function(
    object,
    assay       = assayNames(object)[1],
    subgroupvar = 'subgroup',
    transforms  = c('center', 'invnorm', 'quantnorm', 'vsn' , 'zscore'),
    ...,
    fixed       = list(na.rm = TRUE, show.legend = FALSE, verbose = FALSE),
    verbose     = TRUE
){
    if (!requireNamespace('ggridges', quietly = TRUE)){
      message("`BiocManager::install('ggridges')`. Then re-run.")
      return(NULL)
    }
    . <- transfo <- NULL
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, assayNames(object))
    assert_scalar_subset(subgroupvar, svars(object))
    assert_is_subset(
      transforms,
      c('center', 'center_mean', 'center_median', 'invnorm', 'quantnorm', 'vsn',
        'zscore'))
    assert_is_a_bool(verbose)
    value <- sample_id <- NULL

    dt <- .ldt_transforms(
      object, assay, transforms, subgroupvar, verbose = verbose)

    plot_data(dt, ggridges::geom_density_ridges, x = value, y = sample_id,
            color = NULL, fill = !!sym(subgroupvar), ..., fixed = fixed) +
      facet_grid(
        rows   = vars(!!sym(subgroupvar)),
        cols   = vars(transfo),
        scales = "free") +
      labs(title = paste("Assay:", assay))
}

#' @rdname explore-transforms
#' @author Johannes Graumann
#' @export
plot_violins_transforms <- function(
    object,
    assay       = assayNames(object)[1],
    subgroupvar = 'subgroup',
    transforms  = c('center', 'invnorm', 'quantnorm', 'vsn', 'zscore'),
    ...,
    fixed       = list(
      na.rm=TRUE, trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75),
      show.legend = FALSE),
    verbose     = TRUE
){
    . <- transfo <- NULL
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, assayNames(object))
    assert_scalar_subset(subgroupvar, svars(object))
    assert_is_subset(
      transforms,
      c('center', 'center_mean', 'center_median', 'invnorm', 'quantnorm', 'vsn',
        'zscore'))
    assert_is_a_bool(verbose)
    
    value <- sample_id <- NULL
    
    dt <- .ldt_transforms(
      object, assay, transforms, subgroupvar, verbose = verbose)
    
    plot_data(dt, geom_violin, x = sample_id, y = value, 
                         color = NULL, fill = !!sym(subgroupvar), ...,
                         fixed = fixed) +
      facet_grid(
        rows   = vars(!!sym(subgroupvar)),
        cols   = vars(transfo),
        scales = "free") +
      coord_flip() +
      labs(title = paste("Assay:", assay))
}

#' @author Johannes Graumann
.ldt_transforms <- function(
    object,
    assay = assayNames(object)[1], transforms, subgroupvar, verbose = TRUE)
{
  dt <- lapply(
    c('input', transforms),
    function(tf)
    {
      tmpdt <- if (tf == 'input') {
        sumexp_to_longdt(object, assay = assay, svars = subgroupvar)
      } else {
        tmpobject <- object
        assays(tmpobject) %<>% extract(assay)
        sumexp_to_longdt(
          get(tf)(tmpobject, verbose = verbose),
          assay = assay, svars = subgroupvar)
        }
      tmpdt$transfo <- tf
      tmpdt
    }) %>%
    rbindlist()
  dt$transfo %<>% factor(unique(.))
  dt$assay <- assay
  dt
}

#' @author Johannes Graumann
.ldt_transforms_dimred <- function(
    object, assay = assayNames(object)[1],
    transforms, method, methodname, dims = 1:2,
    verbose = TRUE)
{
    annot <- transfo <- transfoannot <- NULL
    
    # Process each transformation and return a merged data.table
    result_data <- lapply(
      c('input', transforms),
      function(transform) {
        tmpobj <- object
        assays(tmpobj) <- assays(tmpobj)[
          c(assay, setdiff(assayNames(tmpobj), assay))]
        
        # Apply transformation if it's not the 'input' transformation
        if (transform != 'input') tmpobj %<>% get(transform)(verbose = verbose)
        
        # Perform requested dimensional reduction
        tmpobj %<>% get(method)(dims = dims, verbose = verbose)
        variances <- variances(tmpobj, method = method, dims = dims)
        
        # Process the transformed data
        transformed_data <- sdt(tmpobj)
        transformed_data$variance.annot <- if (method %in% DIMREDUN) {
          sprintf('%d & %d%%', variances[1], variances[2])
        } else { sprintf('%d%%', variances[1]) }
        transformed_data$variance.annot.ext <- paste0(
          transform, ': ', transformed_data$variance.annot)
        transformed_data$transform = transform
        transformed_data$assay = assay
        
        return(transformed_data)
      }) %>%
      rbindlist(use.names = TRUE, fill = TRUE)
  
    # Convert relevant columns to factors with unique levels
    result_data$assay %<>% factor(unique(.))
    result_data$transform %<>% factor(unique(.))
    result_data$variance.annot.ext %<>% factor(unique(.))
    result_data$variance.annot %<>% factor(unique(.))
    
    return(result_data)
}

.dimredstrings <- function(object, method, dims)
{
    tmp_object <- get(method)(object, dims = dims, verbose = FALSE)
    by         <- biplot_by(tmp_object, method)[1]
    methodname <- methodname(method, by)
    sep        <- guess_fitsep(fdt(tmp_object))
    x          <- scorenames(method, by = by, dims = dims[[1]], sep = sep)
    y          <- scorenames(method, by = by, dims = dims[[2]], sep = sep)
    c(
      by         = by,
      methodname = methodname,
      sep        = sep,
      x          = x,
      y          = y)
}

#' @rdname explore-transforms
#' @author Johannes Graumann
#' @export
biplot_transforms <- function(
    object,
    assay       = assayNames(object)[1],
    subgroupvar = 'subgroup',
    transforms  = TRANSFORMSTRICT,  # without 'center_mean' & 'center_median'
    method      = DIMREDENGINES[1], # 'pca'
    dims        = 1:2,
    color       = subgroupvar, 
    shape       = NULL, 
    size        = NULL, 
    alpha       = NULL,
    group       = NULL,
    label       = NULL,
    ncol        = NULL,
    nrow        = NULL,
    ... ,
    fixed       = list(shape = 15, size = 3),
    verbose     = FALSE
){
    variance.annot.ext <- NULL
    
    assert_is_valid_sumexp(object)
    assert_scalar_subset(assay, setdiff(assayNames(object), "imputed"))
    assert_scalar_subset(subgroupvar, svars(object))
    assert_is_subset(transforms, TRANSFORMENGINES)
    assert_scalar_subset(method, DIMREDENGINES)
    assert_is_numeric(dims)
    assert_are_same_length(dims, numeric(2))
    assert_is_a_bool(verbose)

    drstrings <- .dimredstrings(object, method, dims)

    scoredt <- .ldt_transforms_dimred(
      object, assay, transforms, method, drstrings['methodname'], dims, verbose)
    
    .sdata_biplot(
      scoredt,
      x        = drstrings['x'],
      y        = drstrings['y'],
      method   = method, 
      dims     = dims,
      color    = color,
      shape    = shape,
      size     = size,
      alpha    = alpha,
      group    = group,
      label    = label,
      fixed    = fixed,
      ... ) +
      facet_wrap(
        vars(variance.annot.ext), ncol = ncol, nrow = nrow, scales = "free") +
      labs(title = paste("Assay:", assay))
}

#' @rdname explore-transforms
#' @author Johannes Graumann
#' @export
biplot_transforms_assays <- function(
    object,
    assays      = assayNames(object)[1],
    subgroupvar = 'subgroup',
    transforms  = TRANSFORMSTRICT,  # without 'center_mean' & 'center_median'
    method      = DIMREDENGINES[1], # 'pca'
    dims        = 1:2,
    color       = subgroupvar, 
    shape       = NULL, 
    size        = NULL, 
    alpha       = NULL,
    group       = NULL,
    label       = NULL,
    ... ,
    verbose     = FALSE,
    fixed       = list(shape = 15, size = 3)
){
    assay <- variance.annot <- NULL
    
    assert_is_valid_sumexp(object)
    assert_is_subset(assays, setdiff(assayNames(object), "imputed"))
    assert_scalar_subset(subgroupvar, svars(object))
    assert_is_subset(transforms, TRANSFORMENGINES)
    assert_scalar_subset(method, DIMREDENGINES)
    assert_is_numeric(dims)
    assert_are_same_length(dims, numeric(2))
    assert_is_a_bool(verbose)
    
    drstrings <- .dimredstrings(object, method, dims)
    
    scoredt <- lapply(
      assays,
      function(as){
        tmpobj <- object
        .ldt_transforms_dimred(
          object, as, transforms, method, drstrings['methodname'], dims, verbose)
      }) %>%
      rbindlist()
    scoredt$assay %<>% factor(unique(.))
    scoredt$transform %<>% factor(unique(.))
    scoredt$variance.annot.ext %<>% factor(unique(.))
    scoredt$variance.annot %<>% factor(unique(.))
    
    .sdata_biplot(
      scoredt,
      x        = drstrings['x'],
      y        = drstrings['y'],
      method   = method, 
      dims     = dims,
      color    = color,
      shape    = shape,
      size     = size,
      alpha    = alpha,
      group    = group,
      label    = label,
      ... ) +
      geom_text(
        data = scoredt[, .SD[1], by = .(assay, transform)],
        aes(x = -Inf, y = Inf, label = variance.annot),
        color = "black",vjust = "inward", hjust = "inward") +
      facet_grid(rows = vars(assay), cols = vars(transform), scales = "free")
}

.sdata_biplot <- function(
    sdata,
    assay,
    x        = scorenames(
                 method, by = by, dims = dims[[1]], sep = guess_fitsep(sdata)),
    y        = scorenames(
                 method, by = by, dims = dims[[2]], sep = guess_fitsep(sdata)),
    method   = DIMREDENGINES[1], # 'pca'
    dims     = 1:2,
    color    = if (method %in% DIMREDSUPER) by else 'subgroup', 
    shape    = NULL, 
    size     = NULL, 
    alpha    = NULL,
    group    = NULL,         # Use 'feature_id' (not 'gene')
    label    = NULL,         # Which use 'feature_id' to guarantee uniqueness
    fixed    = list(shape = 15, size = 3), 
    colorpalette = make_colors(sort(unique(sdata[[color]]))),
    alphapalette = if (is.null(alpha)) NULL
                     else .make_alpha_palette(sdata[[alpha]]), 
    title = paste0(method, guess_fitsep(sdata, by)), 
    theme = ggplot2::theme(plot.title = element_text(hjust = 0.5), 
                           panel.grid = element_blank()),
    ...
){
    # Assert / Process
    assert_is_all_of(sdata, 'data.table')
    assert_scalar_subset(x, names(sdata))
    xsym <- symbolize(x)
    assert_scalar_subset(y, names(sdata))
    ysym <- symbolize(y)
    if (!is.null(color)) assert_scalar_subset(color, names(sdata))
    colorsym <- symbolize(color)
    if (!is.null(shape)){
      assert_scalar_subset(shape, names(sdata))
      fixed %<>% extract(names(.) %>% setdiff('shape'))
    }
    shapesym <- symbolize(shape)
    if (!is.null(size)){
      assert_scalar_subset(size, names(sdata))
      fixed %<>% extract(names(.) %>% setdiff('size'))
    }
    sizesym <- symbolize(size)
    if (!is.null(alpha)) assert_scalar_subset(alpha, names(sdata))
    alphasym <- symbolize(alpha)
    if (!is.null(group)) assert_scalar_subset(group, names(sdata))
    groupsym <- symbolize(group)
    if (!is.null(label)) assert_scalar_subset(label, names(sdata))
    labelsym <- symbolize(label)
 
    # Plot
    mapped_aes <- aes(
      x = !!xsym, y = !!ysym, color = !!colorsym, shape = !!shapesym,
      size = !!sizesym, alpha = !!alphasym, group = !!groupsym)
    p <- ggplot() + theme_bw() + theme +
      layer(
        geom     = 'point',
        stat     = 'identity',
        position = 'identity',
        data     = sdata,
        mapping  = mapped_aes,
        params = fixed,
        ...)
    if (!is.null(colorpalette))  p <- p +
      scale_color_manual(values = colorpalette, na.value = 'gray80')
    if (!is.null(alphapalette))  p <- p +
      scale_alpha_manual(values = alphapalette)
    if (!is.null(label))  p <- p + 
      geom_text_repel(
        data        = sdata,
        mapping     = utils::modifyList(mapped_aes, aes(label = !!labelsym)),
        show.legend = FALSE)
    if (!is.null(shape)){
      n <- if (is.factor(sdata[[shape]]))
        { length(levels(sdata[[shape]])) }
      else 
        { length(unique(sdata[[shape]])) }
      if (n > 6)  p <- p + scale_shape_manual(values = seq(15, 15+n-1))
      # Warning messages: The shape palette can deal with a maximum 
      # of 6 discrete values
      # https://stackoverflow.com/questions/16813278
    }
    # Return
    p
}

#' Data Transformation Methods
#' @details  \itemize{
#'     \item \code{TRANSFORMENGINES: c('center', 'center_mean', 'center_median', 'invnorm', 'quantnorm', 'vsn', 'zscore')}
#'     \item \code{TRANSFORMSTRICT:  c('center', 'invnorm', 'quantnorm', 'vsn', 'zscore')}
#' }
#' @export
TRANSFORMENGINES <- c('center', 'center_mean', 'center_median', 'invnorm',
                      'quantnorm', 'vsn', 'zscore')

#' @rdname TRANSFORMENGINES
#' @export
TRANSFORMSTRICT <- c('center', 'invnorm', 'quantnorm', 'vsn', 'zscore')
