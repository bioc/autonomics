
#' Bin/Factorize
#' @param x       vector, matrix or SummarizedExperiment
#' @param probs   numeric
#' @param k       number of bins/levels
#' @param assay   string
#' @param verbose TRUE or FALSE
#' @param ... (S3 dispatch)
#' @return  vector, matrix or SummarizedExperiment
#' @examples 
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file, impute = TRUE)
#' fdt(object)
#' bin(fdt(object)$imputed)    # bin.logical: unchanged
#' bin(fdt(object)$pepcounts)  # bin.numeric: binned, specific
#' bin(values(object))         # bin.matrix:  binned, agnostic
#' bin(object)                 # bin.SummarizedExperiment: binned, agnostic
#' @export
bin <- function(x, ...)  UseMethod('bin')


#' @rdname bin
#' @export
bin.logical <- function(x, ...)   x


#' @rdname bin
#' @export
bin.character <- function(x, ...) x


#' @rdname bin
#' @export
bin.factor <- function(x, ...)    x


#' @rdname bin
#' @export
bin.numeric <- function(x, probs = c(0, 0.33, 0.66, 1), ...){
    breaks <- quantile(x, probs = probs)
    breaks[1] %<>% subtract(1e-7)           # avoid smallest number from falling outside of bin
    x %<>% cut(breaks)                 # explicit breaks avoid negative bin
    levels(x) %<>% substr(2, nchar(.)) #    https://stackoverflow.com/questions/47189232
    levels(x) %<>% split_extract_fixed(',', 1)
    levels(x) %<>% paste0('>', .)
    x
}

#' @rdname bin
#' @export
bin.matrix <- function(x, k = 3, verbose = TRUE, ...){
    y <- x
    y %<>% apply(1, dplyr::ntile, n = k) %>% t()
    colnames(y) <- colnames(x)
    y
}


#' @rdname bin
#' @export
bin.SummarizedExperiment <- function(x, assay = assayNames(x)[1], k = 3, verbose = TRUE){
# Assert
    assert_scalar_subset(assay, assayNames(x))
    assert_is_a_number(k)
    assert_is_a_bool(verbose)
# Bin
    mat <- assays(x)[[assay]]
    mat %<>% bin.matrix(k = k)
# Add
    newassayname <- sprintf('%s%dbins', assay, k)
    if (verbose)   cmessage('%sAdd  `%s`', spaces(14), newassayname)  # Align with Code `exprs2levels``
    assays(x)[[newassayname]] <- mat
    x
}


#' @export
factorize <- function(x, ...)  UseMethod('factorize')

#' @export
factorize.SummerizedExperiment <- function(x, assay = assayNames(object)[1], k = 3, verbose = TRUE){
# Bin (assertions done during binning)
    object %<>% bin.SummarizedExperiment(assay = assay, k = k, verbose = verbose)
# Factorize
    binnedassay <- sprintf('%s%dbins', assay, k)
    mat <- assays(object)[[binnedassay]]
    mode(mat) <- 'character'
    mat %<>% paste0('bin', .)
    dim(mat) <- dim(object)
    dimnames(mat) <- dimnames(object)
# Add
    newassayname <- sprintf('%s%dlevels', assay, k)
    if (verbose)   cmessage('%sAdd  `%s`', spaces(14), newassayname)  # Align with Code `exprs2levels`
    assays(object)[[newassayname]] <- mat
    object
}
 

