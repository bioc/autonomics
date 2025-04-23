
#' Factorize/Bin
#' @details 
#'          `bin` transforms into numeric bins : c(1,2,3,4,5,6) -> c( 1,  1,  2,  2,  3,  3 )
#'    `factorize` transforms into factor levels: c(1,2,3,4,5,6) -> c('1','1','2','2','3','3')
#' @param x       vector, matrix or SummarizedExperiment
#' @param probs   numeric
#' @param k       number of bins/levels
#' @param assay   string
#' @param verbose TRUE or FALSE
#' @param ... (S3 dispatch)
#' @return  vector, matrix or SummarizedExperiment
#' @examples 
#' # data 
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     object <- read_maxquant_proteingroups(file, impute = TRUE)
#'     fdt(object)
#' 
#' # logical
#'     fdt(object)$imputed
#'     fdt(object)$imputed %>% factorize()
#'     fdt(object)$imputed %>% bin()
#'     
#' # character
#'     as.character(fdt(object)$imputed)
#'     as.character(fdt(object)$imputed) %>% factorize()
#'     as.character(fdt(object)$imputed) %>% bin()
#' 
#' # factor
#'     factor(fdt(object)$imputed)
#'     factor(fdt(object)$imputed) %>% factorize()
#'     factor(fdt(object)$imputed) %>% bin()
#'     
#' # numeric
#'     fdt(object)$pepcounts
#'     fdt(object)$pepcounts %>% factorize()
#'     fdt(object)$pepcounts %>% bin()
#' 
#' # Matrix/SummarizedExperiment
#'     values(object)
#'     values(object) %>% factorize()
#'            object  %>% factorize()
#'     values(object) %>% bin()
#'            object  %>% bin()
#' @export
factorize <- function(x, ...)  UseMethod('factorize')


#' @rdname factorize
#' @export
bin <- function(x, ...)  UseMethod('bin')


#' @rdname factorize
#' @export
factorize.logical <- function(x, ...) as.factor(x)


#' @rdname factorize
#' @export
bin.logical <- function(x, ...)    as.numeric(x)


#' @rdname factorize
#' @export
factorize.character <- function(x, ...)  as.factor(x)


#' @rdname factorize
#' @export
bin.character <- function(x, ...)  as.numeric(as.factor(x))


#' @rdname factorize
#' @export
factorize.factor <- function(x, ...)  x


#' @rdname factorize
#' @export
bin.factor <- function(x, ...)    as.numeric(x)


#' @rdname factorize
#' @export
factorize.numeric <- function(x, k = 3, probs = seq_len(k-1)/k, numericlevels = TRUE, ...){
    assert_is_a_number(k)
    breaks <- quantile(x, probs = probs, na.rm = TRUE)
    breaks %<>% c(  `0%` = min(x, na.rm = TRUE)-1e-7, .)
    breaks %<>% c(`100%` = max(x, na.rm = TRUE)+1e-7   )
    y <- cut(x, breaks)
    if (numericlevels){  levels(y) %<>% seq_along()
    } else {             levels(y) %<>% substr(2, nchar(.))
                         levels(y) %<>% split_extract_fixed(',', 1)
                         levels(y) %<>% paste0('>', .)
    }
    y
}


#' @rdname factorize
#' @export
bin.numeric <- function(x, k = 3, probs = seq_len(k-1)/k, ...){
    as.numeric(factorize.numeric(x, k = k, probs = probs, numericlevels = TRUE))
}


#' @rdname factorize
#' @export
factorize.matrix <- function(x, k = 3, probs = seq_len(k-1)/k, numericlevels = TRUE, ...){
    y <- x
    y %<>% apply(1, factorize.numeric, k = k, probs = probs, numericlevels = numericlevels) %>% t()
    colnames(y) <- colnames(x)
    y
}


#' @rdname factorize
#' @export
bin.matrix <- function(x, k = 3, probs = seq_len(k-1)/k, ...){
    y <- x
    y %<>% apply(1, bin.numeric, k = k, probs = probs, numericlevels = numericlevels) %>% t()
  # y %>% apply(1, dplyr::ntile, n = k) %>% t()    # differs a bit
    colnames(y) <- colnames(x)
    y
}


#' @rdname factorize
#' @export
factorize.SummarizedExperiment <- function(
    x, assay = assayNames(object)[1], k = 3, probs = seq_len(k-1)/k, numericlevels = TRUE, verbose = TRUE, ...
){
    # Assert
    assert_scalar_subset(assay, assayNames(x))
    assert_is_a_bool(verbose)
    
    # Bin
    mat <- assays(x)[[assay]]
    mat %<>% factorize.matrix(k = k, probs = probs, numericlevels = numericlevels)
    
    # Add
    newassayname <- sprintf('%s%dlevels', assay, k)
    if (verbose)   cmessage('%sAdd  `%s`', spaces(14), newassayname)  # Align with Code `exprs2levels``
    assays(x)[[newassayname]] <- mat
    x
}


#' @rdname factorize
#' @export
factorize_assay <- function(object, assay = assayNames(object)[1], k = 3, verbose = TRUE){
    .Deprecated('factorize') # factorize.SummarizedExperiment
    factorize.SummerizedExperiment(object, assay = assay, k = k, verbose = verbose)
}



#' @rdname factorize
#' @export
bin.SummarizedExperiment <- function(
    x, assay = assayNames(x)[1], k = 3, probs = seq_len(k-1)/k, verbose = TRUE
){
    # Assert
    assert_scalar_subset(assay, assayNames(x))
    assert_is_a_bool(verbose)
    
    # Bin
    mat <- assays(x)[[assay]]
    mat %<>% bin.matrix(k = k, probs = probs)
    
    # Add
    newassayname <- sprintf('%s%dbins', assay, k)
    if (verbose)   cmessage('%sAdd  `%s`', spaces(14), newassayname)  # Align with Code `exprs2levels``
    assays(x)[[newassayname]] <- mat
    x
}

#' @rdname factorize
#' @export
bin_assay <- function(object, assay = assayNames(object)[1], k = 3, verbose = TRUE){
    .Deprecated('bin') # bin.SummarizedExperiment
    bin.SummarizedExperiment(object, assay = assay, k = k, verbose = verbose)
}


