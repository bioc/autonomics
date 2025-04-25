#========================================================================================
#
#    Mixture Modeling
#
#========================================================================================


#' Mixture k
#' 
#' Mixture modeling k 
#' 
#' Mixture modeling k (default number of components)
#' 
#' @param engine 
#' @return NULL or number
#' @examples
#' mixk('none')
#' mixk('mclust')
#' mixk('mixtools')
#' @export
mixk <- function(engine){
    assert_scalar_subset(engine, c('none', 'mclust', 'mixtools'))
    if (engine == 'none')      return(3)     # default for non-mixture binning
    if (engine == 'mclust')    return(NULL)  # default for mclust
    if (engine == 'mixtools')  return(2)     # default for mixtools
}


#' Mclust mixture modeling 
#' 
#' @param x  numeric vector
#' @param k  number of mixture components to look for  
#' @return data.table: component, mean, sd, weight
#' @examples
#' set.seed(1)
#' x <- c( rnorm(20,3) , rnorm(20,7) , rnorm(20,11))
#' mixmod_mclust(x)
#' @export
mixmod_mclust <- function(x, k = mixk('mclust')){
    if (!installed('mclust'))  return(NULL)
    mclustBIC <- mclust::mclustBIC
    fit <- mclust::Mclust(x, verbose = FALSE, G = k)
    means <- fit$parameters$mean
    sds <- fit$parameters$variance$sigmasq
    sds %<>% sqrt()    # can be scalar (common variance) or vector (different variances)
    if (length(sds)==1)  sds %<>% rep(length(means))
    weights <- fit$parameters$pro
    data.table( component = seq_along(means), mean = means,  sd = sds, weight = weights )
}


#' Mixtools mixture modeling 
#' 
#' @param x  numeric vector
#' @param k  number of mixture components to look for  
#' @return data.table: component, mean, sd, weight
#' @examples
#' set.seed(1)
#' x <- c( rnorm(20,3) , rnorm(20,7) , rnorm(20,11))
#' mixmod_mixtools(x)
#' mixmod_mixtools(x, k = 3)
#' @export
mixmod_mixtools <- function(x, k = mixk('mixtools')){
    if (!installed('mixtools'))  return(NULL)
        fit <- mixtools::normalmixEM(x, k = k)  # verbose parameter seems to be not working
      means <- fit$mu
        sds <- fit$sigma
    weights <- fit$lambda 
    data.table( component = seq_along(means), mean = means,  sd = sds, weight = weights )
}



#' Mixture modeling
#' @param x       numeric vector
#' @param engine 'mclust' or 'mixtools' 
#' @param k       number of components
#' @param plot    whether to plot
#' @return data.table: component, mean, sd, weight
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, 3), rnorm(20,7), rnorm(20, 11))
#' mixmod(x)
#' mixmod(x, engine = 'mixtools')
#' @export
mixmod <- function(x, engine = 'mclust', k = mixk(engine), plot = FALSE, color = '#F8766D'){
    assert_is_numeric(x)
    assert_scalar_subset(engine, c('mclust', 'mixtools'))
    mixdt <- if (engine == 'mclust'  ){  mixmod_mclust(  x, k = k)
      } else if (engine == 'mixtools'){  mixmod_mixtools(x, k = k) }
    if (plot)  print(mixplot(x, mean = mixdt$mean, 
                                  sd = mixdt$sd, 
                              weight = mixdt$weight, 
                              engine = engine, 
                               color = color))
    mixdt
}


wnorm <- function(x, mean, sd, weight)   weight*dnorm(x, mean = mean, sd = sd)


#' Mixture plot
#' @param x      data points
#' @param mean   component means
#' @param sd     component sds
#' @param weight component weights
#' @param engine 'none', 'mclust' or 'mixtools'
#' @param color string
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, 3), rnorm(20,7), rnorm(20, 11))
#' mixdt <- mixmod(x)
#' mixplot(x, mixdt$mean, mixdt$sd, mixdt$weight)
#' @export
mixplot <- function(x, mean, sd, weight, engine = '', color = '#F8766D'){
    xcurve <- seq(min(x), max(x), length.out = 100)
    ycurve <- mapply(wnorm, mean = mean, sd = sd, weight = weight, MoreArgs = list(x = xcurve), SIMPLIFY = FALSE)
    ycurve %<>% Reduce(cbind, .)
    ycurve %<>% rowSums()
    
    pointdt <- data.table(x = x,     y = .densities(x))
    curvedt <- data.table(x = xcurve, y = .densities(x, xcurve))
    mixdt   <- data.table(x = xcurve, y = ycurve, engine = engine)
    
    p <- ggplot() + theme_bw() + theme(panel.grid = element_blank())
    p <- p + geom_point(aes(x = x, y = y), pointdt, color = color)
    p <- p + geom_line( aes(x = x, y = y), curvedt, color = color)
    p <- p + geom_line( aes(x = x, y = y, linetype = engine), mixdt, color = color)
    p <- p + scale_linetype_manual(values = 'dotted')
    
    idx <- 1+which(diff(sign(diff(ycurve))) > 0)
    segmentdt <- data.table( x = xcurve[idx], 
                          xend = xcurve[idx],
                             y = min(ycurve), 
                          yend = .densities(x, xcurve[idx]))
    p <- p + geom_segment(aes(x = x, xend = xend, y = y, yend = yend), segmentdt, color = color)
    p <- p + geom_label(  aes(x = x, y = y+(yend-y)/2, label = formatC(x, 2)), segmentdt, color = color)
    p
}


#' Quadratic roots
#' 
#' Solves ax^2+bx+c = 0
#' 
#' Computes roots of quadratic equation
#' 
#'     D = b^2-4ac
#'     
#'           -b +- sqrt(D)
#'     x =  -------------
#'              2a
#' 
#' @param a coefficient of x^2
#' @param b coefficient of x^1
#' @param c coefficient of x^0
#' @examples
#' quadrroots(a = 1, b =-5, c = 6)  # two real roots
#' quadrroots(a = 1, b =-4, c = 4)  # one real root
#' quadrroots(a = 1, b = 1, c = 1)  # imaginary root
#' @return vector
#' @export
quadroots <- function(a,b,c){
    D <- b^2 - 4*a*c
    if (a == 0) return(-c/b)
    if (D > 0 ) return( c( (-b + sqrt(D)) / (2 * a),
                           (-b - sqrt(D)) / (2 * a) ) )
    if (D == 0) return( -b / (2*a) )
    # No real roots (two complex conjugate roots)
    real_part <- -b / (2 * a)
    imaginary_part <- sqrt(abs(D)) / (2 * a)
    root1 <- complex(real = real_part, imaginary =  imaginary_part)
    root2 <- complex(real = real_part, imaginary = -imaginary_part)
    return(c(root1 = root1, root2 = root2))
}


.mixbreaks <- function(mean1, mean2, sd1, sd2){
    # https://stats.stackexchange.com/a/311596
    var1 <- sd1^2
    var2 <- sd2^2
    a <- -1/var1 + 1/var2
    b <- 2*(-mean2/var2 + mean1/var1)
    c <- mean2^2/var2 - mean1^2/var1 + log(var2/var1)
    quadroots(a, b, c)
}


#' Mixture modeling breaks
#' @param x numeric vector
#' @return vector 
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, 3), rnorm(20,7), rnorm(20, 11))
#' mixbreaks(x)
#' @export
mixbreaks <- function(x, engine = 'mclust', k = mixk(engine)){
    mixdt <- mixmod(x, engine = engine, k = k)
    y <- lapply(  seq(1, nrow(mixdt)-1), 
                  function(i)  mixdt[ , .mixbreaks(mean[i], mean[i+1], sd[i], sd[i+1] ) ]  )
    y %<>% Reduce(c, .)
    y
} 


#========================================================================================
#
#    Factorize/Bin
#
#========================================================================================


#' Factorize/Bin
#' @details 
#'          `bin` transform into numeric bins : c(1,2,3,4,5,6) -> c( 1,  1,  2,  2,  3,  3 )
#'    `factorize` transform into factor levels: c(1,2,3,4,5,6) -> c('1','1','2','2','3','3')
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
factorize.numeric <- function(x, mixmod = 'none', k = mixk(mixmod), numericlevels = TRUE, ...){
    assert_is_a_number(k)
    assert_scalar_subset(mixmod, c('none', 'mclust', 'mixtools'))
     probs <- if (mixmod == 'none')  seq_len(k-1)/k  else NULL
    breaks <- if (mixmod == 'none'){  unname(quantile(x, probs = probs, na.rm = TRUE))
              } else {          autonomics::mixbreaks(x, engine = {{mixmod}}, k = k)   }
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
bin.numeric <- function(x, mixmod = 'none', k = mixk(mixmod), numericlevels = TRUE, ...){
    y <- factorize.numeric(x, mixmod = mixmod, k = k, numericlevels = TRUE)
    y %<>% as.numeric()
    y
}


#' @rdname factorize
#' @export
factorize.matrix <- function(x, mixmod = 'none', k = mixk(mixmod), numericlevels = TRUE, ...){
    y <- x
    y %<>% apply(1, factorize.numeric, k = k, numericlevels = numericlevels) %>% t()
    colnames(y) <- colnames(x)
    y
}


#' @rdname factorize
#' @export
bin.matrix <- function(x, mixmod = 'none', k = mixk(mixmod), numericlevels = TRUE, ...){
    y <- x
    y %<>% apply(1, bin.numeric, k = k, numericlevels = numericlevels) %>% t()
  # y %>% apply(1, dplyr::ntile, n = k) %>% t()    # differs a bit
    colnames(y) <- colnames(x)
    y
}


#' @rdname factorize
#' @export
factorize.SummarizedExperiment <- function(x, assay = assayNames(object)[1], 
    mixmod = 'none', k = mixk(mixmod), numericlevels = TRUE, verbose = TRUE, ...
){
    # Assert
    assert_scalar_subset(assay, assayNames(x))
    assert_is_a_bool(verbose)
    
    # Bin
    mat <- assays(x)[[assay]]
    mat %<>% factorize.matrix(mixmod = mixmod, k = k, numericlevels = numericlevels)
    
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


