#========================================================================================
#
#    Mixture Modeling
#
#========================================================================================


mixmod_mclust <- function(x, k = NULL){
    if (!installed('mclust'))    return( mixmod_none(x) )
    mclustBIC <- mclust::mclustBIC
    fit <- mclust::Mclust(x, verbose = FALSE, G = k)
    means <- fit$parameters$mean
    sds <- fit$parameters$variance$sigmasq
    sds %<>% sqrt()    # can be scalar (common variance) or vector (different variances)
    if (length(sds)==1)  sds %<>% rep(length(means))
    weights <- fit$parameters$pro
    data.table( component = seq_along(means), mean = means,  sd = sds, weight = weights )
}


mixmod_mixtools <- function(x, k = 2){
    if (!installed('mixtools'))  return( mixmod('none') )
        fit <- mixtools::normalmixEM(x, k = k)  # verbose parameter seems to be not working
      means <- fit$mu
        sds <- fit$sigma
    weights <- fit$lambda 
    data.table( component = seq_along(means), mean = means,  sd = sds, weight = weights )
}



mixmod_none <- function(x)   return( data.table(  component = 1, 
                                                       mean = mean(x, na.rm = TRUE),
                                                         sd = sd(  x, na.rm = TRUE), 
                                                     weight = 1  ) )

#' Mixture model
#' @param x        numeric vector
#' @param engine  'mclust', 'mixtools', 'none'
#' @param k        number of components
#' @param color    string
#' @return  data.table (mixmod), ggplot (mixplot), vector (mixbreaks)
#' @examples
#' set.seed(1)
#' x <- c(rnorm(20, 3), rnorm(20,7), rnorm(20, 11))
#' mixmod(x)
#' mixplot(x)
#' mixbreaks(x)
#' @export
mixmod <- function(
         x, 
    engine = 'mclust', 
         k = switch(engine, none = 3, mclust = NULL, mixtools = 3) 
){
    assert_scalar_subset(engine, c('none', 'mclust', 'mixtools'))
    switch(engine, mclust = mixmod_mclust(x, k = k), 
                 mixtools = mixmod_mixtools(x, k = k), 
                   single = mixmod_none(x))
}


wnorm <- function(x, mean, sd, weight)   weight*dnorm(x, mean = mean, sd = sd)


#' @rdname mixmod
#' @export
mixplot <- function(
         x, 
    engine = 'mclust', 
         k = switch(engine, none = 3, mclust = NULL, mixtools = 3),
     color = '#F8766D'
){

# Model
    assert_scalar_subset(engine, c('none', 'mclust', 'mixtools'))
    y <- xend <- yend <- NULL
     mixdt <- mixmod(x, engine = engine, k = k)
      mean <- mixdt$mean
        sd <- mixdt$sd
    weight <- mixdt$weight
# Prep
    xcurve <- seq(min(x), max(x), length.out = 100)
    ycurve <- mapply(wnorm, mean = mean, sd = sd, weight = weight, MoreArgs = list(x = xcurve), SIMPLIFY = FALSE)
    ycurve %<>% Reduce(`+`, .)

    pointdt <- data.table(x = x,     y = .densities(x))
    curvedt <- data.table(x = xcurve, y = .densities(x, xcurve))
    mixdt   <- data.table(x = xcurve, y = ycurve, engine = engine)
# Plot
    p <- ggplot() + theme_bw() + theme(panel.grid = element_blank())
    p <- p + geom_point(aes(x = x, y = y), pointdt, color = color)
    p <- p + geom_line( aes(x = x, y = y), curvedt, color = color)
    p <- p + geom_line( aes(x = x, y = y, linetype = engine), mixdt, color = color)
    p <- p + scale_linetype_manual(values = 'dotted')
    mbreaks <- mixbreaks(x)    
    if (length(mbreaks) == 0)   return(p)
    segmentdt <- data.table( x = mbreaks, 
                          xend = mbreaks,
                             y = min(ycurve), 
                          yend = .densities(x, mbreaks))
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
#' @noRd
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


#' @rdname mixmod
#' @export
mixbreaks <- function(
         x, 
    engine = 'mclust', 
         k = switch(engine, none = 3, mclust = NULL, mixtools = 3) 
){
    mixdt <- mixmod(x, engine = engine, k = k)
    if (nrow(mixdt) == 1)  return(c())
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
#' @param mixmod  'none', 'mclust', or 'mixtools' (mixture modeling method to use)
#' @param numericlevels TRUE (levels: 1,2, ...) or FALSE (levels: 2.1+, 3.2+, ...)
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
factorize.logical <- function(x, ...) as.factor(x)


#' @rdname factorize
#' @export
factorize.character <- function(x, ...)  as.factor(x)


#' @rdname factorize
#' @export
factorize.factor <- function(x, ...)  x


#' @rdname factorize
#' @export
factorize.numeric <- function(
                x, 
           mixmod = 'none', 
                k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
    numericlevels = TRUE, ...
){
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
factorize.matrix <- function(
                x, 
           mixmod = 'none', 
                k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
    numericlevels = TRUE, 
                 ...
){
    y <- x
    y %<>% apply(1, factorize.numeric, k = k, numericlevels = numericlevels) %>% t()
    colnames(y) <- colnames(x)
    y
}



#' @rdname factorize
#' @export
factorize.SummarizedExperiment <- function(
                x, 
            assay = assayNames(x)[1],
           mixmod = 'none',
                k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
    numericlevels = TRUE,
             drop = TRUE,
          verbose = TRUE,
                 ...
){
    # Assert
    assert_scalar_subset(assay, assayNames(x))
    assert_is_a_bool(verbose)
    
    # Bin
    mat <- assays(x)[[assay]]
    mat %<>% factorize.matrix(mixmod = mixmod, k = k, numericlevels = numericlevels)
    if (!drop)  mat[] %<>% paste0(assay, .)

    # Add
    newassayname <- sprintf('%s%dlevels', assay, k)
    if (verbose)   cmessage('%sAdd  `%s`', spaces(14), newassayname)  # Align with Code `exprs2levels``
    assays(x)[[newassayname]] <- mat
    x
}


#' @rdname factorize
#' @export
factorize_assay <- function(
         x, 
          assay = assayNames(x)[1], 
         mixmod = 'none',
              k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
        verbose = TRUE, 
               ...
){
    .Deprecated('factorize') # factorize.SummarizedExperiment
    factorize.SummarizedExperiment(x, assay = assay, k = k, verbose = verbose)
}



#' @rdname factorize
#' @export
bin <- function(x, ...)  UseMethod('bin')


#' @rdname factorize
#' @export
bin.logical <- function(x, ...)    as.numeric(x)


#' @rdname factorize
#' @export
bin.character <- function(x, ...)  as.numeric(as.factor(x))


#' @rdname factorize
#' @export
bin.factor <- function(x, ...)    as.numeric(x)



#' @rdname factorize
#' @export
bin.numeric <- function(
                x, 
           mixmod = 'none',
                k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
    numericlevels = TRUE, 
    ...
){
    y <- factorize.numeric(x, mixmod = mixmod, k = k, numericlevels = TRUE)
    y %<>% as.numeric()
    y
}


#' @rdname factorize
#' @export
bin.matrix <- function(
                x, 
           mixmod = 'none', 
                k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
    numericlevels = TRUE, 
                 ...
){
    y <- x
    y %<>% apply(1, bin.numeric, k = k, numericlevels = numericlevels) %>% t()
  # y %>% apply(1, dplyr::ntile, n = k) %>% t()    # differs a bit
    colnames(y) <- colnames(x)
    y
}



#' @rdname factorize
#' @export
bin.SummarizedExperiment <- function(
          x, 
      assay = assayNames(x)[1],
     mixmod = 'none',
          k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
      probs = seq_len(k-1)/k,
    verbose = TRUE, 
           ...
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
bin_assay <- function(
     x, 
      assay = assayNames(x)[1], 
     mixmod = 'none',
          k = switch(mixmod, none = 3, mclust = NULL, mixtools = 3),
    verbose = TRUE
){
    .Deprecated('bin') # bin.SummarizedExperiment
    bin.SummarizedExperiment(x, assay = assay, k = k, verbose = verbose)
}


get_density <- function(x, y, ...) {
    # Kamil Slowikowski
    # https://slowkow.com/notes/ggplot2-color-by-density/
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}


#' @rdname densities
#' @export
.densities <- function(x, xpred = x)  approxfun(density(x, na.rm = TRUE))(xpred)



#' Densities
#' @param x      numeric vector: data points
#' @param xpred  numeric vector: prediction points
#' @param plot   whether to plot
#' @param color  string
#' @return numeric vector with same length as xpred
#' @examples
#'  set.seed(1)
#'  x <- c(rnorm(20, 3), rnorm(20,7), rnorm(20, 11))
#'  xpred <- seq(min(x), max(x), length.out = 100)
#' .densities(x, xpred)  # innerfun
#'  densities(x, xpred)  # outerfun
#' @export
densities <- function(
        x, 
    xpred = x, 
     plot = TRUE, 
    color = "#F8766D"
){
# Assert
    assert_is_numeric(x)
    assert_is_numeric(xpred)
    assert_is_a_bool(plot)
    assert_all_are_colors(color)
# Run
        y <- .densities(x)
    ypred <- .densities(x, xpred)
# Plot/Return
    if (plot){
        pointdt <- data.table(x = x, y = y)
        linedt  <- data.table(x = xpred, y = ypred, method = '')
        p <- ggplot() + theme_bw() + theme(panel.grid = element_blank())
        p <- p + geom_point(aes(x = x, y = y), pointdt, color = color)
        p <- p + geom_line( aes(x = x, y = y), linedt,  color = color)
        print(p)
    }
    ypred
}







#' @rdname plot_xy_density
#' @export
plot_x_density <- function(
                   x,
                   y = NULL,
             xbreaks = mixbreaks(x),
               title = NULL,
               color = '#F8766D',
                xlab = NULL,
                ylab = NULL,
          transcolor = '00000000'
){
# Prep
    x %<>% sort()
    densityfun <- approxfun(density(x, na.rm = TRUE))
    xpath <- seq(min(x), max(x), length.out = 100)
# Plot
    p <- ggplot() + theme_bw()
    p <- p + annotate('point', x = x,     y = densityfun(x),     color = color)
    p <- p + annotate('path',  x = xpath, y = densityfun(xpath), color = color)
# Breaks
    if (length(xbreaks)>0)  p <- p + annotate('segment', x = xbreaks, 
                                                      xend = xbreaks, 
                                                         y = 0.95*min(densityfun(xpath)),
                                                      yend = densityfun(xbreaks), 
                                                     color = color, 
                                                  linetype = 'dashed')
# Finishing
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
    p <- p + theme(panel.grid   = element_blank(), 
                   panel.border = element_blank())
    p <- p + theme(plot.title = element_text(color = color, hjust = 0.5))
    p <- p + scale_x_continuous(position = 'top')
    p <- p + theme(axis.line.x  = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.title.x = element_text(color = color))
    p <- p + theme(axis.line.y  = element_line(color = transcolor), 
                   axis.ticks.y = element_line(color = transcolor), 
                   axis.text.y  = element_text(color = transcolor), 
                   axis.title.y = element_text(color = transcolor))
    p <- p + theme(plot.margin = unit(c(5.5,0,0,5.5), 'points'))
# Align with xyplot
    if (!is.null(y)){
          digits <- max(nchar(scales::extended_breaks()(range(y)))) - 1
        labelfun <- function(br) round(br, digits = digits)
        p <- p + scale_y_continuous(labels = labelfun, sec.axis = sec_axis(~., labels = labelfun))
    }
    p
}


#' @rdname plot_xy_density
#' @export
plot_y_density <- function(
                   y,
                   x = NULL,
             ybreaks = mixbreaks(y),
               title = NULL,
               color = '#F8766D',
                xlab = NULL,
                ylab = NULL,
          transcolor = '00000000'
){
# Prep
    y %<>% sort()
    densityfun <- approxfun(density(y, na.rm = TRUE))
    ypath <- seq(min(y), max(y), length.out = 100)
# Plot    
    p <- ggplot() + theme_bw()
    p <- p + annotate('point',   y = y,     x = densityfun(y),     color = color)
    p <- p + annotate('path',    y = ypath, x = densityfun(ypath), color = color)
# Breaks
    if (length(ybreaks)>0)   p <- p + annotate('segment', y = ybreaks, 
                                                       yend = ybreaks, 
                                                          x = 0.95*min(densityfun(ypath)), 
                                                       xend = densityfun(ybreaks), 
                                                      color = color, 
                                                   linetype = 'dashed')
# Finishing
    p <- p + scale_y_continuous(position = 'right')
    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
    p <- p + theme(panel.grid = element_blank(), panel.border = element_blank())
    p <- p + theme(plot.title = element_text(color = color, hjust = 0.5))
    p <- p + theme(plot.margin = unit(c(0,5.5,5.5,0), 'points'))
    p <- p + theme(axis.line.y.left   = element_blank(), 
                   axis.line.y.right  = element_blank(),
                   axis.ticks.y.left  = element_blank(),
                   axis.ticks.y.right = element_blank(),
                   axis.text.y.left   = element_blank(),
                   axis.text.y.right  = element_blank(),
                   axis.title.y.left  = element_text(color = color, angle = 0, vjust = 0.5),
                   axis.title.y.right = element_text(color = color, angle = 0, vjust = 0.5)
    )
    p <- p + theme(axis.line.x  = element_line(color = transcolor), 
                   axis.ticks.x = element_line(color = transcolor), 
                   axis.text.x  = element_text(color = transcolor), 
                   axis.title.x = element_text(color = transcolor))
# Align with xy plot
    if (!is.null(x)){
          digits <- max(nchar(scales::extended_breaks()(range(x)))) - 1
        labelfun <- function(br) round(br, digits = digits)
        p <- p + scale_x_continuous(labels = labelfun, sec.axis = sec_axis(~., labels = labelfun))
    }
    p
}


#' @rdname plot_xy_density
#' @export
plot_xy_scatter <- function(
          x,
          y,
    xbreaks = mixbreaks(x),
    ybreaks = mixbreaks(y),
      color = c('#F8766D', '#00BFC4'),
    contour = FALSE,
     smooth = FALSE,
       xlab = NULL,
       ylab = NULL
){
    p <- ggplot(data.table(x = x, y = y), aes(x = x, y = y)) + theme_bw()
    if (contour) p <- p + geom_density2d(color = 'gray80')
    if (smooth ) p <- p + geom_smooth(   color = 'gray80', se = FALSE, method = 'lm', formula = y ~ x)
    p <- p + theme(plot.margin = unit(c(0,0,5.5,5.5), 'points'))
    p <- p + geom_point()
    if (length(xbreaks)>0)  p <- p + geom_vline(aes(xintercept = xbreaks), color = color[[1]], linetype = 'dashed')
    if (length(ybreaks)>0)  p <- p + geom_hline(aes(yintercept = ybreaks), color = color[[2]], linetype = 'dashed')
    p <- p + theme(panel.grid = element_blank())
    p <- p + scale_x_continuous(sec.axis = sec_axis(~.))
    p <- p + scale_y_continuous(sec.axis = sec_axis(~.))
    p <- p + theme(axis.line.x  = element_line(color = color[[1]]),
                   axis.ticks.x = element_line(color = color[[1]]),
                   axis.text.x  = element_text(color = color[[1]]), 
                   axis.title.x = element_text(color = color[[1]]))
    p <- p + theme(axis.line.y  = element_line(color = color[[2]]),
                   axis.ticks.y = element_line(color = color[[2]]),
                   axis.text.y  = element_text(color = color[[2]]), 
                   axis.title.y = element_text(color = color[[2]]))
    p <- p + theme(panel.border = element_blank())
    p <- p + xlab(xlab) + ylab(ylab)
    #p <- p + theme(plot.margin = margin(c(0,0,0,0), 'points'))
    p
}


#' Plot xy densities
#' @param x           numeric vector
#' @param y           numeric vector
#' @param xbreaks      numeric vector
#' @param ybreaks      numeric vector
#' @param title       NULL or string
#' @param color       vector or string
#' @param contour     TRUE or FALSE: plot density contours ?
#' @param smooth      TRUE or FALSE: plot smooth line ?
#' @param xlab        NULL or string
#' @param ylab        NULL or string
#' @param transcolor  string
#' @return ggplot
#' @examples
#' # Bimodal
#'     set.seed(1)
#'     x <- c(rnorm(10, 3), rnorm(10,7))
#'     y <- c(rnorm(10, 3), rnorm(10,7))
#'     plot_xy_density(x,y)
#'     plot_xy_density(x,y, contour = TRUE)
#'     plot_xy_density(x,y,  smooth = TRUE)
#'     plot_xy_scatter(x,y)
#'     plot_x_density(x)
#'     plot_y_density(y)
#' # Unimodal
#'     set.seed(1)
#'     x <- c(rnorm(20, 3))
#'     y <- c(rnorm(20, 3))
#'     plot_xy_density(x,y)
#'     plot_xy_scatter(x,y)
#'     plot_x_density(x)
#'     plot_y_density(y)
#' @export
plot_xy_density <- function(
              x,
              y,
         xbreaks = mixbreaks(x),
         ybreaks = mixbreaks(y),
            xlab = 'X',
            ylab = 'Y',
           color = c('#F8766D', '#00BFC4' ),
         contour = FALSE,
          smooth = FALSE
){
    px  <- plot_x_density( x, y, color = color[[1]], xlab = xlab, ylab = ylab)
    py  <- plot_y_density( y, x, color = color[[2]], xlab = xlab, ylab = ylab)
    pxy <- plot_xy_scatter(x, y, color = color, contour = contour, smooth = smooth, xlab = xlab, ylab = ylab)
    layout <- matrix(c(1,1,4,
                       2,2,3,
                       2,2,3), nrow = 3, byrow = TRUE)
    grid.arrange(px, pxy,py, layout_matrix = layout)
}


