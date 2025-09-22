#=============================================================================
#
#               subgroup_matrix
#                   split_subgroup_levels
#                       split_subgroup_values
#                           split_values
#
#=============================================================================

split_values <- function(x){
    sep <- guess_sep(x)
    dt <- data.table::data.table(x = x)
    dt[, data.table::tstrsplit(x, sep) ]
}


split_subgroup_values <- function(object, subgroupvar){
    subgroupvalues <- svalues(object, subgroupvar)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}


split_subgroup_levels <- function(object, subgroupvar){
    subgrouplevels <- slevels(object, subgroupvar)
    cbind(subgroup = subgrouplevels, split_values(subgrouplevels))
}


#' @rdname subgroup_matrix
#' @export
subgroup_array <- function(object, subgroupvar){
    . <- NULL
    x <- slevels(object, subgroupvar)
    sep <- guess_sep(object)
    #x %<>% sort()
    dt <- data.table(subgroup = factor(x, x))
    components <- dt[, tstrsplit(subgroup, sep, fixed=TRUE)]
    for (i in seq_len(ncol(components)))   components[[i]] %<>%
                                    factor(., levels=unique(.))
    dt %<>% cbind(components)
    data.table::setorderv(dt, rev(names(components)))
    levels  <- dt[, -1] %>% lapply(unique)
    #levels[1:2] %<>% rev()
    nlevels <- levels %>% vapply(length, integer(1))
    array(dt$subgroup, dim = nlevels, dimnames = levels)
}


#' Get subgroup matrix
#'
#' Arrange (subgroup)levels in matrix
#'
#' @param object SummarizedExperiment
#' @param subgroupvar subgroup svar
#' @return matrix
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object$subgroup <- paste0(object$Diabetes, '.', object$subgroup)
#' subgroup_matrix(object, 'subgroup')
#' @export
subgroup_matrix <- function(object, subgroupvar){
    . <- NULL
    subgroup_array <- subgroup_array(object, subgroupvar)
    if (length(dim(subgroup_array)) == 1){
        return(matrix(subgroup_array, 
                      byrow = TRUE, 
                      nrow = 1, 
                      dimnames = list(NULL, subgroup_array)))  
    }
    otherdims <- names(dim(subgroup_array)) %>% setdiff('V1')
    ncol1   <- Reduce('*', dim(subgroup_array)[otherdims])
    colnames1 <- dimnames(subgroup_array)[otherdims] %>%
                 expand.grid()                       %>%
                 apply(1, paste0, collapse = '.')
    subgroupmat <- matrix(subgroup_array,
                        nrow = nrow(subgroup_array), ncol = ncol1,
                        dimnames=list(rownames(subgroup_array), colnames1))
    subgroupmat %>% extract(nrow(.):1, )
    #dt <- split_subgroup_levels(object)
    #subgroupmat <- as.matrix(data.table::dcast(
    #    dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    #subgroupmat %>% extract(rev(order(rownames(.))), order(colnames(.)))
}


#=============================================================================
#
#           contrast_subgroup_cols
#           contrast_subgroup_rows
#
#==============================================================================


#' Row/Col contrasts
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @return  matrix
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object$subgroup <- paste0(object$Diabetes, '.', object$Time)
#' subgroup_matrix(object, subgroupvar = 'subgroup')
#' contrast_subgroup_cols(object, subgroupvar = 'subgroup')
#' contrast_subgroup_rows(object, subgroupvar = 'subgroup')
#' @export
contrast_subgroup_cols <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s-%s',    # no space: as lm(contr.sdiff)
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s-%s',   # no space: as lm(contr.sdiff)
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s-%s',  # no space: as lm(contr.sdiff)
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s-%s', # no space: as lm(contr.sdiff)
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}


#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#
#=============================================================================


 true_names <- function(x) names(x)[ x]
false_names <- function(x) names(x)[!x]

fit_limma_contrastogram <- function(object, subgroupvar, design){
    colcontrasts <- contrast_subgroup_cols(object, subgroupvar)
    rowcontrasts <- contrast_subgroup_rows(object, subgroupvar)
    contrasts <-  c( c(t(colcontrasts)), c(t(rowcontrasts)))
    object %<>% linmod_limma(design = design, contrasts = contrasts)
    object
}

compute_connections <- function(
    object, subgroupvar, design,
    colors = make_colors(slevels(object, subgroupvar), guess_sep(object))
){
# subgroup matrix, difference contrasts, limma
    fdt(object) %<>% add_adjusted_pvalues('fdr')
    fdrvalues <- fdrmat(object)
    effects <- effectmat(object)
    colnames(fdrvalues) %<>% split_extract_fixed('~', 1)
    colnames(effects)   %<>% split_extract_fixed('~', 1)
    nsignif <- apply(fdrvalues < 0.05, 2, sum, na.rm=TRUE)
                #colSums( fdrvalues < 0.05, na.rm=TRUE)  # BREAKS ON SINGLE CONTR!
    nup     <- apply(fdrvalues < 0.05 & effects > 0, 2, sum, na.rm = TRUE)
    ndown   <- apply(fdrvalues < 0.05 & effects < 0, 2, sum, na.rm = TRUE)
# Create diagram
    sep <- guess_sep(object)
    subgroupmatrix <- subgroup_matrix(object, subgroupvar = subgroupvar)
    subgrouplevels <- c(t(subgroupmatrix))
    arrowsizes <- arrowcolors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    arrowlabels <- matrix("0", nrow = nrow(arrowsizes), ncol = ncol(arrowsizes),
                        dimnames = dimnames(arrowsizes))
# Add contrast numbers
    contrastmat  <- makeContrasts(contrasts = coefs(object), levels = design)
    for (contrastname in colnames(contrastmat)){
        contrastvector <- contrastmat[, contrastname]
        to   <- true_names(contrastvector>0)
        from <- if (any(contrastvector<0)) true_names(contrastvector<0) else to
        ns <- nsignif[[contrastname]]
        nu <- nup[[contrastname]]
        nd <- ndown[[contrastname]]
        arrowsizes[ to, from] <- nu#ns
        arrowsizes[ from, to] <- nd#ns
        arrowcolors[to, from] <- colors[[to]]
        arrowcolors[from, to] <- colors[[from]]
        arrowlabels[to, from] <- if (nu>0) nu else 0
                            #paste0(nu,  " %up% phantom(.)") else "phantom(.)"
        arrowlabels[from, to] <- if (nd>0) nd else 0
                            #paste0(nd," %down% phantom(.)") else "phantom(.)"
    }
# Return
    #arrowlabels[arrowcolors==0] <- "0"
    list(arrowsizes = arrowsizes,
        arrowcolors = arrowcolors,
        arrowlabels = arrowlabels)
}


#' Plot contrastogram
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @param formula      formula
#' @param colors       named color vector (names = subgroups)
#' @param curve        arrow curvature
#' @return list returned by \code{\link[diagram]{plotmat}}
#' @examples
#' if (installed('diagram')){
#'    file <- download_data('halama18.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    plot_contrastogram(object, subgroupvar = 'subgroup')
#' }
#' @export
plot_contrastogram <- function(
    object, 
    subgroupvar,
    formula = as.formula(paste0('~ 0 +', subgroupvar)),
    colors  = make_colors(slevels(object, subgroupvar), guess_sep(object)),
    curve   = 0.1
){
# Initialize
    V2 <- N <- NULL
    if (!installed('diagram'))  return(NULL)
# Fit limma
    design <- create_design(object, formula = formula)
    object %<>% fit_limma_contrastogram(subgroupvar = subgroupvar, design = design)
# Compute connections
    contrastogram_matrices <- compute_connections(object, design = design, subgroupvar = subgroupvar, colors = colors)
    arrowsizes  <- contrastogram_matrices$arrowsizes
    arrowcolors <- contrastogram_matrices$arrowcolors
    arrowlabels <- contrastogram_matrices$arrowlabels
    widths <- scales::rescale(arrowsizes, c(0.01,30))
# Plot
    dt <- split_subgroup_levels(object, subgroupvar)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    if (all(nperrow==1)) nperrow %<>% length()
    # basedir <- file.path(tools::R_user_dir('autonomics', 'cache'), 'contrastogram')
    # dir.create(basedir)
    # pdf(file.path(basedir, 'directed_contrastogram.pdf'), 
    # width = 9, height = 9)
    arrowlabels %<>% as.data.frame()
    diagram::plotmat(A          = arrowlabels,
                    pos         = nperrow,
                    curve       = curve,
                    name        = rownames(arrowsizes),
                    relsize     = 1,
                    box.size    = 0.05,
                    box.col     = colors[rownames(arrowsizes)],
                    box.type    = 'square',
                    box.prop    = 0.8,
                    arr.lwd     = widths,
                    shadow.size = 0, # sqrt(arrowsizes)
                    arr.lcol    = arrowcolors,
                    arr.col     = arrowcolors,
                    arr.type    = 'triangle')
    #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
# Return
    object # limma!
}


