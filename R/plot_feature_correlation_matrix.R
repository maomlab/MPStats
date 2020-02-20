
#' plot feature heatmap
#'
#' @param x matrix to be plotted
#' @return gplots object with a matrix where cell (i,j) is colored by the in the cell
#'
#' @usage For example, to plot the correlation
#'
#' 
#'   heatmap_args <- cell_features %>%
#'     stats::cor() %>%
#'     MPStats::plot_heatmap()
#'   pdf(fname, height=height, width=width)
#'   suppressWarnings(ret <- do.call(gplots::heatmap.2, heatmap_args))
#    dev.off()
#'
#'@export
plot_heatmap <- function(
  x,
  ref_x=x,
  subtitle=NULL){
  if(!is.null(ref_x)){
    dist_row <- dist(ref_x)
    o_row <- seriation::seriate(dist_row, method = "OLO", control = NULL)[[1]]
    Rowv=as.dendrogram(o_row)
  } else {
    Rowv = FALSE
  }
  args <- list()
  if (any(x < 0, na.rm = TRUE)){
    args$col <- seriation::bluered(n=100, bias=1)
  } else {
    args$col <- seriation::greys(n=100, power=1)
  }
  args$trace <- "none"
  args$density.info <- "none"
  args$cexRow <- 1
  args$cexCol <- 1
  args$dendrogram="none"
  args$key = FALSE
  args$keysize = 0.03
  args$colsep = seq(0, ncol(x), by=5000000)
  args$rowsep = seq(0, nrow(x), by=5000000)
  #args$sepwidth = c(0.02, 0.02)
  args$sepwidth = c(0, 0)
  args$margins = c(3,7)
  args <- c(list(x = x, Colv = FALSE, Rowv = Rowv), args)
  return(args)
}