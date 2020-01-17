

#' A diagnostic plot for error analysis for binomial trials
#' 
#' @return ggplot2 plot of probability positive by cell count
#' @export
plot_binomial_variance_by_cell_count <- function() {
  
  data <- expand.grid(
    prob_positive = c(0, .1, .3, .6, 1),
    cell_count = seq(1, 200, length.out=200))
  
  p <- ggplot2::ggplot(data=data) + 
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
      mapping=ggplot2::aes(
        x=cell_count,
        y=prob_positive,
        ymin=binomial_quantile(prob_positive*cell_count, cell_count, .025),
        ymax=binomial_quantile(prob_positive*cell_count, cell_count, .975),
        group=prob_positive),
      stat="identity") +
    ggplot2::ggtitle(
      label="Bayesian 95% credible intervals for Binomial trials") +
    ggplot2::scale_y_continuous("Probability Positive") +
    ggplot2::scale_x_continuous("Cell Count")
  p
}