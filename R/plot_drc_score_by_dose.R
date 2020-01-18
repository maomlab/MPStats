
#' Plot the score as a function of compound dose
#' 
#' For each compound plot a score by log_dose facet panel with
#'  *  the scores for each dose shown as the pooled score across replicates
#'  *  error bars as the standard deviation of the per-dose Bayesian uncertainty with a flat prior
#'  *  drc 4-parameter sigmoid with fixed to 0.
#'  *  labeled with the p-value of a constant-value fit
#' 
#' @param well_scores tibble::tibble as output by read_well_scores
#' @param fits tibble::tibble as output by fit_drc_score_by_dose
#' @param subtitle string subtitle for plot
#' @return ggplot2 plot
#' 
#' @export
plot_drc_score_by_dose <- function(well_scores, fits, subtitle=NULL){
  
  compound_dose_scores <- well_scores %>%
    dplyr::filter(!is_control) %>%
    dplyr::group_by(log_dose, compound) %>%
    dplyr::summarize(
      n_positive = sum(n_positive),
      cell_count = sum(cell_count),
      prob_positive = n_positive/cell_count,
      prob_positive_low = binomial_quantile(n_positive, cell_count, .025),
      prob_positive_high = binomial_quantile(n_positive, cell_count, .975)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
      data=fits,
      mapping=ggplot2::aes(
        x=log_dose,
        y=pred_value),
      color="blue",
      size=1.5) +
    ggplot2::geom_errorbar(
      data=compound_dose_scores,
      mapping=ggplot2::aes(
        x=log_dose,
        ymin=prob_positive_low,
        ymax=prob_positive_high)) +
    ggplot2::geom_point(
      data=compound_dose_scores,
      mapping=ggplot2::aes(
        x=log_dose,
        y=prob_positive)) +
    geom_indicator(
      data=fits %>% dplyr::distinct(compound, p_value),
      mapping=ggplot2::aes(
        indicator=paste0("fit: ", p_value)),
      xpos="left",
      ypos="top",
      group=1) +
    ggplot2::ggtitle(
      label="Score by log dose",
      subtitle=subtitle) +
    ggplot2::scale_x_continuous(
      "log[Compound dose] (uM)") +
    ggplot2::scale_y_continuous(
      "Score",
      limits=c(0,1),
      labels=scales::percent_format()) +
    ggplot2::facet_wrap(~compound, scales="free_x")
}
