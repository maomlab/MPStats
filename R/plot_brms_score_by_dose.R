
#' Plot the score as a function of compound dose for brms models
#' 
#' For each compound plot a score by log_dose facet panel with
#'  *  the scores for each dose shown as the pooled score across replicates
#'  *  error bars as the standard deviation of the per-dose Bayesian uncertainty with a flat prior
#'  *  plot the following models
#'  *     flat model
#'  #     4-parameter sigmoid with bottom fixed to 0.
#'  *  labeled with the p-value of a constant-value fit
#' 
#' @param well_scores tibble::tibble as output by read_well_scores
#' @param fits tibble::tibble as output by fit_drc_score_by_dose
#' @param subtitle string subtitle for plot
#' @return ggplot2 plot
#' 
#' @export
plot_brms_score_by_dose <- function(
  well_scores,
  flat_model,
  hill_model,
  subtitle=NULL){
  
  dose_ranges <- well_scores %>%
    dplyr::filter(!is_control) %>%
    dplyr::group_by(compound) %>%
    dplyr::summarize(
      min_log_dose = min(log_dose),
      max_log_dose = max(log_dose))
  
  flat_model_draws <- flat_model %>%
    dplyr::left_join(dose_ranges, by="compound") %>%
    purrr::pmap_dfr(function(compound, model, min_log_dose, max_log_dose, ...){
      cat("Fittting summary for compound '", compound, "'\n", sep="")
      model %>%
        brms:::posterior_samples(pars=c("Intercept")) %>%
        dplyr::sample_n(50) %>%
        dplyr::mutate(draw_id = dplyr::row_number()) %>%
        purrr::pmap_dfr(function(b_Intercept, draw_id){
          tibble::tibble(
            compound = compound,
            log_dose = seq(min_log_dose, max_log_dose, length.out=2),
            score = brms::inv_logit_scaled(b_Intercept),
            draw_id = draw_id)
        })
    })
  
  hill_model_draws <- flat_model %>%
    dplyr::left_join(dose_ranges, by="compound") %>%
    purrr::pmap_dfr(function(compound, model, min_log_dose, max_log_dose, ...){
      cat("Fittting summary for compound '", compound, "'\n", sep="")
      model %>%
        brms:::posterior_samples(pars=c("Intercept")) %>%
        dplyr::sample_n(50) %>%
        dplyr::mutate(draw_id = dplyr::row_number()) %>%
        purrr::pmap_dfr(function(b_top_Intercept, b_ic50_Intercept, b_hill_Intercept, draw_id){
          tibble::tibble(
            compound = compound,
            log_dose = seq(min_log_dose, max_log_dose, length.out=2),
            score = hill_model_forward(b_top_),
            draw_id = draw_id)
        })
    })
  
  
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
      data=flat_model_draws,
      mapping=ggplot2::aes(
        x=log_dose,
        y=score,
        group=draw_id),
      alpha=.4) +
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
