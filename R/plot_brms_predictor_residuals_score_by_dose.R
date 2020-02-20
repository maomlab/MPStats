
# predictor residual plots
# counter factual plots
# posterior predictive plots



#' Plot the score as a function of compound dose for brms models
#' 
#' @param well_scores tibble::tibble as output by read_well_scores
#' @param flat_model tibble::tibble as output by fit_brms_flat_score_by_dose
#' @param hill_model tibble::tibble as output by fit_brms_hill_score_by_dose
#' @param subtitle string subtitle for plot
#' @return ggplot2 plot
#' 
#' @export
plot_brms_predictor_residual_score_by_dose <- function(
  well_scores,
  flat_model,
  hill_model,
  subtitle=NULL){
  
 
  
  
   
}