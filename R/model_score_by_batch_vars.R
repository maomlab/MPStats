
#' Model score by batch variables using a weighted least squares with linear predictors
#' 
#' 
#' 
#' @param well_scores result of read_well_scores()
#' @return "lm" object of model
#' 
#' @export
model_score_by_batch_vars_lm <- function(well_scores){
  scores <- well_scores %>% dplyr::filter(!is_control)
  model <- lm(
    prob_positive ~ week +  plate + row,
    weights = 1/binomial_variance(scores$n_positive, scores$cell_count),
    data=scores)
}




