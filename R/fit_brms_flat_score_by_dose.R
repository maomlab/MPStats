


#' brms flat model for score by compound dose
#'
#' Use the R BRMS package to fit dose response data to with constant function
#'  
#' @param well_scores tibble::tibble as output by read_well_scores
#' @return list of brms::brmsfit objects one for each compound
#' 
#' @export
fit_brms_flat_score_by_dose <- function(well_scores){
  
  grouped_data <- well_scores %>%
    dplyr::filter(!is_control) %>%
    dplyr::group_by(compound) %>%
    tidyr::nest()
  
  # fit the model
  model <- brms::brm_multiple(
    formula= n_positive | trials(cell_count) ~ 0 + Intercept,
    data=grouped_data$data,
    family = binomial("logit"),
    prior = c(brms::prior(normal(0, 10), coef="Intercept")),
    combine=FALSE)
  
  # evaluate model fit
  model <- model %>%
    purrr::imap(function(compound_model, i){
      compound <- grouped_data$compound[i]
      cat("Evaluating model fit for compound '", compound, "' ...\n", sep="")
      compound_model %>% brms::add_criterion(
        criterion=c("loo", "bayes_R2"),
        model_name=paste0("flat_", compound),
        reloo=TRUE)
    })
  
  tibble::tibble(
    compound = grouped_data$compound,
    model = model)
}
