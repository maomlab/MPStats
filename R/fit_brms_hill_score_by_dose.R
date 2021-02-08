#' Simulate hill model from the parameters
#'
#' @param log_dose log10(doses) at which to evaluate the model
#' @param top top parameter
#' @param bottom bottom parameter
#' @param ic50 ic50 parameter
#' @param hill hill parameter
#' @return vector of values
#'
#' @export
generate_hill_effects <- function(
  log_dose,
  top,
  bottom,
  ic50,
  hill) {

  brms::inv_logit_scaled(
    x = hill * 4 / top * (log_dose - ic50),
    lb = bottom,
    ub = top)

}


#' brms hill model for score by compound dose
#'
#' Use the R BRMS package to fit dose response data to with a sigmoid hill function with bottom = 0
#'
#' @param well_scores tibble::tibble as output by read_well_scores
#' @return list of brms::brmsfit objects one for each compound
#'
#' @export
fit_brms_hill_score_by_dose <- function(scores){

  grouped_data <- scores %>%
    dplyr::filter(!is_control) %>%
    dplyr::group_by(compound) %>%
    tidyr::nest()

  model <- brms::brm_multiple(
    formula = brms::brmsformula(
      n_positive | trials(cell_count) ~ top * inv_logit(hill*4/top*(log_dose - ic50)),
      top + ic50 + hill ~ 1,
      nl=TRUE),
    data=grouped_data$data,
    prior = c(
      brms::prior(uniform(0, 1), nlpar="top", lb=0, ub=1),
      brms::prior(normal(2, 10), nlpar="ic50"),
      brms::prior(normal(1, 10), nlpar="hill")),
    family=binomial("identity"),
    inits=function(){
      list(
        top=as.array(runif(1, 0, 1)),
        ic50=as.array(rnorm(1, 2, 10)),
        hill=as.array(rnorm(1, 1, 10)))},
    iter=8000,
    control=list(
      adapt_delta=0.99,
      max_treedepth=12),
    combine=FALSE)

  # evalate fits
  model <- model %>%
    purrr::imap(function(model, i){
      compound <- group_data$compound[i]
      cat("Evaluating model fit for compound '", compound, "' ...\n", sep="")
      model <- model %>% brms::add_criterion(
        criterion=c("loo", "bayes_R2"),
        model_name=paste0("hill_", compound),
        reloo=TRUE)
      model
    })

  tibble::tibble(
    compound=grouped_data$compound,
    model=model)

}
################

#' brms model for the counts of positive and negative cells by compound dose
#'
#' Use the R BRMS package to fit dose response data to with a multiple poisson model
#'
#' W. Scott Comulada and Robert E. Weiss, On Models for Binomial Data with Random Numbers of Trials
#' Biometrics, 2007, 63(2): 610–617. doi:10.1111/j.1541-0420.2006.00722.x.
#'
#' @param well_scores tibble::tibble as output by read_well_scores
#' @return list of brms::brmsfit objects one for each compound
#'
#' @export
fit_brms_npos_nneg_by_dose <- function(well_scores){

  dose_response_data_by_compound <- well_scores %>%
    dplyr::mutate(n_negative = cell_count - n_positive) %>%
    dplyr::group_by(compound) %>%
    tidyr::nest()

  hill_model <- brms::brm_multiple(
    formula = brms::brmsformula(
      brms::mvbind(n_positive, n_negative) ~ (top-bottom) * inv_logit(hill*4/(top-bottom)*(log_dose - ic50)) + bottom,
      top + bottom + ic50 + hill ~ 1,
      nl=TRUE),
    data=dose_response_data_by_compound$data,
    prior = c(
      brms::prior(uniform(0, 1), nlpar="top", lb=0, ub=1),
      brms::prior(normal(2, 10), nlpar="ic50"),
      brms::prior(normal(1, 10), nlpar="hill")),
    family=binomial("identity"),
    inits=function(){
      list(
        top=as.array(runif(1, 0, 1)),
        ic50=as.array(rnorm(1, 2, 10)),
        hill=as.array(rnorm(1, 1, 10)))},
    iter=8000,
    control=list(
      adapt_delta=0.99,
      max_treedepth=12),
    combine=FALSE)

  # assign the compound id to each model
  hill_model <- hill_model %>%
    purrr::imap(function(model, i){
      compound <- dose_response_data_by_compound$compound[i]
      cat("Evaluating model fit for compound '", compound, "' ...\n", sep="")
      model$compound <- compound
      model <- model %>% brms::add_criterion(
        criterion=c("loo", "bayes_R2"),
        model_name=paste0("hill_", compound),
        reloo=TRUE)
      model
    })
  hill_model
}
