 
#' Drug Synergy
#' Musyc Drug Synergy model
#'
#' Assume that the response metric decreases with more effective drugs
#' Let E3 be the effect at the maximum concentration of both drugs
#' 
#' 
#' Special cases:
#'    * dose additive model: alpha1 = alpha2 = 0
#'        * loewe: h1 = h2 = 1
#'            * CI:  E0 = 1, E1 = E2 = E3 = 0
#'                   the drug effect is equated with percent inhibition
#'    * bliss drug independence model:
#'        E0 = 1, E1 = E2 = E3 = 0, alpha1 = alpha2 = 1
#' @param d1 Dose of drug 1
#' @param d2 Dose of drug 2
#' 
#' @param E0 effect with no drug treatment
#' 
#' # params for drug 1 by it self
#' @param h1 drug 1 hill slope
#' @param C1 drug 1 EC50
#' @param E1 drug 1 maximum effect
#' 
#' # params for drug 2 by it self
#' @param h2 drug 2 hill slope
#' @param C2 drug 2 EC50
#' @param E2 drug 2 maximum effect
#' 
#' @param beta synergistic efficacy
#'    percent increase in a drug combination’s effect
#'    beyond the most efficacious single drug.
#'     
#'    beta > 0 => synergistic efficacy
#'      the effect at the maximum concentration of both drugs (E3) exceeds the
#'      maximum effect of either drug alone (E1 or E2)
#'          
#'    beta < 0 => antagonistic efficacy
#'      at least one or both drugs are more efficacious as
#'      single agents than in combination
#'      
#' @param alpha1 synergistic potency 
#'    how the effective dose of drug 1
#'    is altered by the presence of drug 2
#' @param alpha2 synergistic potency 
#'    how the effective dose of drug 2
#'    is altered by the presence of drug 1
#'    
#'    alpha > 1 => synergistic potency
#'      the EC50 decreases because of the addition of the other drug,
#'      corresponding to an increase in potency
#'      
#'    0 <= alpha < 1 => antagonistic potency
#'      the EC50 of the drug increases as a result of the other drug,
#'      corresponding to a decrease in potency
#'
#'    alpha1 == alpha2 if detailed balance
#' @export
generate_MuSyC_effects <- function(
  d1,
  d2,
  E0,
  h1, C1, E1, 
  h2, C2, E2,
  alpha,
  E3 = NULL,
  beta = NULL) {

  if(!is.null(beta)){
    E3 <- min(E1, E2) - beta * min(E1, E2)
  } else if(is.null(E3)){
    stop("either E3 or beta must be non-null")
  }

  numerator <-
    C1^h1 * C2^h2 * E0 +
    d1^h1 * C2^h2 * E1 +
    C1^h1 * d2^h2 * E2 +
    d1^h1 * d2^h2 * E3 * alpha
  denominator <- 
      C1^h1 * C2^h2 +
      d1^h1 * C2^h2 +
      C1^h1 * d2^h2 +
      d1^h1 * d2^h2 * alpha
  response <- numerator / denominator
}


fit_MuSyC_score_by_dose <- function(
  well_scores,
  group_vars = vars(compound),                                  
  C1_prior = brms::prior(normal(0, 5), nlpar = "C1"),
  C2_prior = brms::prior(normal(0, 5), nlpar = "C2"),
  h1_prior = brms::prior(normal(-1, 1), nlpar = "h1", ub = -.1),
  h2_prior = brms::prior(normal(-1, 1), nlpar = "h2", ub = -.1),
  log10_alpha = brms::prior(normal(0, 2), nlpar = "alpha", lb=0),
  E0_prior = brms::prior(beta(1, 1), nlpar = "E0", lb = 0, ub = 1),
  E1_prior = brms::prior(beta(1, 1), nlpar = "E1", lb = 0, ub = 1),  
  E2_prior = brms::prior(beta(1, 1), nlpar = "E2", lb = 0, ub = 1),
  E3_prior = brms::prior(beta(1, 1), nlpar = "E3", lb = 0, ub = 1),
  C1_init = function(){rnorm(1, 0.5, 5)},
  C2_init = function(){rnorm(1, 0.5, 5)},
  h1_init = function(){rnorm(1, -1, 1)},
  h2_init = function(){rnorm(1, -1, 1)},
  log10_alpha_init = function(){rnorm(0, 2)},
  E0_init = function(){rbeta(1, 1, 1)},
  E1_init = function(){rbeta(1, 1, 1)},  
  E2_init = function(){rbeta(1, 1, 1)},
  E3_init = function(){rbeta(1, 1, 1)},
  combine = FALSE,
  ...){

  dots <- list(...)

  if(is.data.frame(well_scores)){
    grouped_data <- well_scores %>%
      dplyr::group_by(!!!group_vars) %>%
      tidyr::nest()
  }

  model <- brms::brm_multiple(
    formula = brms::brmsformula(
      n_positive | trials(cell_count) ~ log10(
        C1^h1 * C2^h2 * E0 +
        d1^h1 * C2^h2 * E1 +
        C1^h1 * d2^h2 * E2 +
        d1^h1 * d2^h2 * E3 * 10^log10_alpha
      ) - log10(
        C1^h1 * C2^h2 +
        d1^h1 * C2^h2 +
        C1^h1 * d2^h2 +
        d1^h1 * d2^h2 * 10^log10_alpha),
      C1 + C2 + h1 + h2 + log10_alpha + E0 + E1 + E2 + E3 ~ 1,
      nl = TRUE),
    data = grouped_data,
    prior = c(
      C1_prior,
      C2_prior,
      h1_prior,
      alpha_prior,
      E0_prior,
      E1_prior,
      E2_prior,
      E3_prior),
    inits = function(){
      list(
          C1 = as.array(C1_init()),
          C2 = as.array(C2_init()),          
          h1 = as.array(h1_init()),
          h2 = as.array(h2_init()),
          alpha = as.array(alpha_init()),
          E0 = as.array(E0_init()),
          E1 = as.array(E1_init()),          
          E2 = as.array(E2_init()),
          E3 = as.array(E3_init()))},
    combine = FALSE,
    dots)

  # evalate fits
  model <- model %>%
    purrr::imap(function(model, i){
      group_index <- grouped_data[i,] %>% dplyr::select(-data)
      group_index_label <- paste0(names(group_index), ": ", group_index, collapse = ", ")
      cat("Evaluating model fit for ", group_index_label, "...\n", sep = "")
      model <- model %>% brms::add_criterion(
        criterion=c("loo", "bayes_R2"),
        model_name=paste0("MuSyC:", group_index_label),
        reloo=TRUE)
      model
    })

  grouped_data %>%
    dplyr::mutate(
      model = model)
  }
  
