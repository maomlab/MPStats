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

#' Fit the MuSyC synergy model by dose
#'
#' @param well_scores data.frame
#'        with columns: dose1, dose2, n_positive, count, [<group_vars>]
#' @param group_vars quosures list
#'        dplyr::vars(...) columns to over when fitting synergy model
#' @param C1_prior prior distribution for Ed when d1=d1_IC50, d2=0
#' @param C2_prior prior distribution for Ed when d1=0, d2=d2_IC50
#' @param s1_prior prior distribution for d(Ed)/d(d1) when d1=d1_IC50, d2=0
#' @param s2_prior prior distribution for d(Ed)/d(d2) when d1=0, d2=d2_IC50
#' @param alpha_prior prior distribution for alpha synergy parameter
#' @param E0_prior prior distribution for Ed when d1=0, d2=0
#' @param E1_prior prior distribution for Ed when d1=Inf, d2=0
#' @param E2_prior prior distribution for Ed when d1=0, d2=Inf
#' @param E3_prior prior distribution for Ed when d1=Inf, d2=Inf
#' @param C1_init initial sampling distribution for the C1 parameter
#' @param C2_init initial sampling distribution for the C2 parameter
#' @param s1_init initial sampling distribution for the s1 parameter
#' @param s2_init initial sampling distribution for the s2 parameter
#' @param alpha_init intial sampling distribution for the alpha parameter
#' @param E0_init initial sampling distribution for the E0 parameter
#' @param E1_init initial sampling distribution for the E1 parameter
#' @param E2_init initial sampling distribution for the E2 parameter
#' @param E3_init initial sampling distribution for the E3 parameter
#' @param combine combine the grouped models into a single brms model
#' @param verbose verbose output
#'
#' @param iter number of stan NUTS sampling steps
#' @param cores number of cores used for sampling
#' @param stan_model_args stan model arguments
#' @param control stan control arguments
#'
#' The
#'
#'    bernoulli_inf(n_positive / count) = Ed ~ MuSyC(d1, d2, C_params, E_params, s_params, alpha)
#'
#'    To improve numeric stability, the d1 and d2 and C1 and C2 variables
#'    are scaled to improve numeric stability:
#'
#'       d1 = dose1/max(dose1)
#'       d2 = dose2/max(dose2)
#'       drug1_IC50 = C1 * max(dose1)
#'       drug2_IC50 = C2 * max(dose2)
#'
#' Functional form:
#' Ed ~ (
#'        C1^h1 * C2^h2 * E0 +
#'        d1^h1 * C2^h2 * E1 +
#'        C1^h1 * d2^h2 * E2 +
#'        d1^h1 * d2^h2 * E3 * alpha
#'      ) / (
#'        C1^h1 * C2^h2 +
#'        d1^h1 * C2^h2 +
#'        C1^h1 * d2^h2 +
#'        d1^h1 * d2^h2 * alpha
#'      )
#'
#'
#'
#'
#'
#' ##############################################
#' # Proof of the definitions of the parameters #
#' ##############################################
#'
#' Claim: When d1=0 and d2=0 then Ed = E0
#' Ed = (C1^h1 * C2^h2 * E0) / (C1^h1 * C2^h2)
#'    = E0
#'
#' Claim: When d1=0 and d2 -> Inf then Ed = E2
#' Ed = (C2^h2 * E0 + d2^h2 * E2) / (C2^h2 + d2^h2)
#'    = (d2^h2 * E2) / (d2^h2)
#'    = E2
#'
#; Claim: When d1=0 and d2=C2 then Ed = (E0 + E2) / 2
#' When d1>0 and d2 -> Inf then Ed
#' Ed = (C1^h1 * C2^h2 * E0 + C1^h1 * C2^h2 * E2) / (C1^h1 * C2^h2 + C1^h1 * C2^h2)
#'    = (E0 + E2) / 2
#'
#' Claim: When d1=0 and d2=C2 then d(Ed)/d(d2) = s2
#'        where h2 = s2 * (4 * C2) / (E0 + E2))
#'
#' d(Ed)/d(d2)
#'   =  d/d(d2)
#'      (C1^h1 * C2^h2 * E0 + C1^h1 * d2^h2 * E2) /
#'      (C1^h1 * C2^h2      + C1^h1 * d2^h2)
#'
#'Cancle the C1^h1 terms:
#'   =  d/d(d2)
#'      (C2^h2 * E0 + d2^h2 * E2) /
#'      (C2^h2      + d2^h2)
#'
#'
#' distribute the derivative across the terms in the numerator
#'   =  E0 * C2^h2 * [d/d(d2) 1     / (C2^h2 + d2^h2)] +
#'      E2         * [d/d(d2) d2^h2 / (C2^h2 + d2^h2)]
#'
#'   =  E0 * C2^h2 * [h2 * d2^(h2-1) / (C2^h2 + d2^h2)^2] +
#'      E2 * [C2^h2 * h2 * d2^(h2-1) / (C2^h2 + d2^h2)^2]
#'
#'   =  (E0 + E2) * C2^h2 * h2 * d2^(h2-1)/(C2^h2 + d2^h2)^2
#'
#' Evaluate at d2 = C2:
#'   =  (E0 + E2) * h2 * C2^(2*h2-1) / [4*C2^(2*h2))]
#'   =  h2 * (E0 + E2) / (4 * C2)


fit_MuSyC_score_by_dose <- function(
  well_scores,
  group_vars = vars(compound),
  C1_prior = brms::prior(normal(0.5, 0.5), nlpar = "C1", lb = 0),
  C2_prior = brms::prior(normal(0.5, 0.5), nlpar = "C2", lb = 0),
  s1_prior = brms::prior(normal(1, 3), nlpar = "s1", lb = -.1),
  s2_prior = brms::prior(normal(1, 3), nlpar = "s2", lb = -.1),
  alpha_prior = brms::prior(normal(0, 2), nlpar = "alpha", lb = 0),
  E0_prior = brms::prior(beta(1, 1), nlpar = "E0", lb = 0, ub = 1),
  E1_prior = brms::prior(beta(1, 1), nlpar = "E1", lb = 0, ub = 1),
  E2_prior = brms::prior(beta(1, 1), nlpar = "E2", lb = 0, ub = 1),
  E3_prior = brms::prior(beta(1, 1), nlpar = "E3", lb = 0, ub = 1),
  C1_init = function() {as.array(rnorm(1, 0.5, 5))},
  C2_init = function() {as.array(rnorm(1, 0.5, 5))},
  s1_init = function() {as.array(rnorm(1, .1, 1.5))},
  s2_init = function() {as.array(rnorm(1, .1, 1.5))},
  alpha_init = function() {as.array(runif(1, 0, 3))},
  E0_init = function() {as.array(rbeta(1, 1, 1))},
  E1_init = function() {as.array(rbeta(1, 1, 1))},
  E2_init = function() {as.array(rbeta(1, 1, 1))},
  E3_init = function() {as.array(rbeta(1, 1, 1))},
  combine = FALSE,
  verbose = FALSE,
  iter = 8000,
  cores = 4,
  stan_model_args = list(verbose = FALSE),
  control = list(
      adapt_delta = .99,
      max_treedepth = 12),
  ...) {

  if (is.data.frame(well_scores)) {
    grouped_data <- well_scores %>%
      dplyr::group_by(!!!group_vars) %>%
      dplyr::mutate(
        d1_scale_factor = max(dose1),
        d2_scale_factor = max(dose2)) %>%
      tidyr::nest()
  }

  if(verbose){
      cat("Fitting MuSyC model\n")
  }

  model <- brms::brm_multiple(
    formula = brms::brmsformula(
      n_positive | trials(count) ~ (
        C1^h1 * C2^h2 * E0 +
        d1^h1 * C2^h2 * E1 +
        C1^h1 * d2^h2 * E2 +
        d1^h1 * d2^h2 * E3 * 10^alpha
      ) / (
        C1^h1 * C2^h2 +
        d1^h1 * C2^h2 +
        C1^h1 * d2^h2 +
        d1^h1 * d2^h2 * 10^alpha),
      brms::nlf(d1 ~ dose1 / d1_scale_factor),
      brms::nlf(d2 ~ dose2 / d2_scale_factor),
      brms::nlf(h1 ~ s1 * (4 * C1) / (E0 + E1)),
      brms::nlf(h2 ~ s2 * (4 * C2) / (E0 + E2)),
      E0 + C1 + E1 + s1 + C2 + E2 + s2 + alpha + E3 ~ 1,
      nl = TRUE),
    data = grouped_data$data,
    family = binomial("identity"),
    prior = c(
      C1_prior,
      C2_prior,
      s1_prior,
      s2_prior,
      alpha_prior,
      E0_prior,
      E1_prior,
      E2_prior,
      E3_prior),
    inits = function() {
      list(
        b_C1 = C1_init,
        b_C2 = C2_init(),
        b_s1 = s1_init(),
        b_s2 = s2_init(),
        b_alpha = alpha_init(),
        b_E0 = E0_init(),
        b_E1 = E1_init(),
        b_E2 = E2_init(),
        b_E3 = E3_init())},
#    stanvars = c(
#        brms::stanvar(
#            scode = "  real d1_scale_factor = max(dose1));",
#            block = "tdata",
#            position = "end"),
#        brms::stanvar(
#            scode = "  real d2_scale_factor = max(dose2));",
#            block = "tdata",
#            position = "end"),
#        brms::stanvar(
#            scode = "  real drug1_IC50 = b_C1 * d1_scale_factor);",
#            block = "genquant",
#            position = "end"),
#        brms::stanvar(
#            scode = "  real drug2_IC50 = b_C2 * d2_scale_factor;",
#            block = "genquant",
#            position = "end")),
    combine = FALSE,
    data2 = NULL,
    iter = iter,
    cores = cores,
    stan_model_args = stan_model_args,
    control = control)

  # evalate fits
  model <- model %>%
    purrr::imap(function(model, i) {
      group_index <- grouped_data[i, ] %>% dplyr::select(-data)
      group_index_label <- paste0(
          names(group_index), ": ", group_index, collapse = ", ")
      cat("Evaluating model fit for ", group_index_label, "...\n", sep = "")
      model <- model %>% brms::add_criterion(
        criterion = c("loo", "bayes_R2"),
        model_name = paste0("MuSyC:", group_index_label),
        reloo = TRUE)
      model
    })

  grouped_data %>%
    dplyr::mutate(
      model = model)
  }
