
#' Customized trace plot
#'
#' Customise bayesplot::mcmc_traceplot by
#'   * adding a red smoothing line for each trace.
#'   * adding a title that uses the model$name field
#'
#' @param model brmsfit model
#' @return ggplot2::ggplot object
#'
#' @export
traceplot <- function(model, ...) {
  bayesplot::mcmc_trace(model, ...) +
    ggplot2::geom_smooth(
      mapping = ggplot2::aes(x = iteration, y = value, group = chain),
      color = "red",
      size = 1,
      method = 'gam',
      formula = y ~ s(x, bs = "cs")) +
    ggplot2::ggtitle(
      label = paste0("Trace plot: ", model$name, " model"))
}


#' Customized rank plot
#'
#' Customise bayesplot::mcmc_rank_overlay by
#'   * adding a title that uses the model$name field
#'
#' @param model brmsfit model
#' @return ggplot2::ggplot object
#'
#' @export
rankplot <- function(model) {
  bayesplot::mcmc_rank_overlay(model, ...) +
    ggplot2::ggtitle(
      label = paste0("Rank plot: ", model$name, " model"))
}

#' Customized pairs plot
#'
#' Customise pairs by
#'   * adjusting the default off-diagonal points to be smaller and more transperent
#'   * adding a title that uses the model$name field
#'
#' @param model brmsfit model
#' @return ggplot2::ggplot object
#'
#' @export
pairsplot <- function(model, size = 0.2, alpha = 0.1, ...) {
  model %>%
  pairs(off_diag_args = list(
    size = size,
    alpha = alpha),
    ...) %>%
    ggplot2::ggtitle(
      label = paste0("Pairs plot: ", model$name, " model"))
}

#' Summarize the loo_R2 as a plot indicator
#'
#'   Produces a label on the plot using the coordinates of the viewport rather than
#'   the plot coordinates. The indicator is with the indicated number of digits as sig figs.
#'
#'             R^2 = mean [<Q2.5>, <Q97.5>]
#' 
#' @param model brmsfit
#' @param digits number of significant figure digits
#' @return MPStats::geom_indicator object to that can be added to a ggplot2::ggplot
#'
#' @export
loo_R2_indicator <- function(model, digits = 2, ...) {
  R2_scores <- data.frame(
    indicator = paste0(
        'R^2~"="~',
        model$criteria$loo_R2 %>% mean() %>% signif(digits),
        '~" ["*', model$criteria$loo_R2 %>% quantile(0.025) %>% signif(digits),
        '*", "*', model$criteria$loo_R2 %>% quantile(0.975) %>% signif(digits),
        '*"]"'))
  geom_indicator(
    data = R2_scores,
    mapping = ggplot2::aes(indicator = indicator),
    parse = TRUE,
    ...)
}

#' Create a plot that compares the distribution of prior and posterior
#' marginal distributions for each variable
#'
#' To make the plot managable, the distribution for each variable
#' is trimmed at the 95% IQR
#' 
#' The plot title uses `model$name`
#' 
#' @param model brmsfit model
#' @return ggplot2::ggplot object
#'
#' @export
prior_posterior_plot <- function(model) {
  model_prior <- model %>%
    brms:::update.brmsfit(sample_prior = "only")

  draws <- dplyr::bind_rows(
    model %>%
      tidybayes::tidy_draws() %>%
      tidybayes::gather_variables() %>%
      dplyr::mutate(sample_type = "Posterior"),
    model_prior %>%
      tidybayes::tidy_draws() %>%
      tidybayes::gather_variables() %>%
      dplyr::mutate(sample_type = "Prior")) %>%
    dplyr::filter(!stringr::str_detect(.variable, "__$")) %>%
    dplyr::group_by(.variable) %>%
      dplyr::filter(
        .value < quantile(.value, 0.975),
        .value > quantile(.value, 0.025)) %>%
    dplyr::ungroup()

  ggplot2::ggplot(data = draws) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_density(
      mapping = ggplot2::aes(
        x = .value,
        fill = sample_type,
        group = sample_type),
      color = "black",
      alpha = .9) +
    ggplot2::ggtitle(
      label = paste0("Prior vs Posterior Marginal Distribution: ", model$name, " model")) +
    ggplot2::facet_wrap(
      facets = dplyr::vars(.variable),
      scales = "free") +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::scale_x_continuous("Parameter Value") +
    ggplot2::scale_fill_discrete("Distribution")
}
