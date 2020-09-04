


#' DRC 4-param dose response fits with bottom fixed = 0
#'
#' Use the R drc package to fit dose response data
#' 
#' The model that is fit is the 4-parameter sigmoid with bottom fixed at 0.
#' To fit the curve, drc uses weighted least squares. Assuming errors are not
#' correlated the BLUE has weight w_i = 1/V_i, where V_i is the variance of the
#' bayesian probability positive estimator using a prior of Beta(10,10)
#' and then normalized so each replica has on average weight 1.
#' 
#' @param well_scores tibble::tibble as output by read_well_scores
#' @return data.frame with information about fits
#' 
#' @export
fit_drc_score_by_dose <- function(well_scores){
  fits <- well_scores %>%
    plyr::ddply(c("compound"), function(curve_data){
      tryCatch({
        weights <- 1/binomial_variance(curve_data$n_positive, curve_data$cell_count, prior_positive=10, prior_negative=10)
        weights <- weights * nrow(curve_data)/sum(weights)
        fit <- drc::drm(
          formula=prob_positive ~ log_dose,
          weights=weights,
          data=curve_data,
          fct=drc::L.4(fixed=c(NA, 0, NA, NA)))
        log_dose <- seq(min(curve_data$log_dose), max(curve_data$log_dose), length.out=100)
        pred_value <- predict(fit, expand.grid(log_dose, 1))
        test_no_fit <- drc::noEffect(fit) %>% data.frame
        data.frame(log_dose, pred_value) %>%
          dplyr::mutate(
            slope=fit$coefficients[1],
            bottom=0,
            top=fit$coefficients[2],
            ic50=fit$coefficients[3],
            chi_squared_test = test_no_fit$.[1],
            degrees_of_freedom = test_no_fit$.[2],
            p_value = test_no_fit$.[3] %>% signif(2)) %>%
          return()
      }, error=function(e){
        cat("Failed to fit curve for compound: ", curve_data$compound[1], "\n", sep="")
        return(data.frame())
      })
    })
}
