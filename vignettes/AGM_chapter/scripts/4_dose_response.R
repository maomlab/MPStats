# dose resposne effects

library(plyr)
library(tidyverse)
library(drc)

# this is a custom ggplot2 element that allows putting the p-values directly on each facet
source("scripts/geom_indicator.R")

load("intermediate_data/well_scores.Rdata")

dose_response_data <- well_scores %>%
  dplyr::filter(!control) %>%
  dplyr::transmute(
    compound = compound,
    log_dose = log10(concentration*1000),
    value = prob_positive,
    well_cell_count=well_cell_count,
    variance = .25/(1+well_cell_count))

# DRC 4-param dose response fits with bottom fixed = 0
fits <- dose_response_data %>%
  plyr::ddply(c("compound"), function(curve_data){
    tryCatch({
      fit <- drc::drm(
        formula=value ~ log_dose,
        weights=1/curve_data$variance,
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
      cat("Failed to fit curve for compound: ",curve_data$compound[1], "\n", sep="")
      return(data.frame())
    })
  })


dr_summary <- dose_response_data %>%
  dplyr::group_by(log_dose, compound) %>%
  dplyr::summarize(
    value = mean(value),
    # the variance of an average of iid variables:
    variance = sum(variance) / dplyr::n()^2)


p <- ggplot2::ggplot(data=dr_summary) +
  ggplot2::theme_bw() +
  ggplot2::geom_line(
    data=fits,
    mapping=ggplot2::aes(
      x=log_dose,
      y=pred_value)) +
  ggplot2::geom_errorbar(
    data=dr_summary,
    mapping=ggplot2::aes(
      x=log_dose,
      ymin=value - sqrt(variance),
      ymax=value + sqrt(variance))) +
  ggplot2::geom_point(
    data=dr_summary,
    mapping=ggplot2::aes(
      x=log_dose,
      y=value)) +
  geom_indicator(
    data=fits %>% dplyr::distinct(compound, p_value),
    mapping=ggplot2::aes(
      indicator=paste0("fit: ", p_value)),
    xpos="left",
    ypos="top",
    group=1) +
  ggplot2::ggtitle("BBBC021 screen") +
  ggplot2::scale_x_continuous(
    "log[Compound dose] (nM)") +
  ggplot2::scale_y_continuous(
    "Well-score (% taxol-like)",
    limits=c(0,1),
    labels=scales::percent_format()) +
  ggplot2::facet_wrap(~compound, scales="free_x")


ggplot2::ggsave(filename="product/dose_response_curves_191125.pdf", height=20, width=20)
ggplot2::ggsave(filename="product/dose_response_curves_191125.png", height=20, width=20)


fits_summary <- fits %>%
  dplyr::distinct(compound, .keep_all=T) %>%
  dplyr::select(compound, top, bottom, ic50, slope, p_value)



###
max_value <- dose_response_data %>%
  dplyr::group_by(compound) %>%
  dplyr::summarize(max_value=max(value)) %>%
  dplyr::arrange(max_value) %>%
  data.frame

dose_response_data <- dose_response_data %>%
  dplyr::semi_join(max_value %>% dplyr::filter(max_value > .5), by="compound")

fits <- fits %>%
  dplyr::semi_join(max_value %>% dplyr::filter(max_value > .5), by="compound")

p <- ggplot2::ggplot(data=dose_response_data) +
  ggplot2::theme_bw() +
  ggplot2::geom_line(
    data=fits,
    mapping=ggplot2::aes(
      x=log_dose,
      y=pred_value),
    size=1.5,
    color="blue") +
  ggplot2::geom_point(
    mapping=ggplot2::aes(
      x=log_dose,
      y=value),
    size=.8) +
  ggplot2::scale_x_continuous(
    "log[Compound dose] (nM)") +
  ggplot2::scale_y_continuous(
    "Activity",
    limits=c(0,1),
    labels=scales::percent_format()) +
  ggplot2::facet_wrap(~compound, scales="free_x")

ggplot2::ggsave(filename="product/dose_response_curves_clean_191214.pdf", height=10, width=10)
ggplot2::ggsave(filename="product/dose_response_curves_clean_191214.png", height=10, width=10)
