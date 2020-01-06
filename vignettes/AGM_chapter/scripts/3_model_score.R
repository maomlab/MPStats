library(plyr)
library(tidyverse)


load("intermediate_data/well_scores.Rdata")


model <- lm(
  prob_positive ~ week_id +  plate_ind + column,
  data=well_scores %>% dplyr::filter(!control))

summary(model)
summary(model) %>%
  readr::write_lines(path="product/well_score_by_batch_vars_model_summary_191126.txt")

data <- well_scores %>%
  dplyr::filter(!control) %>%
  dplyr::select(
    week_id, plate_ind, column,        # batch variables
    prob_positive, prob_negative) %>%  # response values
  tidyr::pivot_longer(
    cols=c("week_id", "plate_ind", "column"),
    names_to="batch_variable",
    values_to="batch_value")
    
p <- ggplot2::ggplot(data=data) +
  ggplot2::theme_bw() +
  ggplot2::geom_jitter(
    mapping=ggplot2::aes(
      x=batch_value,
      y=prob_positive),
    size=.8,
    height=0,
    width=.2) +
  ggplot2::geom_smooth(
    mapping=ggplot2::aes(
      x=batch_value,
      y=prob_positive),
    method="lm") +
  ggplot2::facet_wrap(
    facets=vars(batch_variable),
    scales="free_x") +
  ggplot2::scale_x_continuous(
    "Batch Value") +
  ggplot2::scale_y_continuous(
    "Probability of Taxol-like phenotype",
    limits=c(0,1),
    labels=scales::percent_format()) +
  ggplot2::ggtitle(
    label="Classifier Score by Batch Dimensions",
    subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
p
ggplot2::ggsave(
  filename="product/score_by_batch_var_191126.pdf",
  width=8, height=5)
ggplot2::ggsave(
  filename="product/score_by_batch_var_191126.png",
  width=8, height=5)

