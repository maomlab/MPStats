library(plyr)
library(tidyverse)


load("intermediate_data/well_scores.Rdata")


model <- lm(
  sqrt(well_cell_count) ~ week_id +  plate_ind + column,
  data=well_scores %>% dplyr::filter(!control))
summary(model)
summary(model) %>%
  readr::write_lines(path="product/well_cell_count_by_batch_vars_191126.txt")


model_week_factor <- lm(
  sqrt(well_cell_count) ~ factor(week_id) +  plate_ind + column,
  data=well_scores %>% dplyr::filter(!control))
summary(model_week_factor)

summary(model_week_factor) %>%
  readr::write_lines(path="product/well_cell_count_by_batch_vars_week_factor_191126.txt")


data <- well_scores %>%
  dplyr::filter(!control) %>%
  dplyr::select(
    week_id, plate_ind, column,        # batch variables
    well_cell_count) %>%               # response value
  tidyr::pivot_longer(
    cols=c("week_id", "plate_ind", "column"),
    names_to="batch_variable",
    values_to="batch_value")

p <- ggplot2::ggplot(data=data) +
  ggplot2::theme_bw() +
  ggplot2::geom_jitter(
    mapping=ggplot2::aes(
      x=batch_value,
      y=sqrt(well_cell_count)),
    size=.8,
    height=0,
    width=.2) +
  ggplot2::geom_smooth(
    mapping=ggplot2::aes(
      x=batch_value,
      y=sqrt(well_cell_count)),
    method="lm") +
  ggplot2::facet_wrap(
    facets=vars(batch_variable),
    scales="free_x") +
  ggplot2::scale_x_continuous(
    "Batch Value") +
  ggplot2::scale_y_continuous(
    "Well Cell Count",
    breaks=c(0,10,20,30,40),
    labels=c(0,100,400,900,1600)) +
  ggplot2::ggtitle(
    label="Well Cell Count by Batch Dimensions",
    subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")

ggplot2::ggsave(
  filename="product/well_cell_count_by_batch_var_191126.pdf",
  width=8, height=5)
ggplot2::ggsave(
  filename="product/well_cell_count_by_batch_var_191126.png",
  width=8, height=5)



#####################
p <- ggplot2::ggplot(data=data) +
  ggplot2::theme_bw() +
  ggplot2::geom_density(
    mapping=ggplot2::aes(
      x=sqrt(well_cell_count),
      color=batch_value,
      group=batch_value)) +
  ggplot2::facet_wrap(
    facets=vars(batch_variable)) +
  ggplot2::scale_x_continuous(
    "Well Cell Count",
    breaks=c(0, 10, 20, 30, 40),
    labels=c(0, 100, 400, 900, 1600)) +
  ggplot2::scale_y_continuous("Density") +
  ggplot2::scale_color_continuous("Batch Value") +
  ggplot2::ggtitle(
    label="Well Cell Count Density by Batch Dimensions",
    subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)") +
  ggplot2::theme(
    legend.position=c(.84,1.13),
    legend.direction="horizontal",
    legend.margin=ggplot2::margin(t=0,r=0,b=0,l=0),
    legend.box.background=ggplot2::element_blank())
p
ggplot2::ggsave(
  filename="product/well_cell_count_density_by_batch_var_191126.pdf",
  width=8, height=5)
ggplot2::ggsave(
  filename="product/well_cell_count_density_by_batch_var_191126.png",
  width=8, height=5)




