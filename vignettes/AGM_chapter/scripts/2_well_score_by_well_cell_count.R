library(plyr)
library(tidyverse)
library(ggplot2)

load("intermediate_data/well_scores.Rdata")

p <- ggplot2::ggplot(data=well_scores) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping=ggplot2::aes(
      x=sqrt(well_cell_count),
      y=prob_positive)) +
  ggplot2::geom_smooth(
    mapping=ggplot2::aes(
      x=sqrt(well_cell_count),
      y=prob_positive)) +
  ggplot2::ggtitle(
    label="Classifier Score by Well Cell Count",
    subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)") +
  ggplot2::scale_x_continuous(
    "Well Cell Count",
    breaks=c(0,10,20,30,40),
    labels=c(0,100,400,900,1600)) +
  ggplot2::scale_y_continuous(
    "Probability of Taxol-like phenotype",
    limits=c(0,1),
    labels=scales::percent_format())
p

ggplot2::ggsave(
  filename="product/score_by_well_cell_count_191127.pdf",
  width=8, height=5)
ggplot2::ggsave(
  filename="product/score_by_well_cell_count_191127.png",
  width=8, height=5)
