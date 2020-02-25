library(plyr)
library(tidyverse)
library(MPStats)
library(ggplot2)

cat("Analyzing Zprime scores by features and dyes\n")

load("intermediate_data/Zprime_by_plate.R")


plot <- MPStats::plot_Zprime_by_dye(Zprime_by_plate)
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/Zprime_by_dye_", MPStats::date_code(), ".pdf"),
  width=4, height=4,
  useDingbats=FALSE)
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/Zprime_by_dye_", MPStats::date_code(), ".png"),
  width=4, height=4)