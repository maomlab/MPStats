library(plyr)
library(tidyverse)
library(MPStats)
library(brms)

cat("Modeling score by batch variables\n")

load("intermediate_data/well_scores.Rdata")


model <- MPStats::model_score_by_batch_vars_lm(well_scores)
summary(model) %>%
  readr::write_lines(
    path=paste0("product/model_score_by_batch_vars_summary_", MPStats::date_code(), ".txt"))


compound_moa %>%
  plyr::d_ply("moa2", function(compounds){
    moa_type = compounds$moa2[1]
    cat("Fitting batch variables conditional on known MOA type: '", moa_type, "' ...\n", sep="")
    scores <- well_scores %>% dplyr::semi_join(compounds, by="compound")
    model <- MPStats::model_score_by_batch_vars_lm(well_scores=scores)
    model %>% summary() %>% print()
  })


plot <- MPStats::plot_score_by_batch_vars(
  well_scores=well_scores,
  subtitle="Human MCF7 cells -- compound-profiling experiment (BBBC021v1)")
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/score_by_batch_vars_", MPStats::date_code(), ".pdf"),
  width=8, height=5)
ggplot2::ggsave(
  plot=plot,
  filename=paste0("product/score_by_batch_vars_", MPStats::date_code(), ".png"),
  width=8, height=5)