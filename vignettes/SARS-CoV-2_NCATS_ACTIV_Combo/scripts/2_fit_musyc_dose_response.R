

library(plyr)
library(tidyverse)
library(MPStats)

field_scores <- readr::read_tsv("intermediate_data/field_scores.tsv")

well_scores <- field_scores %>%
    dplyr::filter(!is.na(drug_combo)) %>%
    dplyr::group_by(plate_id, drug_combo, well_id) %>%
    dplyr::summarize(
        dose1 = SampleID1_stock_conc_in_uM,
        dose2 = SampleID1_stock_conc_in_uM,
        n_positive = sum(syn_nucs_count),
        count = sum(nuclei_count),
        .groups = "drop")

source("../../R/fit_MuSyC_score_by_dose.R")
model <- well_scores %>%
    fit_MuSyC_score_by_dose(
        group_vars = vars(drug_combo),
        verbose = TRUE,
        stan_model_args = list(verbose = TRUE))

