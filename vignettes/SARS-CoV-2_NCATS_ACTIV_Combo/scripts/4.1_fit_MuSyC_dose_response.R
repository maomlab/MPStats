library(plyr)
library(tidyverse)
library(MPStats)
library(brms)


use_future <- FALSE

use_future <- TRUE
library(future)
Sys.setenv(DEBUGME = "batchtools")
library(batchtools)
library(future.batchtools)

future::plan(
    list(
        # for each drug combo
        future::tweak(
            future.batchtools::batchtools_slurm,
                resources = list(
                    account = "maom99",
                    ntasks = 1,
                    ncpus = 1L,
                    memory = "10GB"),
            template = "../../inst/batchtools.greatlakes.tmpl",
            registry = list(
                packages = "MPStats")),
        future::tweak(
            future.batchtools::batchtools_slurm,
                resources = list(
                    account = "maom99",
                    ntasks = 1,
                    ncpus = 1L,
                    memory = "2GB"),
            template = "../../inst/batchtools.greatlakes.tmpl",
            registry = list(
                packages = "MPStats"))))


########################################
# fit the MuSyC synergy model separate #
########################################


well_scores <- readr::read_tsv("intermediate_data/well_scores.tsv") %>%
    dplyr::mutate(
        dose1_nM = dose1,
        dose2_nM = dose2,
        dose1 = dose1_nM * 1e-9,
        dose2 = dose2_nM * 1e-9,
        d1_scale_factor = 1e-6,
        d2_scale_factor = 1e-6,
        sample1 = sample_label_1,
        sample2 = sample_label_2)

devtools::load_all()
synergy_model_combined_v1 <- well_scores %>%
    MPStats::fit_MuSyC_score_by_dose_robust(
        formula = "MuSyC_separate",
        group_vars = dplyr::vars(drug_combo),
        stan_model_args = list(verbose = TRUE),
        open_progress = FALSE,
        silent = FALSE,
        future = use_future)

save(synergy_model_v10, file = "intermediate_data/synergy_model_v10.Rdata")
load("intermediate_data/synergy_model_v10.Rdata")


estimated_parameters <- synergy_model_v10 %>%
    dplyr::rowwise() %>%
    dplyr::do({
        combo_model <- .
        combo_model$model %>%
            tidybayes::spread_draws(
                E0,
                E1, C1, s1,
                E2, C2, s2,
                E3, alpha) %>%
            tidybayes::median_qi() %>%
            dplyr::mutate(
                drug_combo = combo_model$drug_combo)
    }) %>%
    dplyr::ungroup()
estimated_parameters %>%
    dplyr::transmute(
        drug_combo,
        E0 = E0 * 100,
        C1, E1 = E1 * 100, s1,
        C2, E2 = E2 * 100, s2,
        E3 = E3 * 100,
        alpha) %>%
    data.frame

estimated_parameters %>%
    readr::write_tsv("product/estimated_parameters_20201127.tsv")


########################################
# fit the MuSyC synergy model combined #
########################################
DMSO_recode_sample1 <- "GS-441524"
DMSO_recode_sample2 <- "Tizoxanide"

well_scores <- readr::read_tsv("intermediate_data/field_scores.tsv") %>%
    dplyr::mutate(
        dose1 = SampleID1_stock_conc_in_uM * 1e-9,
        dose2 = SampleID2_stock_conc_in_uM * 1e-9) %>%
    dplyr::group_by(
        plate_id,
        well_id,
        sample_label_1,
        sample_label_2,
        dose1,
        dose2) %>%
    dplyr::summarize(
        n_positive = sum(syn_nucs_count),
        count = sum(nuclei_count),
        score = n_positive / count,
        .groups = "drop") %>%
    dplyr::mutate(
        sample1 = ifelse(sample_label_1 == "DMSO", DMSO_recode_sample1, sample_label_1),
        sample2 = ifelse(sample_label_2 == "DMSO", DMSO_recode_sample2, sample_label_2),
        d1_scale_factor = 1e-6,
        d2_scale_factor = 1e-6)

devtools::load_all()
synergy_model_combined_od_v1 <- well_scores %>%
    MPStats::fit_MuSyC_score_by_dose_robust(
        formula = "MuSyC_combined_od",
        stan_model_args = list(verbose = TRUE),
        open_progress = FALSE,
        control = list(),
        silent = FALSE,
        future = use_future)

# seems to work running loo now
# rel error vs predicted % infected shows some of the error can be
# attributed to plate batch effects
save(synergy_model_combined_v3, file = "intermediate_data/synergy_model_combined_v3.Rdata")

# seems to work running loo now
# rel error vs predicted % infected shows some of the error can be
# attributed to plate batch effects
save(synergy_model_combined_v4, file = "intermediate_data/synergy_model_combined_v4.Rdata")

# same as above but modeling zero inflated counts
# rel error vs predicted % infected shows some of the error can be
# attributed to plate batch effects
save(synergy_model_combined_zi_v1, file = "intermediate_data/synergy_model_combined_zi_v1.Rdata")

save(synergy_model_combined_od_v1, file = "intermediate_data/synergy_model_combined_od_v1.Rdata")

