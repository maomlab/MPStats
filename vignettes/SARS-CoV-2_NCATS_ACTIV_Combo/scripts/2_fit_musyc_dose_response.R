library(plyr)
library(tidyverse)
library(MPStats)
library(brms)


library(future)
Sys.setenv(DEBUGME = "batchtools")
library(batchtools)
library(future.batchtools)

future::plan(
    list(
        # for each drug combo
        tweak(
            future.batchtools::batchtools_slurm,
                resources = list(
                    account = "maom99",
                    ntasks=1,
                    ncpus=1L,
                    memory="10GB"),
            template = "../../inst/batchtools.greatlakes.tmpl"),
        tweak(
            future.batchtools::batchtools_slurm,
                resources = list(
                    account = "maom99",
                    ntasks=1,
                    ncpus=1L,
                    memory="2GB"),
            template = "../../inst/batchtools.greatlakes.tmpl")))


###############################
# fit the MuSyC synergy model #
###############################


well_scores <- readr::read_tsv("intermediate_data/well_scores.tsv") %>%
    dplyr::mutate(
        dose1_nM = dose1,
        dose2_nM = dose2,
        dose1 = dose1_nM * 1e-9,
        dose2 = dose2_nM * 1e-9)

devtools::load_all()
synergy_model_v10 <- well_scores %>%
    dplyr::filter(drug_combo == "GS-441524_Nitazoxanide") %>%
    MPStats::fit_MuSyC_score_by_dose_robust(
        group_vars = vars(drug_combo),
        stan_model_args = list(verbose = TRUE),
        model_evaluation_criteria = NULL,
        iter = 16000,
        chains = 8,
        control = list(
            adapt_delta = .999,
            max_treedepth = 15),
        open_progress = FALSE,
        silent = FALSE,
        future = FALSE)

save(synergy_model_v10, file = "intermediate_data/synergy_model_v10.Rdata")
load("intermediate_data/synergy_model_v10.Rdata")

estimated_parameters <- synergy_model_v10 %>%
    dplyr::rowwise() %>%
    dplyr::do({
        combo_model <- .
        combo_model$model %>%
            tidybayes::spread_draws(
                b_logE0_Intercept,
                b_logC1_Intercept,
                b_logE1_Intercept,
                b_h1_Intercept,
                b_logC2_Intercept,
                b_logE2_Intercept,
                b_h2_Intercept,
                b_logalpha_Intercept,
                b_logE3_Intercept) %>%
            tidybayes::median_qi() %>%
            dplyr::rename_with(
                ~stringr::str_replace(., "^b_", "") %>%
                    stringr::str_replace("_Intercept", "")) %>%
            dplyr::mutate(
                drug_combo = combo_model$drug_combo)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        E0 = exp(logE0),
        C1 = exp(logC1 - 13.8),
        E1 = exp(logE1),
        s1 = h1 * (E0 + E1) / (4 * exp(logC1)),
        C2 = exp(logC2 - 13.8),
        E2 = exp(logE2),
        s2 = h2 * (E0 + E2) / (4 * exp(logC2)),
        E3 = exp(logE3),
        alpha = exp(logalpha))
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


#############################
# Analyze the fitted models #
#############################

# interactively look for problems with the model fit
synergy_model_v10$model[[5]] %>%
    shinystan::launch_shinystan()
