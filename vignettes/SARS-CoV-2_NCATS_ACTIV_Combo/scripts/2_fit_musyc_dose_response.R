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

well_scores <- readr::read_tsv("intermediate_data/well_scores.tsv")

synergy_model_v5 <- well_scores %>%
    MPStats::fit_MuSyC_score_by_dose(
        group_vars = vars(drug_combo),
        E0_prior = brms::prior(student_t(200, 0, .2), nlpar = "E0", lb=0, ub=1),
        E1_prior = brms::prior(student_t(200, 0, .2), nlpar = "E1", lb=0, ub=1),
        E2_prior = brms::prior(student_t(200, 0, .2), nlpar = "E2", lb=0, ub=1),
        E3_alpha_prior = brms::prior(student_t(200, 0, .4), nlpar = "E3alpha", lb=0, ub=3),
        E0_init = function() {as.array(brms::rstudent_t(1, 200, 0, .2))},
        E1_init = function() {as.array(brms::rstudent_t(1, 200, 0, .2))},
        E2_init = function() {as.array(brms::rstudent_t(1, 200, 0, .2))},
        E3_alpha_init = function() {as.array(brms::rstudent_t(1, 200, 0, .))},
        chains = 5,
        #control = list(
        #    adapt_delta = .99,
        #    max_treedepth = 12),
        stan_model_args = list(verbose = TRUE),
        model_evaluation_criteria = NULL,
        open_progress = FALSE,
        seed = 1234,
        silent = FALSE,
        future = TRUE)


save(synergy_model_v5, file = "intermediate_data/synergy_model_v5.Rdata")

estimated_parameters <- synergy_model_v3 %>%
    dplyr::rowwise() %>%
    dplyr::do({
        combo_model <- .
        combo_model$model %>%
            tidybayes::spread_draws(
                b_E0_Intercept,
                b_C1_Intercept,
                b_E1_Intercept,
                b_s1_Intercept,
                b_C2_Intercept,
                b_E2_Intercept,
                b_s2_Intercept,
                b_log10alpha_Intercept,
                b_E3_alpha_Intercept) %>%
            tidybayes::median_qi() %>%
            dplyr::rename_with(
                ~stringr::str_replace(., "^b_", "") %>%
                    stringr::str_replace("_Intercept", "")) %>%
            dplyr::mutate(
                drug_combo = combo_model$drug_combo)
    }) %>%
    dplyr::ungroup()


estimated_parameters %>%
    readr::write_tsv("product/estimated_parameters_20201121.tsv")


#############################
# Analyze the fitted models #
#############################

# interactively look for problems with the model fit
synergy_model$model[[3]] %>%
    shinystan::launch_shinystan()
