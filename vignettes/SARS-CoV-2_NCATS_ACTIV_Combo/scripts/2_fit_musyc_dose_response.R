

library(plyr)
library(tidyverse)
library(MPStats)
library(future)
library(batchtools)

well_scores <- readr::read_tsv("intermediate_data/well_scores.Rdata")


###############################
# fit the MuSyC synergy model #
###############################

source("../../R/fit_MuSyC_score_by_dose.R")
synergy_model_c <- well_scores %>%
    fit_MuSyC_score_by_dose(
        group_vars = vars(drug_combo),
        control = list(
            adapt_delta = .99,
            max_treedepth = 12),
        stan_model_args = list(verbose = TRUE),
        model_evaluation_criteria = NULL,
        open_progress = FALSE)

synergy_model_v4.4 <- well_scores %>%
    dplyr::filter(drug_combo == "NCGC00686694-02_NCGC00388427-03") %>%
    fit_MuSyC_score_by_dose(
        group_vars = vars(drug_combo),
        E0_prior = brms::prior(student_t(200, 0, .2), nlpar = "E0", lb=0, ub=1),
        E1_prior = brms::prior(student_t(200, 0, .2), nlpar = "E1", lb=0, ub=1),
        E2_prior = brms::prior(student_t(200, 0, .2), nlpar = "E2", lb=0, ub=1),
        E3_prior = brms::prior(student_t(200, 0, .2), nlpar = "E3", lb=0, ub=1),
        E0_init = function() {as.array(brms::rstudent_t(200, 0, .2))},
        E1_init = function() {as.array(brms::rstudent_t(200, 0, .2))},
        E2_init = function() {as.array(brms::rstudent_t(200, 0, .2))},
        E3_init = function() {as.array(brms::rstudent_t(200, 0, .2))},
        control = list(
            adapt_delta = .99,
            max_treedepth = 12),
        stan_model_args = list(verbose = TRUE),
        model_evaluation_criteria = NULL,
        open_progress = FALSE)


# NCGC00686670-01_NCGC00388427-03
# 
# Population-Level Effects: 
#                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# E0_Intercept        0.09      0.00     0.09     0.09 1.00    11378    11549
# C1_Intercept        0.21      0.00     0.20     0.21 1.00     9660    10386
# E1_Intercept        0.00      0.00     0.00     0.00 1.00     9997     9891
# s1_Intercept        0.51      0.01     0.50     0.53 1.00    10231     9325
# C2_Intercept        0.52      0.02     0.51     0.57 1.00     4892     6563
# E2_Intercept        0.07      0.00     0.06     0.07 1.00    12693    11212
# s2_Intercept        2.10      1.48     0.56     6.04 1.00     4882     8430
# alpha_Intercept     0.46      0.04     0.38     0.54 1.00     9402     9478
# E3_Intercept        0.00      0.00     0.00     0.00 1.00    14399     9655

# control = NULL
# combos 5 and 6 converged
save(synergy_model, file = "intermediate_data/synergy_model.Rdata")

# control = list(adapt_delta = .99, max_treedepth = 12)
# combos 1, 2, and 5 converged
save(synergy_model_c, file = "intermediate_data/synergy_model_c.Rdata")

# narrow priors
# control = NULL
# combos 1, 2, 3, 5, 6 converged
save(synergy_model_v3, file = "intermediate_data/synergy_model_v3.Rdata")

# narrow priors
# control = list(adapt_delta = .99, max_treedepth = 12)
# combo 4 did not converge 
save(synergy_model_v4.4, file = "intermediate_data/synergy_model_v4.4.Rdata")


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
                b_alpha_Intercept,
                b_E3_Intercept) %>%
            tidybayes::median_qi() %>%
            dplyr::rename_with(
                ~stringr::str_replace(., "^b_", "") %>%
                    stringr::str_replace("_Intercept", "")) %>%
            dplyr::mutate(
                drug_combo = combo_model$drug_combo)
    }) %>%
    dplyr::ungroup()
    

estimated_parameters %>%
    readr::write_tsv("product/estimated_parameters_20201118.tsv")



#############################
# Analyze the fitted models #
#############################

# interactively look for problems with the model fit
synergy_model$model[[3]] %>%
    shinystan::launch_shinystan()




