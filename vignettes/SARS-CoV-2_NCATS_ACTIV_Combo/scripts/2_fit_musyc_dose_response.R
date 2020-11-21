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


field_scores <- readr::read_tsv("intermediate_data/field_scores.tsv")


##########################################
# Aggregate field data to the well level #
##########################################
# one trick thing is that DMSO controls and
# single drug data can be re-used for different
# combo screens.

well_scores <- field_scores %>%
    # assign zero and single agent treatments to the combo treatments
    dplyr::filter(SampleID1 != "DMSO", SampleID2 != "DMSO") %>%
    dplyr::group_by(SampleID1, SampleID2) %>%
    dplyr::do({
        data <- .
        DMSO_treatments <- field_scores %>%
            dplyr::filter(SampleID1 == "DMSO", SampleID2 == "DMSO") %>%
            dplyr::mutate(
                SampleID1_stock_conc_in_uM = 0,
                SampleID2_stock_conc_in_uM = 0,
                drug_combo = data$drug_combo[1])
         drug1_treatments <- field_scores %>%
                dplyr::filter(
                    SampleID1 == "DMSO",
                    SampleID2 == data$SampleID2[1]) %>%
                dplyr::mutate(
                    SampleID1_stock_conc_in_uM = 0,
                    drug_combo = data$drug_combo[1])
         drug2_treatments <- field_scores %>%
                dplyr::filter(
                    SampleID1 == data$SampleID1[1],
                    SampleID2 == "DMSO") %>%
                dplyr::mutate(
                    SampleID2_stock_conc_in_uM = 0,
                    drug_combo = data$drug_combo[1])
        combo_treatments <- data
        cat(data$drug_combo[1], "\n")
        cat("  n DMSO treatments: ", nrow(DMSO_treatments) / 16, "\n", sep = "")
        cat("  n drug1 treatments: ", nrow(drug1_treatments) / 16, "\n", sep = "")
        cat("  n drug2 treatments: ", nrow(drug2_treatments) / 16, "\n", sep = "")
        cat("  n combo treatments: ", nrow(combo_treatments) / 16, "\n", sep = "")
        dplyr::bind_rows(
            DMSO_treatments,
            drug1_treatments,
            drug2_treatments,
            combo_treatments)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
        dose1 = SampleID1_stock_conc_in_uM,
        dose2 = SampleID2_stock_conc_in_uM) %>%
    dplyr::group_by(plate_id, drug_combo, well_id, dose1, dose2) %>%
    dplyr::summarize(
        n_positive = sum(syn_nucs_count),
        count = sum(nuclei_count),
        .groups = "drop")

#check per combo number of wells
#   combo: (plate replicates) * (drug1 doses) * (drug2 doses) +
#   drug1: (plate_replicates) * (drug1_doses) * 3
#   drug2: (plate_replicates) * (drug2_doses) * 2
#   DMSO:  (plate_replicates) * (n combos)
#
#   5*5*5 + 5*5*3 + 5*5*2 + 5*6 = 280
well_scores %>% dplyr::count(drug_combo)

###############################################
# plot checkerboards for all the combinations #
###############################################
treatment_scores <- well_scores %>%
    dplyr::group_by(drug_combo, dose1, dose2) %>%
    dplyr::summarize(
        score = median(n_positive / count),
        .groups = "drop")

treatment_scores %>%
    dplyr::group_by(drug_combo) %>%
    dplyr::do({
        data <- .
        drug_combo <- data$drug_combo[1]
        drug1 <- drug_combo %>% stringr::str_replace("_.+$", "")
        drug2 <- drug_combo %>% stringr::str_replace("^.+_", "")
        plot <- MPStats::plot_checkerboard_score_by_dose(
            treatment_scores = data,
            treatment_1_label = drug1,
            treatment_2_label = drug2,
            treatment_1_units = "nM",
            treatment_2_units = "nM") +
            ggplot2::scale_fill_continuous(
                "% Infected",
                limits = c(0, .1),
                breaks = c(.0, .05, .1),
                labels = scales::percent(c(0, .05, .1)))
        ggplot2::ggsave(
            paste0("product/figures/checkerboard_", data$drug_combo[1], ".pdf"),
            width = 6,
            height = 6)
        data.frame()
        })



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

synergy_model_v2 <- well_scores %>%
    dplyr::filter(drug_combo == "NCGC00686694-02_NCGC00090774-08") %>%
    MPStats::fit_MuSyC_score_by_dose(
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
        chains = 50,
        iter = 2000,
        stan_model_args = list(verbose = TRUE),
        model_evaluation_criteria = NULL,
        open_progress = FALSE,
        silent = FALSE,
        future = TRUE)


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



estimated_parameters <- synergy_model_c %>%
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
    readr::write_tsv("product/estiamted_parameters_20201118.tsv")



#############################
# Analyze the fitted modesl #
#############################

# interactively look for problems with the model fit
synergy_model$model[[3]] %>%
    shinystan::launch_shinystan()



source("../../R/fit_MuSyC_score_by_dose.R")
treatment_scores <- well_scores %>%
    dplyr::group_by(drug_combo, dose1, dose2) %>%
    dplyr::summarize(
        score = median(n_positive / count),
        .groups = "drop") %>%
    dplyr::group_by(drug_combo) %>%
    dplyr::mutate(
        d1_scale_factor = max(dose1),
        d2_scale_factor = max(dose2)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
        estimated_parameters,
        by = c("drug_combo")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        prior_score = generate_MuSyC_effects(
            d1 = dose1 / d1_scale_factor,
            d2 = dose2 / d2_scale_factor,
            E0 = E0,
            C1 = 0.5, E1 = E1, s1 = 0.95,
            C2 = 0.5, E2 = E2, s2 = 0.95,
            alpha = 1.5, E3 = E3),
        fitted_score = generate_MuSyC_effects(
            d1 = dose1 / d1_scale_factor,
            d2 = dose2 / d2_scale_factor,
            E0 = E0,
            C1 = C1, E1 = E1, s1 = s1,
            C2 = C2, E2 = E2, s2 = s2,
            alpha = alpha, E3 = E3)) %>%
    dplyr::ungroup()

source("../../R/plot_checkerboard_score_by_dose.R")
treatment_scores %>%
    dplyr::group_by(drug_combo) %>%
    dplyr::do({
        data <- .
        drug_combo <- data$drug_combo[1]
        drug1 <- drug_combo %>% stringr::str_replace("_.+$", "")
        drug2 <- drug_combo %>% stringr::str_replace("^.+_", "")
        plot <- plot_checkerboard_score_by_dose(
            treatment_scores = data,
            treatment_1_label = drug1,
            treatment_2_label = drug2,
            treatment_1_units = "nM",
            treatment_2_units = "nM") +
            ggplot2::geom_contour(
                mapping = ggplot2::aes(
                    x = log10(dose1),
                    y = log10(dose2),
                    z = prior_score),
                color = "green") +
            ggplot2::geom_contour(
                mapping = ggplot2::aes(
                    x = log10(dose1),
                    y = log10(dose2),
                    z = fitted_score),
                color = "purple") +
            ggplot2::scale_fill_continuous(
                "% Infected",
                limits = c(0, .1),
                breaks = c(.0, .05, .1),
                labels = scales::percent(c(0, .05, .1)))
        ggplot2::ggsave(
            paste0("product/figures/checkerboard_fitted_", data$drug_combo[1], ".pdf"),
            width = 6,
            height = 6)
        data.frame()
        })
