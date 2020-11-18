

library(plyr)
library(tidyverse)
library(MPStats)
library(future)
library(batchtools)

reg = makeRegistry(NA)





field_scores <- readr::read_tsv("intermediate_data/field_scores.tsv")

well_scores <- field_scores %>%
    # assign zero and single agent treatments to the combo treatmetns
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
         drug2_treatments <-field_scores %>%
                dplyr::filter(
                    SampleID1 == data$SampleID1[1],
                    SampleID2 == "DMSO") %>%
                dplyr::mutate(
                    SampleID2_stock_conc_in_uM = 0,
                    drug_combo = data$drug_combo[1])
        combo_treatments <- data
        cat(data$drug_combo[1], "\n")
        cat("  n DMSO treatments: ", nrow(DMSO_treatments)/16, "\n", sep = "")
        cat("  n drug1 treatments: ", nrow(drug1_treatments)/16, "\n", sep = "")
        cat("  n drug2 treatments: ", nrow(drug2_treatments)/16, "\n", sep = "")
        cat("  n combo treatments: ", nrow(combo_treatments)/16, "\n", sep = "")        
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



source("../../R/fit_MuSyC_score_by_dose.R")

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




synergy_model_c <- well_scores %>%
    fit_MuSyC_score_by_dose(
        group_vars = vars(drug_combo),
        control = list(
            adapt_delta = .99,
            max_treedepth = 12),
        stan_model_args = list(verbose = TRUE),
        model_evaluation_criteria = NULL,
        open_progress = FALSE)


