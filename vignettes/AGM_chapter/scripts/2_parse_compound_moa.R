library(plyr)
library(tidyverse)

# Load known compound mechanism of action
# and merge it with the well scores


load("intermediate_data/well_scores.Rdata")


compound_moa <- readr::read_csv(
  file="raw_data/AGM_moa.csv",
  col_types=readr::cols(
    compound = readr::col_character(),
    concentration = readr::col_double(),
    moa = readr::col_character(),
    moa2 = readr::col_character())) %>%
  dplyr::distinct(compound, moa, moa2)

unannotated <- well_scores %>%
  dplyr::distinct(compound) %>%
  dplyr::anti_join(compound_moa, by="compound")


save(compound_moa, file="intermediate_data/compound_moa_191126.Rdata")


