library(plyr)
library(tidyverse)


plate_map <- readr::read_tsv("intermediate_data/plate_map.tsv")

field_scores <- readr::read_csv("raw_data/NCATS_combos_20201115.csv") %>%
    dplyr::select(
        plate_id = PlateID,
        well_id = Image_Metadata_WellID,
        field_id = Image_Metadata_FieldID,
        nuclei_count = Image_Count_Nuclei,
        syn_nucs_count = Image_Count_syn_nucs) %>%
    dplyr::left_join(
        plate_map,
        by = c("plate_id", "well_id"))

field_scores %>%
    readr::write_tsv("intermediate_data/field_scores.tsv")
