# parse well_scores

library(plyr)
library(tidyverse)

# load well scores

well_scores <- readr::read_csv(
  file="raw_data/ProbPos_CellCount_55Plates.csv",
  col_types=readr::cols(
    Metadata_PlateID_Nuclei = readr::col_character(),
    Metadata_WellID_Nuclei = readr::col_character(),
    Children_Cells_Count_Nuclei = readr::col_double(),
    Metadata_Image_Metadata_Compound = readr::col_character(),
    Metadata_Image_Metadata_Concentration = readr::col_double(),
    probPos = readr::col_double(),
    probNeg = readr::col_double())) %>%
  dplyr::transmute(
    plate_id = Metadata_PlateID_Nuclei,
    well_id = Metadata_WellID_Nuclei,
    well_cell_count = Children_Cells_Count_Nuclei,
    compound = Metadata_Image_Metadata_Compound,
    concentration = Metadata_Image_Metadata_Concentration,
    prob_positive = probPos,
    prob_negative = probNeg) %>%
  dplyr::filter(!is.na(plate_id)) %>%
  dplyr::mutate(
    week_id = plate_id %>%
      stringr::str_extract("Week[0-9]+") %>%
      stringr::str_replace("Week", "") %>%
      as.integer(),
    plate_num = plate_id %>%
      stringr::str_extract("[0-9]+$") %>%
      as.integer(),
    column = well_id %>%
      stringr::str_extract("^[A-Z]") %>%
      purrr::map_int(~which(LETTERS==., arr.ind=T)),
    row = well_id %>%
      stringr::str_extract("[0-9]+$") %>%
      as.integer(),
    control = row %in% c(2, 11))

plate_indexes <- well_scores %>%
  dplyr::distinct(week_id, plate_num) %>%
  dplyr::arrange(week_id, plate_num) %>%
  dplyr::group_by(week_id) %>%
  dplyr::mutate(plate_ind = dplyr::row_number()) %>%
  dplyr::ungroup()

well_scores <- well_scores %>%
  dplyr::left_join(
    plate_indexes,
    by=c("week_id", "plate_num"))

### example plate ####
well_scores %>%
  dplyr::filter(week_id == 5, plate_ind==1) %>%
  dplyr::mutate(well_name = paste0(compound, "|", concentration)) %>%
  dplyr::arrange(row, column) %>%
  tidyr::pivot_wider(
    id_cols="row",
    names_from="column",
    values_from="well_name") %>%
  readr::write_tsv("product/example_plate_map_191126.tsv")
  

# example plate map for plate_id: Week10_40111
#   row `2`                             `3`        `4`       `5`                `6`       `7`           
#  <dbl> <chr>                           <chr>      <chr>     <chr>              <chr>     <chr>         
#1     2 DMSO|0                          DMSO|0     DMSO|0    taxol|0.3          taxol|0.3 taxol|0.3     
#2     3 Cdk1/2 inhibitor (NU6102)|10    AZ138|30   AZ-U|30   temozolomide|20    TKK|10    monastrol|100 
#3     4 Cdk1/2 inhibitor (NU6102)|3     AZ138|10   AZ-U|10   temozolomide|6     TKK|3     monastrol|30  
#4     5 Cdk1/2 inhibitor (NU6102)|1     AZ138|3    AZ-U|3    temozolomide|2     TKK|1     monastrol|10  
#5     6 Cdk1/2 inhibitor (NU6102)|0.3   AZ138|1    AZ-U|1    temozolomide|0.6   TKK|0.3   monastrol|3   
#6     7 Cdk1/2 inhibitor (NU6102)|0.1   AZ138|0.3  AZ-U|0.3  temozolomide|0.2   TKK|0.1   monastrol|1   
#7     8 Cdk1/2 inhibitor (NU6102)|0.03  AZ138|0.1  AZ-U|0.1  temozolomide|0.06  TKK|0.03  monastrol|0.3 
#8     9 Cdk1/2 inhibitor (NU6102)|0.01  AZ138|0.03 AZ-U|0.03 temozolomide|0.02  TKK|0.01  monastrol|0.1 
#9    10 Cdk1/2 inhibitor (NU6102)|0.003 AZ138|0.01 AZ-U|0.01 temozolomide|0.006 TKK|0.003 monastrol|0.03
#10   11 taxol|0.3                       taxol|0.3  taxol|0.3 DMSO|0             DMSO|0    DMSO|0  



save(well_scores, file="intermediate_data/well_scores.Rdata")
