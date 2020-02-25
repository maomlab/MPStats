

#'Read Zprime by plate
#'
#' @param input path to Zprime csv (e.g. Master55.csv) or
#'              data.frame with columns
#'                  Metadata_PlateID_Nuclei               : chr
#'                  Metadata_WellID_Nuclei                : chr
#'                  Children_Cells_Count_Nuclei           : num
#'                  Metadata_Image_Metadata_Compound      : chr
#'                  Metadata_Image_Metadata_Concentration : num
#'                  probPos                               : num
#'                  probNeg                               : num
#'                      
#' @return well_scores data.frame with columns
#'   week_id    : int  10 10 10 10 10 10 10 10 10 10 ...
#'   plate_ind  : int  1 1 1 1 1 1 1 1 1 1 ...
#'   row        : int  2 2 2 2 2 2 2 2 2 2 ...
#'   column     : int  2 3 4 5 6 7 8 9 10 11 ...
#'   compound   : chr  "DMSO" "Cdk1/2 inhibitor (NU6102)" "Cdk1/2 inhibitor (NU6102)" ...
#'   log_dose   : num  -Inf 4 3.48 3 2.48 ...
#'   is_control : logi  TRUE FALSE FALSE FALSE FALSE FALSE ...
#'   cell_count : num  1442 479 743 888 951 ...
#'   n_positive : int  0 0 0 0 0 0 0 0 0 303 ...
#'  
#' @export
read_well_scores <- function(input){
  Zprime_by_plate <- readr::read_csv(
    input,
    col_types=readr::cols(
      `Metadata_PlateID:Nuclei` = readr::col_character(),
      `Positive Control` = readr::col_character(),
      `Negative Control` = readr::col_character())) %>%
    dplyr::mutate(plate_id = `Metadata_PlateID:Nuclei`) %>%
    dplyr::select(-`Metadata_PlateID:Nuclei`, -`Positive Control`, -`Negative Control`) %>%
    tidyr::pivot_longer(
      cols=dplyr::starts_with("Multivariate Z Prime:"),
      names_to="dye_set",
      values_to="Zprime") %>%
    dplyr::mutate(dye_set = dye_set %>% stringr::str_replace("Multivariate Z Prime:", "")) %>%
    dplyr::filter(plate_id != "Mean")
}
  