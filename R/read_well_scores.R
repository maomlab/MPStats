#'Read well scores table
#'
#' @param file path to well scores csv, (e.g. ProbPos_CellCount_55Plates.csv)
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
read_well_scores <- function(file){
  well_scores <- readr::read_csv(
    file=file,
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
      cell_count = Children_Cells_Count_Nuclei,
      compound = Metadata_Image_Metadata_Compound,
      concentration = Metadata_Image_Metadata_Concentration,
      prob_positive = probPos) %>%
    dplyr::filter(!is.na(plate_id)) %>%
    dplyr::mutate(
      week = plate_id %>%
        stringr::str_extract("Week[0-9]+") %>%
        stringr::str_replace("Week", "") %>%
        as.integer(),
      plate_num = plate_id %>%
        stringr::str_extract("[0-9]+$") %>%
        as.integer(),
      row = well_id %>%
        stringr::str_extract("^[A-Z]") %>%
        purrr::map_int(~which(LETTERS==., arr.ind=T)),
      column = well_id %>%
        stringr::str_extract("[0-9]+$") %>%
        as.integer(),
      is_control = column %in% c(2, 11),
      log_dose = log10(concentration*1000),
      n_positive = as.integer(prob_positive * cell_count))
  
  plate_indexes <- well_scores %>%
    dplyr::distinct(week, plate_num) %>%
    dplyr::arrange(week, plate_num) %>%
    dplyr::group_by(week) %>%
    dplyr::mutate(plate = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  well_scores <- well_scores %>%
    dplyr::left_join(
      plate_indexes,
      by=c("week", "plate_num"))
  
  well_scores %>%
    dplyr::select(
      week,
      plate,
      row,
      column,
      compound,
      concentration,
      log_dose,
      is_control,
      cell_count,
      n_positive,
      prob_positive)
}