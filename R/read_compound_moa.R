
#'Read well scores table
#'
#' @param filename path to compound mechanism of action csv, (e.g. AGM_moa.csv)
#' @return compound_moa tibble::tibble with columns
#'   Classes 'spec_tbl_df', 'tbl_df', 'tbl' and 'data.frame':	40 obs. of  3 variables:
#'     $ compound      : chr  "PP-2" "emetine" "AZ258" "cytochalasin B" ...
#'     $ concentration : num  3 0.3 1 10 3 0.01 0.3 3 0.3 0.3 ...
#'     $ moa           : chr  "Epithelial" "Protein synthesis" "Aurora kinase inhibitors" "Actin disruptors" ...
#'     $ moa2          : chr  "other" "other" "other" "cytoskeleton" ...
#'
#' @export
read_compound_moa <- function(file){
  compound_moa <- readr::read_csv(
    file=file,
    col_types=readr::cols(
      compound = readr::col_character(),
      concentration = readr::col_double(),
      moa = readr::col_character(),
      moa2 = readr::col_character())) %>%
    dplyr::distinct(compound, moa, moa2, .keep_all=TRUE)
}