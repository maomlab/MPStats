
#'Read compound mechanism of action table
#'
#' @param input path to compound mechanism of action csv, (e.g. AGM_moa.csv) or
#'              data.frame with columns
#'                  compound      : chr    # name of compound
#'                  concentration : num    # 
#'                  moa           : chr    # detailed mechanism of action
#'                  moa2          : chr    # class of mechanism of action
#'       
#' @return compound_moa tibble::tibble with columns
#'   Classes 'spec_tbl_df', 'tbl_df', 'tbl' and 'data.frame':	40 obs. of  3 variables:
#'     $ compound      : chr  "PP-2" "emetine" "AZ258" "cytochalasin B" ...
#'     $ concentration : num  3 0.3 1 10 3 0.01 0.3 3 0.3 0.3 ...
#'     $ moa           : chr  "Epithelial" "Protein synthesis" "Aurora kinase inhibitors" "Actin disruptors" ...
#'     $ moa2          : chr  "other" "other" "other" "cytoskeleton" ...
#'
#' @export
read_compound_moa <- function(input){
  if(is.character(input)){
    cat("Reading compound MOA table from '", input, "' ...\n", sep="")
    compound_moa <- readr::read_csv(
      file=input,
      col_types=readr::cols(
        compound = readr::col_character(),
        concentration = readr::col_double(),
        moa = readr::col_character(),
        moa2 = readr::col_character()))
  } else if(is.data.frame(input)){
    cat("Parsing compound MOA table ...\n")
    compound_moa <- input
  } else {
    cat("Unable to read compound MOA table\n")
    stop()
  }
  
  compound_moa <- compound_moa %>%
    dplyr::distinct(compound, moa, moa2, .keep_all=TRUE)
}