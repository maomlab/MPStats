#'Read cell features
#'
#' @param input path to cell features csv (e.g. AGM_FeatureSet_RAW55.csv) or
#'                      
#' @return cell_features tibble::tibble with columns
#'     COND
#'     Texture*
#'     RadialDistribution*
#'     Intensity*
#'     Granularity*
#'     Correlation*
#'     AreaShape*

#' @export
read_cell_features <- function(input){
  if(is.character(input)){
    cell_features <- readr::read_csv(file=input)
  }
}

                                       
