
#' Model cell count by batch variables using a linear model
#' 
#' Here using a sqrt() transform to 
#' 
#' @param well_scores result of read_well_scores()
#' @return "lm" object of model
#' 
#' @export
model_cell_count_by_batch_vars_lm <- function(well_scores){
  model <- lm(
    sqrt(cell_count) ~ week +  plate + row,
    data=well_scores %>% dplyr::filter(!is_control))
}

