
#' plot feature correlation matrix
#'
#' @param cell_features tibble::tibble with a column for each feature and a row for each cell
#' @param subtitle plot subtitle, typically the study identifier
#' @return ggplot2 object with a matrix where cell (i,j) is colored by the correlation of
#'         feature i and feature j
#'
#'@export
