
#' plot cell count by batch variables as a scatterplot
#' 
#' For each batch variable, show a facet panel with sqrt(cell_count) by batch variable
#'
#' 
#'@param well_scores tibble::tibble as output by read_well_scores
#'@param subtitle plot subtitle, typically the study identifier
#'@return ggplot2 object with the lattice plot with a panel for each batch variable
#' 
#'@export 
plot_cell_count_by_batch_vars_scatter <- function(well_scores, subtitle=NULL){
  data <- well_scores %>%
    dplyr::filter(!is_control) %>%
    dplyr::select(
      week, plate, row,           # batch variables
      cell_count) %>%        # response value
    tidyr::pivot_longer(
      cols=c("week", "plate", "row"),
      names_to="batch_variable",
      values_to="batch_value")
  
  p <- ggplot2::ggplot(data=data) +
    ggplot2::theme_bw() +
    ggplot2::geom_jitter(
      mapping=ggplot2::aes(
        x=batch_value,
        y=sqrt(cell_count)),
      size=.8,
      height=0,
      width=.2) +
    ggplot2::geom_smooth(
      mapping=ggplot2::aes(
        x=batch_value,
        y=sqrt(cell_count)),
      method="lm") +
    ggplot2::facet_wrap(
      facets=ggplot2::vars(batch_variable),
      scales="free_x") +
    ggplot2::scale_x_continuous(
      "Batch Value") +
    ggplot2::scale_y_continuous(
      "Well Cell Count",
      breaks=c(0,10,20,30,40),
      labels=c(0,100,400,900,1600)) +
    ggplot2::ggtitle(
      label="Well Cell Count by Batch Variables",
      subtitle=subtitle)
}