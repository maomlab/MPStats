#' plot cell count by batch variables as a scatter plot
#' 
#'@param well_scores tibble::tibble as output by read_well_scores
#'@param subtitle plot subtitle, typically the study identifier
#'@return ggplot2 object with the lattice plot with a panel for each batch variable
#' 
#'@export 
plot_score_by_batch_vars_scatter <- function(well_scores, subtitle=NULL){
  data <- well_scores %>%
    dplyr::filter(!is_control) %>%
    dplyr::select(
      week, plate, row,                  # batch variables
      prob_positive)%>%                  # response values
    tidyr::pivot_longer(
      cols=c("week", "plate", "row"),
      names_to="batch_variable",
      values_to="batch_value")
  
  p <- ggplot2::ggplot(data=data) +
    ggplot2::theme_bw() +
    ggplot2::geom_jitter(
      mapping=ggplot2::aes(
        x=batch_value,
        y=prob_positive),
      size=.8,
      height=0,
      width=.2) +
    ggplot2::geom_smooth(
      mapping=ggplot2::aes(
        x=batch_value,
        y=prob_positive),
      method="lm") +
    ggplot2::facet_wrap(
      facets=ggplot2::vars(batch_variable),
      scales="free_x") +
    ggplot2::scale_x_continuous(
      "Batch Value") +
    ggplot2::scale_y_continuous(
      "Score",
      limits=c(0,1),
      labels=scales::percent_format()) +
    ggplot2::ggtitle(
      label="Classifier Score by Batch Dimensions",
      subtitle=subtitle)
}