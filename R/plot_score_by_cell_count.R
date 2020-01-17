#' plot score by cell count
#' 
#' @param well_scores tibble::tibble as output by read_well_scores
#' @param subtitle plot subtitle, typically the study identifier
#' @return ggplot2 object 
#' 
#' @export
plot_score_by_cell_count <- function(well_scores, subtitle=NULL){
  p <- ggplot2::ggplot(data=well_scores) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
      mapping=ggplot2::aes(
        x=sqrt(cell_count),
        y=prob_positive)) +
    ggplot2::geom_smooth(
      mapping=ggplot2::aes(
        x=sqrt(cell_count),
        y=prob_positive)) +
    ggplot2::ggtitle(
      label="Classifier Score by Well Cell Count",
      subtitle=subtitle) +
    ggplot2::scale_x_continuous(
      "Well Cell Count",
      breaks=c(0,10,20,30,40),
      labels=c(0,100,400,900,1600)) +
    ggplot2::scale_y_continuous(
      "Score",
      limits=c(0,1),
      labels=scales::percent_format())
}