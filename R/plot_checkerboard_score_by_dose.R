

#' Plot dose-response checkerboard
#'
#' @param well_scores data.frame with columns [dose_1, dose_2, score]
#' @param contour_color
#'
#' @export
plot_checkerboard_score_by_dose <- function(
   well_scores,
   contour_color = "gold"){

  ggplot2::ggplot(data = well_scores) + 
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      panel.background=element_rect(fill="grey40", colour="grey40")) +
    ggplot2::geom_tile(
      mapping = ggplot2::aes(
        x = log10(dose_1),
        y = log10(dose_2),
        fill = score)) +
    ggplot2::geom_contour(
      mapping = ggplot2::aes(
        x = log10(dose_1),
        y = log10(dose_2),
        z = score),
      color = contour_color) +
    ggplot2::coord_fixed()    
}                                            
