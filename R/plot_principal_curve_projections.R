
#' Plot embedded cluster principal cluves
#'
#'    This makes a UMAP density plot and adds principal curves to each cluster
#' 
#' @param cluster_curves  tibble::tibble with columns
#'        UMAP_1  :         coordinate 1 of the umap embedding
#'        UMAP_2  :         coordiante 2 of the umap embedding
#'        cluster_id :      cluster identifier
#'        princurve_1 :     coordinate 1 of the principal curve projection
#'        princurve_2 :     coordiante 2 of the principal curve projection
#'        princurve_order : order of points along the principal curve projection
#' @param subtitle plot subtitle
#'
#' @usage
#'
#'  cluster_principal_curves <- dplyr::bind_cols(
#'    arrow::read_parquet("input/umap_embedding.parquet"),
#'    arrow::read_parquet("input/hbscan_clustering.parquet")) %>%
#'    MPStats::fit_principal_curves() %>%
#'    dplyr::arrange(princurve_order)
#' 
#' plot <- MPStats::plot_embedded_principal_curves(cluster_principal_curves)
#' ggplot2::ggsave(
#'   plot=plot,
#'   filename=paste0("product/embedded_principal_curves_", MPStats::date_code(), ".pdf"),
#'   width=6, height=6,
#'   useDingbats=FALSE)
#' ggplot2::ggsave(
#'   plot=dye_plot,
#'   filename=paste0("product/embedded_principal_curves_", MPStats::date_code(), ".png"),
#'   width=6, height=6)
#'
#' @export
plot_principal_curve_projections <- function(
  cluster_principal_curves,
  subtitle=NULL){

    plot <- ggplot2::ggplot() +
      ggplot2::theme_bw() +
      geom_split_violin(
        data=cluster_principal_curves,
        mapping=ggplot2::aes(
          x=factor(cluster_id),
          y=lambda,
          fill=condition)) +
      ggplot2::coord_flip() +
      ggplot2::scale_x_discrete("Cluster") +
      ggplot2::scale_y_continuous("Principal Curve") +
      ggplot2::scale_fill_manual(
        "Control",
        values=c(
            Positive="#224EDE",
            Negative="#DEB23E")) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::ggtitle("Principal Curve Projection Density", subtitle=subtitle)
}
