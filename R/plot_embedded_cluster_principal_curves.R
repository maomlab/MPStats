
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
plot_embedded_cluster_principal_curves <- function(
  cluster_principal_curves,
  thin=2000,                                         
  subtitle=NULL){

  UMAP_1_limits <- range(cluster_principal_curves$UMAP_1)
  UMAP_2_limits <- range(cluster_principal_curves$UMAP_2)

  embedding <- cluster_principal_curves %>%
    dplyr::select(UMAP_1, UMAP_2, cluster_label)

  curves <- cluster_principal_curves %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::group_modify(~if(thin && (thin < nrow(.))) dplyr::sample_n(., thin) else .) %>%
    dplyr::arrange(princurve_lambda) %>%
    dplyr::ungroup() %>%  
    dplyr::select(princurve_1, princurve_2, cluster_label) 

  # plotting with ggplot2
  plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_density_2d(
      data=embedding,
      mapping=ggplot2::aes(
        x=UMAP_1, y=UMAP_2, color=factor(cluster_label))) +
    ggplot2::geom_line(
      data=curves,
      mapping=ggplot2::aes(
        x=princurve_1, y=princurve_2, group=factor(cluster_label)),
      size=.1) +
    ggplot2::scale_x_continuous(
      "UMAP_1",
      limits=UMAP_1_limits) +
    ggplot2::scale_y_continuous(
      "UMAP_2",
      limits=UMAP_2_limits) +
    ggplot2::scale_color_discrete("Cluster Label") +
    ggplot2::ggtitle(paste0("Cluster principal curves"), subtitle=subtitle)
}
