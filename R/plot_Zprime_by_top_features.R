
#' Compute Zprime score
#'
#' J. H. Zhang, T. D. Chung, K. R. Oldenburg. A Simple Statistical Parameter for Use in Evaluation and Validation of High Throughput Screening Assays. J Biomol Screening, 1999.
#'
#' @export
Zprime <- function(positives, negatives) {
    1 - 3 * (sd(positives) + sd(negatives))/abs(mean(positives) - mean(negatives))
}


#' Plot Zprime by top features
#'
#' For each object, plot the top 5 feature by Zprime colored by the
#' dye used to identify it
#'
#' @export
plot_Zprime_by_top_features <- function(
  cell_features,
  top_k_features=4,
  plate_id_column=`Metadata_PlateID:Nuclei`,
  condition_column=COND,
  group_column=dye,
  group_label=deparse(substitute(group_column)),
  group_colors=c(
        Hoechst="#4472C4",  # blue
        Tubulin="#00B050",  # green
        Actin="#ED3833")    # red
){

  cell_feature_columns <- cell_features %>%
    colnames() %>%
    tibble::tibble(column_name=.) %>%
    dplyr::mutate(col_number = dplyr::row_number()) %>%
      dplyr::filter(
               !stringr::str_detect(column_name, "^Metadata"),
               !stringr::str_detect(column_name, "^ImageNumber"),
               !stringr::str_detect(column_name, "^Parent"),
               !stringr::str_detect(column_name, "^Children"),
               !stringr::str_detect(column_name, "^Number"),
               !(column_name %in% c(
                   deparse(substitute(condition_column)),
                   "Group_Index",
                   "ObjectNumber:Nuclei")),
               !stringr::str_detect(column_name, "^Location"))

    # compute over plate and then average
    Zprime_scores <- cell_features %>%
      dplyr::filter(!!enquo(condition_column) %in% c("PC", "NC")) %>%
      plyr::ddply(plate_id_column, function(plate_features){
        plate_id <- plate_features[1, plate_id_column %>% substitute %>% deparse]
        cat("Computing Zprime scores for plate: ", plate_id, "\n", sep="")
        plyr::ldply(cell_feature_columns$column_name, function(feature_id){   
          positives <- plate_features %>%
            dplyr::filter(!!enquo(condition_column) == "PC") %>%
            magrittr::extract2(feature_id)
          negatives <- plate_features %>%
            dplyr::filter(!!enquo(condition_column) == "NC") %>%
            magrittr::extract2(feature_id)
          data.frame(
            feature_id=feature_id,
            Zprime=Zprime(positives, negatives))
        })
      })

    Zprime_summary <- Zprime_scores %>%
      dplyr::group_by(feature_id) %>%
      dplyr::summarize(
        Zprime_mean=mean(Zprime),
        Zprime_std_err=sd(Zprime)/sqrt(dplyr::n()))

    Zprime_summary <- Zprime_summary %>%
      tidyr::separate(
        feature_id,
        into=c("feature_name", "object", "dye", "feature_type"),
        sep=":") %>%
      dplyr::mutate(feature_name = feature_name %>% stringr::str_replace("_", "\n"))

    plot <- ggplot2::ggplot(
      data=Zprime_summary %>%
        dplyr::filter(!is.na(Zprime_mean)) %>%
        dplyr::group_by(object) %>%
        dplyr::arrange(dplyr::desc(Zprime_mean)) %>%
        dplyr::slice(1:top_k_features) %>%
        dplyr::ungroup()) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position=c(-0.95,0.1)) +  
      ggplot2::geom_rect(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf, fill="gray85") +
      ggplot2::geom_rect(xmin=0, xmax=0.5, ymin=-Inf, ymax=Inf, fill="gray90") +
      ggplot2::geom_rect(xmin=0.5, xmax=1, ymin=-Inf, ymax=Inf, fill="gray95") +
      ggplot2::geom_hline(yintercept=c(1:top_k_features), size=.5, color="grey92") +
      ggplot2::geom_vline(xintercept=c(-6, -4, -2, 0, 0.5, 1), size=.5, color="grey92") +
      ggplot2::geom_errorbarh(
        mapping=ggplot2::aes(
          y=reorder(feature_name, Zprime_mean),
          xmin=Zprime_mean - Zprime_std_err,
          xmax=Zprime_mean + Zprime_std_err,
          color=!!enquo(group_column))) +
      ggplot2::geom_point(
        mapping=ggplot2::aes(
          y=reorder(feature_name, Zprime_mean),
          x=Zprime_mean,
          color=!!enquo(group_column)),
        size=3) +
      ggplot2::facet_wrap(                 
        facets=dplyr::vars(object),
        scales="free_y",
        ncol=1,
        strip.position="right") +
      ggplot2::scale_y_discrete("Feature Type") +
      ggplot2::scale_x_continuous(
        "Zprime Score",
        limits=c(-6, 1),
        breaks=c(-6, -4, -2, 0, 1),
        expand=c(0,0))  +
      ggplot2::scale_color_manual(
        group_label,
        values=group_colors)+
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
}    
