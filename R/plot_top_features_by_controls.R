

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#' Plot top features by controls as a split-violin plot
#'
#' For each object, plot the top 4 feature by Zprime colored by the
#' dye used to identify it
#'
#' @export
plot_Zprime_by_top_features <- function(cell_features, top_k_features=4){
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
               !(column_name %in% c("COND", "Group_Index","ObjectNumber:Nuclei")),
               !stringr::str_detect(column_name, "^Location"))

    cell_features <- cell_features %>%
      dplyr::rename(
        plate_id = `Metadata_PlateID:Nuclei`,
        well_id = `Metadata_WellID:Nuclei`,                 
        condition = COND) %>%
      dplyr::filter(condition %in% c("PC", "NC")) %>%        
      dplyr::mutate(condition=ifelse(condition=="PC", "Positive", "Negative"))  
  
    # compute over plate and then average
    Zprime_scores <- cell_features %>%
      plyr::ddply("plate_id", function(plate_features){
        cat("Computing Zprime scores for plate: ", plate_features$plate_id[1], "\n", sep="")
        plyr::ldply(cell_feature_columns$column_name, function(feature_id){            
          positives <- plate_features %>%
            dplyr::filter(condition == "Positive") %>%
            magrittr::extract2(feature_id)
          negatives <- plate_features %>%
            dplyr::filter(condition == "Negative") %>%
            magrittr::extract2(feature_id)
          data.frame(
            feature_id=feature_id,
            Zprime=Zprime(positives, negatives))
        })
      })

    top_features <- Zprime_scores %>%
      dplyr::group_by(feature_id) %>%
      dplyr::summarize(
        Zprime_mean=mean(Zprime),
        Zprime_std_err=sd(Zprime)/sqrt(dplyr::n())) %>%
      tidyr::separate(
        col=feature_id,
        into=c("feature_name", "object", "dye", "feature_type"),
        sep=":",
        remove=FALSE) %>%
      dplyr::filter(!is.na(Zprime_mean)) %>%
      dplyr::group_by(object) %>%
      dplyr::arrange(dplyr::desc(Zprime_mean)) %>%
      dplyr::slice(1:top_k_features) %>%
      dplyr::ungroup()

    top_cell_features <- cell_features %>%
      dplyr::select(which(colnames(cell_features) %in% c(as.character(top_features$feature_id), "condition"))) %>%
      tidyr::pivot_longer(
        cols=-condition,
        names_to="feature_id") %>%
      dplyr::group_by(feature_id) %>%
        dplyr::mutate(value = robustHD::standardize(value)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(y=top_features, by="feature_id") %>%
      dplyr::mutate(
        feature_name=feature_name %>% stringr::str_replace("_", "\n"))

#    top_well_features <- cell_features %>%
#      dplyr::count(plate_id, well_id, condition, name="value") %>%
#      dplyr::mutate(
#        feature_id = "well_cell_count",
#        feature_name = "Well Cell Count",
#        object="Wells",
#        feature_type="Count") %>%
#      dplyr::group_by(feature_id) %>%
#      dplyr::mutate(
#        value=robustHD::standardize(value))
  
    plot <- ggplot2::ggplot(data=top_cell_features) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position=c(-.6,0.05)) +
      ggplot2::geom_hline(yintercept=0, size=.3) +
      geom_split_violin(
        mapping=ggplot2::aes(
          x=reorder(feature_name, Zprime_mean),
          y=value,
          fill=condition),
        width=1.5) +
      ggplot2::facet_wrap(                 
        facets=dplyr::vars(object),
        scales="free_y",
        ncol=1,
        strip.position="right") +
      ggplot2::coord_flip() +
      ggplot2::scale_x_discrete("Feature Type") +
      ggplot2::scale_y_continuous(
        "Standardized Feature Value",
        breaks=c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
        limits=c(-5, 5),
        expand=c(0,0)) +
      ggplot2::scale_fill_manual(
        "Control",
        values=c(
            Positive="#224EDE",
            Negative="#DEB23E")) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
}    
