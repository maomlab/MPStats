#' Plot Zprime score by Dye combination
#'
#'    This makes a plot like from the UpSetR package with two sections
#'    The upper section is a dot plot with the Zprime score on the
#'    y-axis and the dye set on the x-axis. The value for each dye set
#'    is represented by a dot at the mean for each plate and error
#'    bars as the standard error in the mean. The lower section
#'    indicates which dyes are part of the dye_set. The y-axis give
#'    each dye and and the x-axis is the dye_set lining up with the
#'    dye_set on the upper section The values are represented as dots
#'    for each dye that is in the dye set and a line from the first to
#'    the last dye in the set. Columns across both sections alternate
#'    and the rows in the bottom section alternate shading.
#'
#' @param Zprime_by_plate tibble::tibble with columns
#'        dye_set  : factor this should be in the order that they will be displayed
#'                   the name of the dye set should be a string of the individual dyes
#'                   separated by the <dye_set_separator>.
#'        Zprime   : Score to be plotted should be in (-Inf, 1) where
#'                   1 is perfect, .5-1 is good, 0-0.5 is poor, and < 1 is bad
#'        plate_id : Identifiers for different plates. The Zscores should be computed
#'                   separately for each plate and then these Zprime scores will be
#'                   aggregated to make the plot
#' @param dye_set_separator separator to figure out what dyes make up the dye_set variable.
#'        Default so "_"
#' @param subtitle string subtitle for plot
#' @return ggplot2 plot
#'
#' dye_plot <- MPStats::plot_Zprime_by_dye(Zprime_by_plate)
#' ggplot2::ggsave(
#'   plot=dye_plot,
#'   filename=paste0("product/Zprime_by_dye_", MPStats::date_code(), ".pdf"),
#'   width=4, height=6,
#'   useDingbats=FALSE)
#' ggplot2::ggsave(
#'   plot=dye_plot,
#'   filename=paste0("product/Zprime_by_dye_", MPStats::date_code(), ".png"),
#'   width=4, height=6)
#'
#' @export
plot_Zprime_by_dye <- function(
  Zprime_by_plate,
  dye_set_separator="_",
  subtitle=NULL){

  Zprime_summary <- Zprime_by_plate %>%
    dplyr::group_by(dye_set) %>%
    dplyr::summarize(
      mean = mean(Zprime),
      n_plates = dplyr::n(),
      std_err = sd(Zprime) / sqrt(dplyr::n())) %>%
    dplyr::mutate(set_size = dye_set %>% stringr::str_count(dye_set_separator) + 1) %>%
    dplyr::select(-set_size) %>%
    dplyr::mutate(set_id = dplyr::row_number())

  Zprime_limits <- c(
    min(Zprime_by_plate %>% magrittr::extract2("Zprime") %>% min(), 0) - 0.05,
    1)

  sets_1hot <- Zprime_summary %>%
    dplyr::select(dye_set, set_id) %>%
    dplyr::mutate(count=1) %>%
    tidyr::separate_rows(dye_set, sep=dye_set_separator) %>%
    tidyr::pivot_wider(id_cols=set_id, names_from=dye_set, values_from=count) %>%
    tibble::column_to_rownames("set_id") %>%
    t

  # this adapts the bottom panel from UpSetR

  name_size_scale <- 2
  shade_alpha <- 0.2
  shade_color <- "gray80"
  point_size <- 5
  line_size <- 2
  dot_color <- "gray10"
  dot_alpha <- 1
  n_sets <- ncol(sets_1hot)
  n_items <- nrow(sets_1hot)

  sets_data <- expand.grid(
    y = 1:nrow(sets_1hot),
    x = 1:ncol(sets_1hot)) %>%
    dplyr::mutate(
      value = as.vector(sets_1hot),
      color = ifelse(value > 0L, dot_color, "gray83"),
      alpha = ifelse(value > 0L, 1, dot_alpha),
      intersection = ifelse(
        value > 0L,
        paste0(x, "yes"),
        paste0(row_number(), "no")))

  upper_sets_shading_data <- sets_data %>%
    dplyr::distinct(x) %>%
    dplyr::filter(x %% 2 != 0) %>%
    dplyr::mutate(
      ymin=Zprime_limits[1],
      ymax=Zprime_limits[2],
      xmin=x-0.5,
      xmax=x+0.5,
      color=shade_color)

  lower_sets_shading_data <- sets_data %>%
    dplyr::distinct(x) %>%
    dplyr::filter(x %% 2 != 0) %>%
    dplyr::mutate(
      ymin=0 + 0.5,
      ymax=n_items+0.5,
      xmin=x-0.5,
      xmax=x+0.5,
      color=shade_color)

  item_shading_data <- sets_data %>%
    dplyr::distinct(y) %>%
    dplyr::filter(y %% 2 != 0) %>%
    dplyr::mutate(
      xmin=0 + 0.5,
      xmax=max(sets_data$x) + 0.5,
      ymin=y-0.5,
      ymax=y+0.5,
      color=shade_color)

  dots_plot <- ggplot() +
    theme(
      legend.position = c(0.5, 0.9),
      legend.direction = "horizontal",
      panel.background = element_rect(fill = "white"),
      plot.margin=unit(c(-0.2,0.5,-.5,2), "lines"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(colour="grey80", size=0.5),
      panel.grid.minor.y = element_line(colour="grey90", size=0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(
        colour = "gray0", size = 7*name_size_scale, hjust = 0.4),
      axis.title.y = element_text(
        colour = "gray0", size = 8*name_size_scale,
        margin = margin(t = 0, r = -10, b = 0, l = 0))) +
    xlab(NULL) +
    scale_y_continuous(
      "Zprime Score",
      limits=Zprime_limits) +
    scale_x_continuous(
      limits=c(0, n_sets+1),
      expand = c(0, 0)) +
    scale_color_discrete("Dye Set") +
    geom_rect(
      data = upper_sets_shading_data,
      aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax),
      fill = shade_color,
      alpha = shade_alpha) +
    geom_hline(yintercept=0.5, size=0.5) +
    geom_hline(yintercept=0.0, size=0.5) +
    geom_errorbar(
      data=Zprime_summary,
      aes(x=set_id, ymin=mean-std_err, ymax=mean+std_err)) +
    geom_point(
      data=Zprime_summary,
      aes(x=set_id, y=mean),
      size=3)

  sets_plot <- ggplot() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.margin=unit(c(-0.5,0.5,0.5,0.5), "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(
        colour = "gray0", size = 7*name_size_scale, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    xlab(NULL) + ylab("   ") +
    scale_y_continuous(
      breaks = c(1:n_items),
      limits = c(0.5, n_items+0.5),
      labels = rownames(sets_1hot),
      expand = c(0,0)) +
    scale_x_continuous(
      limits = c(0, n_sets+1),
      expand = c(0, 0)) +
    geom_rect(
      data = item_shading_data,
      aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax),
      fill = shade_color,
      alpha = shade_alpha) +
    geom_rect(
      data = lower_sets_shading_data,
      aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax),
      fill = shade_color,
      alpha = shade_alpha) +
    geom_point(
      data=sets_data,
      aes(x=x, y=y),
      colour = sets_data$color,
      size= point_size,
      alpha = sets_data$alpha,
      shape=16) +
    geom_line(
      data= sets_data,
      aes(x=x, y=y,
          group = intersection,
          colour=color),
      size = line_size) +
    scale_color_identity()

    dots_plot <- dots_plot %>% ggplotGrob
    sets_plot <- sets_plot %>% ggplotGrob

    dots_plot$widths[1:5] <- sets_plot$widths[1:5]

    p <- gridExtra::grid.arrange(
      dots_plot,
      sets_plot,
      heights=c(4.5, 1.5),
      ncol=1)
}
