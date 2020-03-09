#' fit principle curves
#'
#' @param cell_embedding tibble::tibble with columns
#'    UMAP_1
#'    UMAP_2
#'    cluster_label
#' @param thin number of samples for each cluster_label to fit curve to
#'        if false, use the full dataset
#' @param princurve_smoother default: lowess (see princurve::principal_curve)
#' @param princurve_maxit default: 300       (see princurve::principal_curve)
#' @param princurve_stretch default: 0       (see princurve::principal_curve)
#' 
#' @return tibble::tibble representing the principal curve for each cluster
#'         for each cluster
#'    UMAP_1                       input embedding coordinate
#'    UMAP_2                       input embedding coordiante
#'    cluster_label                input cluster identifier
#'    princurve_1                  UMAP_1 coordiante for the projection on to the principal curve
#'    princurve_2                  UMAP_2 coordinate for the projection on to the principal curve
#'    princurve_labmda             distance along the principal curve of the projection
#'    princurve_order              order of projection along the principal curve
#'
#' @usage To use this for example,
#' 
#'  cluster_principal_curves <- dplyr::bind_cols(
#'    arrow::read_parquet("input/umap_embedding.parquet"),
#'    arrow::read_parquet("input/hbscan_clustering.parquet")) %>%
#'    MPStats::fit_principal_curves() %>%
#'    dplyr::arrange(princurve_order)
#' 
#' @export
fit_principal_curves <- function(
   cell_embedding,
   thin=2000,                              
   princurve_smoother="lowess",
   princurve_maxit=300,
   princurve_stretch=0,
   verbose=TRUE,
   ...){          

  cell_embedding %>%
    plyr::ddply("cluster_label", function(cluster_embedding){
          
    cluster_label <- cluster_embedding$cluster_label[1] %>% as.numeric()

    if(verbose){
      cat(
        "Fitting curve for cluster: ", cluster_label,
        " with ", nrow(cluster_embedding), " points\n", sep="")
    }  

    if(thin && (thin < nrow(cluster_embedding))){
      fit_data <- cluster_embedding %>% dplyr::sample_n(thin)
    } else{
      fit_data <- cluster_embedding
    }

    fit <-  fit_data %>%
      dplyr::select(UMAP_1, UMAP_2) %>%
      as.matrix() %>%
      princurve::principal_curve(
        smoother=princurve_smoother,
        maxit=princurve_maxit,
        stretch=princurve_stretch,
        trace=verbose,
        ...)

    if(verbose){
      if(!fit$converged) {
        cat("Principal curve did not converge after ", fit$num_iterations , " iterations.\n", sep="")
        
        if(fit$num_iteractions < princurve_maxit){
          cat("Consider increasing the number of maximum number of itertions by setting the princurve_maxit parameter.\n")
        }
      }
    }

    if(thin && (thin < nrow(cluster_embedding))){
      fit <- cluster_embedding %>%
        dplyr::select(UMAP_1, UMAP_2) %>%
        as.matrix() %>%
        princurve::project_to_curve(
          s=fit$s,
          stretch=princurve_stretch)             
    }

        
    cluster_embedding %>%
      dplyr::mutate(
        princurve_1 = fit$s[,1],
        princurve_2 = fit$s[,2],
        princurve_lambda = fit$lambda,
        princurve_order = fit$ord)
  })
}
