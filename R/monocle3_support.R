
#' Popluate a Monocle3 cell data set from cell features
#'
#' The Monocle3 cell data set is a container
#'  
#' @param cell_features data.frame:
#'        rows: cells, columns: features
#' @param cell_metadata_columns data.frame:
#'        rows: features, columns feature metadata
#'        Note that there must be a column named `feature`
#' @param cell_metadata
#'        rows: cells, columns cell metadata
#'
#' @export
populate_cds <- function(
    cell_features,
    cell_feature_columns,
    cell_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = NULL,
    verbose = FALSE) {

    assertthat::assert_that("feature" %in% names(cell_feature_columns))

    n_cells <- nrow(cell_features)
    n_features <- nrow(cell_feature_columns)
    cat("Loading cell dataset with dimensions [<feature>, <cell>] = [", n_features, ", ",  n_cells, "]\n", sep = "")
    
    expression_data <- cell_features %>%
        dplyr::select(
            tidyselect::one_of(cell_feature_columns$feature)) %>%
        as.matrix() %>%
        t()
    gene_metadata <- cell_feature_columns %>%
        dplyr::mutate(
            gene_short_name = feature)
    cell_metadata <- cell_features %>%
        dplyr::select(
            tidyselect::one_of(cell_metadata_columns$feature))

    row.names(expression_data) <- cell_feature_columns$feature
    row.names(gene_metadata) <- cell_feature_columns$feature
    names(expression_data) <- expression_data %>% ncol %>% seq_len
    row.names(cell_metadata) <- expression_data %>% ncol %>% seq_len

    if(verbose){
        cat("Creating a SingleCellExperiment object ...\n")
    }
    # unpack monocle3::new_cell_data_set(...)
    # to not use dgCMatrix for the expression matrix they are dense feature matrices
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = expression_data),
        rowData = gene_metadata,
        colData = cell_metadata)

    if(verbose){
        cat("Creating a Cell Data Set object ...\n")
    }
    cds <- methods::new(
        Class = methods::getClass("cell_data_set", where = "monocle3"),
        assays = SummarizedExperiment::Assays(list(counts = expression_data)),
        colData = colData(sce),
        int_elementMetadata = int_elementMetadata(sce),
        int_colData = int_colData(sce),
        int_metadata = int_metadata(sce),
        metadata = S4Vectors::metadata(sce),
        NAMES = NULL,
        elementMetadata = elementMetadata(sce)[,0],
        rowRanges = rowRanges(sce))

    if(verbose){
        cat("Configuring the cell data set ...\n")
    }

    S4Vectors::metadata(cds)$cds_version <- Biobase::package.version("monocle3")
    clusters <- stats::setNames(S4Vectors::SimpleList(), character(0))
    
    # we want to do this but there is a bug where it fails to detect that it
    # SingleCellExperiment::counts(cds) is not sparse
    # cds <- monocle3::estimate_size_factors(cds)
    if (any(Matrix::colSums(SingleCellExperiment::counts(cds)) == 0)) {
      warning(
        "Your CDS object contains cells with zero reads. ", 
        "This causes size factor calculation to fail. Please remove ", 
        "the zero read cells using ",
        "cds <- cds[,Matrix::colSums(exprs(cds)) != 0] and then ", 
        "run cds <- estimate_size_factors(cds)")
    }
    cell_total <- cds |>
      SingleCellExperiment::counts() |>
      round() |>
      apply(2, sum)
    # "mean-geometric-mean-total"
    sf <- cell_total/exp(mean(log(cell_total)))
    sf[is.na(sf)] <- 1
    SummarizedExperiment::colData(cds)$Size_Factor <- sf

    row.names(SummarizedExperiment::colData(cds)) <- expression_data %>% ncol %>% seq_len
    if (!is.null(embedding)) {
        SingleCellExperiment::reducedDims(cds)[[embedding_type]] <- embedding
    }
    cds
}


#' Write out clusters from monocle3 cell dataset
#'
#' @param cds monocle3 cell dataset
#' @param output_fname path to output parquet file where the clusters should be written
#'                     the data has a single column [cluster_label] with values for each object
#' @param reduction_method the reduction method for which the clusters were computed [default: UMAP]
#' @param verbose verbose output [default: FALSE]
#'
#' @export
serialize_clusters <- function(
    cds,
    output_fname,
    reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
    verbose = FALSE) {

    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(is.logical(verbose))

    assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
        msg = paste("No cell clusters for", reduction_method,
            "calculated.", "Please run cluster_cells with", "reduction_method =",
            reduction_method, "before trying to serialize clusters."))

    if(verbose) {
        message(
            "Writing clusters for cell data set reduced by ",
            "'", reduction_method, "' reduction method to ",
            "'", output_fname, "'\n", sep = "")
    }
    
    cds@clusters[[reduction_method]]$clusters %>%
        data.frame(cluster_label = .) %>%
        arrow::write_parquet(output_fname)
}
    

#' Compute important quality control covariates suggested by (Lueken and Theis 2019)
#'
#' Add count_depth, n_genes, mt_fraction values to
#' the column data for the input cell_data_set
#'
#' @param cds cell_data_set
#'@export
compute_qc_covariates <- function(cds) {
    count_depth <- cds %>%
        SingleCellExperiment::counts() %>%
        Matrix::colSums()
    SummarizedExperiment::colData(cds)[["count_depth"]] <- count_depth
    
    # compute number of non-zero entries per-column for a dgCMatrix
    # https://stackoverflow.com/a/51560622/198401
    SummarizedExperiment::colData(cds)[["n_genes"]] <- SingleCellExperiment::counts(cds)@p %>% diff()
    mt_genes <- SummarizedExperiment::rowData(cds) %>%
        data.frame %>%
        dplyr::mutate(index = dplyr::row_number()) %>%
        dplyr::filter(gene_short_name %>% stringr::str_starts("MT-"))
         
    mt_count_depth <- cds[mt_genes$index, ] %>%
        SingleCellExperiment::counts() %>%
        Matrix::colSums()

    SummarizedExperiment::colData(cds)[["mt_fraction"]] <- mt_count_depth / count_depth
    cds
}
