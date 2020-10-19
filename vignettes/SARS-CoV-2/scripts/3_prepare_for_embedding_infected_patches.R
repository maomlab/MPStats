
library(plyr)
library(tidyverse)
library(arrow)
library(caret)

data_path <- "intermediate_data/infected_patch_1999B_2020A_2021A_20201017"

viral_features <- arrow::read_parquet(
    file = paste0(data_path, "/viral_features.parquet"))
viral_feature_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_feature_columns.tsv"))
viral_metadata_columns <- readr::read_tsv(
    file = paste0(data_path, "/viral_metadata_columns.tsv"))


viral_embed_features <- viral_features %>%
    dplyr::filter(condition != "BLANK") %>%
    dplyr::transmute(
        plate_id,
        dplyr::across(
            .cols = viral_feature_columns$feature,
            .fns = function(feature) {
                transform <- viral_feature_columns %>%
                    dplyr::filter(feature == dplyr::cur_column()) %>%
                    dplyr::pull(transform)
                if (transform == "log") {
                    feature <- log(feature)
                } else if (transform == "log1p") {
                    feature <- log(feature + 1)
                } else if (transform == "identity") {
                    feature <- feature
                } else {
                    stop(paste0("Unrecongized transform '", transform, "'", sep = ""))
                }
            })) %>%
    dplyr::group_by(plate_id) %>%
    dplyr::mutate_at(viral_feature_columns$feature, ~ scale(.)[, 1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-plate_id)

viral_embed_features %>%
    arrow::write_parquet(
        sink = paste0(
            data_path,
            "/viral_plate_scaled_MasterDataTable.parquet"))



system(paste0("
        cd ~/opt/MPLearn/vignettes/SARS-CoV-2/S25 &&
        /home/ubuntu/anaconda3/envs/sextonlab/bin/python \\
            ~/anaconda3/envs/sextonlab/bin/embed_umap \\
            --dataset ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_plate_scaled_MasterDataTable.parquet \\
            --tag UMAP_viral_plate_scaled_neighbors=15_neg_sampling_rate=20_epochs=2000_20201019 \\
            --feature_columns ~/opt/MPStats/vignettes/SARS-CoV-2/", data_path, "/viral_feature_no_transform_columns.tsv \\
            --no_standardize_features \\
	    --umap_low_memory \\
            --umap_negative_sample_rate 20 \\
            --umap_n_epochs 2000 \\
	    --verbose
"))


source("scripts/monocle3_support.R")
infected_cds <- populate_cds(
    cell_features = viral_embed_features,
    cell_feature_columns = viral_feature_columns,
    cell_metadata_columns = viral_metadata_columns,
    embedding_type = c("UMAP"),
    embedding = infected_cell_features %>% dplyr::select(UMAP_1, UMAP_2),
    verbose = TRUE)

infected_cds %>% serialize_clusters(
    output_fname = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_top_hits_infected_plate_scaled_epochs=2000_re_embed_epochs=2000_200730/clusters_leiden_res=5e-7.parquet")
# as resolution gets bigger --> more clusters
infected_cds <- infected_cds %>%
    monocle3::cluster_cells(
        k = 200,
        resolution = .00001,
        num_iter = 10,
        verbose = TRUE)
infected_cds %>% serialize_clusters(
    output_fname = "~/opt/MPLearn/vignettes/SARS-CoV-2/S25/intermediate_data/UMAP_embedding_top_hits_infected_plate_scaled_epochs=2000_re_embed_epochs=2000_200730/clusters_leiden_k=200_res=1e-5.parquet")
