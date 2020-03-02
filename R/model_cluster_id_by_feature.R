

model_cluster_label_by_feature_lm <- function(cell_features, cluster_labels){
    data <- cell_features %>%
        dplyr::mutate(cluster_label = as.factor(cluster_labels$cluster_label))
    model <- lm(cluster_label ~ ., data=data)
}


bstDense <- xgboost(
    data = as.matrix(cell_features),
    label=cluster_labels$cluster_label == 39,
    max.depth = 2,
    eta = 1,
    nthread = 40,
    nrounds = 2,
    objective = "binary:logistic")


bstDense <- xgboost(
    data = as.matrix(cell_features),
    label=cluster_labels$cluster_label == 2,
    max.depth = 2,
    eta = 1,
    nthread = 40,
    nrounds = 2,
    objective = "binary:logistic")
                                    Feature        Gain      Cover Frequency
1:     Intensity_IntegratedIntensity_Hoe_ER 0.837153462 0.34957503 0.1666667
2:             Texture_SumAverage_CMO_20_00 0.119862828 0.15042497 0.1666667
3:          Texture_Correlation_Lipids_8_00 0.022304826 0.04777903 0.1666667
4:          Texture_SumAverage_Hoe_ER_20_00 0.011423058 0.30179601 0.1666667
5:          Texture_Correlation_Hoe_ER_5_00 0.006264655 0.02315743 0.1666667
6: Intensity_IntegratedIntensityEdge_Hoe_ER 0.002991171 0.12726753 0.1666667



importance_matrix <- xgb.importance(model = bstDense)
print(importance_matrix)
