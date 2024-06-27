source("integration_scmap_helper.R")
# ============== CASE 1 ==============
# ORGANOID 8w vs RAWLINGS 6w

organoid <- readRDS(file.path(organoid_object_folder, "R_objects/aggregated/HV/by_organ_age_condition/aggr_filtered_coding_highrp-organ-Lung-age-08w-condition--Wnt_without_oesophagus_mesenchymal.rds"))
organoid_ca <- readRDS(file.path(organoid_object_folder, "R_objects/clustassess/by_organ_age_condition/organ-Lung-age-08w-condition--Wnt_without_oesophagus_mesenchymal.rds"))
organoid_hv_genes <- organoid_ca[["Highly_Variable"]]$feature_list
organoid_ma_genes <- organoid_ca[["Most_Abundant"]]$feature_list
organoid <- RunPCA(organoid)
organoid@reductions$pca@cell.embeddings <- organoid_ca$Most_Abundant$"4000"$pca
organoid$target_clusters <- organoid_ca$Most_Abundant$"4000"$clustering_stability$"split_by_k"$SLM$"4"$partitions[[1]]$mb

rawlings <- readRDS(file.path(rawlings_object_folder, "data/R_objects/6w_epithelial.rds"))
rawlings_ca <- readRDS(file.path(rawlings_object_folder, "data/R_objects/clustassess/6w_epithelial.rds"))
rawlings_hv_genes <- rawlings_ca[["Highly_Variable"]]$feature_list
rawlings_ma_genes <- rawlings_ca[["Most_Abundant"]]$feature_list
rawlings <- RunPCA(rawlings)
rawlings@reductions$pca@cell.embeddings <- rawlings_ca$Most_Abundant$"1500"$pca
rawlings$target_clusters <- rawlings_ca$Most_Abundant$"1500"$clustering_stability$"split_by_k"$SLM$"6"$partitions[[1]]$mb

psocket <- parallel::makePSOCKcluster(50)
doParallel::registerDoParallel(psocket)

k_values <- c(seq(from = 5, to = 50, by = 5), seq(from = 60, to = 100, by = 10))
integration_quality_assessment <- assess_integration_quality(
    dataset_list = list(r6w = rawlings, o8w = organoid),
    k = k_values,
    target_clusters = 20,
    nfeatures = 2000,
    id = "r6w_vs_o8w",
    features_list = list(r6w = rawlings_ma_genes[seq_len(1500)], organoid_ma_genes[seq_len(4000)]),
    umap_arguments = list(seed = 42, min_dist = 0.3, verbose = FALSE),
    output_path = output_path,
    n_reps = 50,
    n_neigh_seq = seq(from = 5, to = 50, by = 5),
    resolution_seq = seq(from = 0.01, to = 0.5, by = 0.01),
    threshold = 5
)

# stop clusters
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()

saveRDS(integration_quality_assessment, file = file.path(output_path, "integration_quality", "r6w_vs_o8w", "dataset_specific_cells.rds"))
# integration_quality_assessment <- readRDS(file = file.path(output_path, "integration_quality", "r6w_vs_o8w","dataset_specific_cells.rds"))

display_matrix <- as.matrix(integration_quality_assessment[, 3:4])
colnames(display_matrix) <- c("r6w", "o8w")

integration_quality_assessment$display_r <- glue::glue("{sprintf('%.04f', integration_quality_assessment$r6w_percentage)}\n({integration_quality_assessment$r6w})")
integration_quality_assessment$display_o <- glue::glue("{sprintf('%.04f', integration_quality_assessment$o8w_percentage)}\n({integration_quality_assessment$o8w})")
# create heatmap of the integration qualtiy and insert the values in the cells
pheatmap::pheatmap(
    display_matrix,
    color = colorRampPalette(c("white", "orange"))(n = 299),
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 10,
    fontsize_number = 10,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = "black",
    cellwidth = 50,
    cellheight = 40,
    display_numbers = as.matrix(integration_quality_assessment)[, 5:6],
    number_format = "%.04f",
    main = glue::glue("Integration quality\nassessment - r6w vs o8w\n"),
    angle_col = 0,
    file = file.path(output_path, "integration_quality", "r6w_vs_o8w", "integration_quality_assessment.pdf"),
    width = 10,
    height = 10
)

# ============== CASE 2 ==============
# ORGANOID 12w vs RAWLINGS 8w

organoid <- readRDS(file.path(organoid_object_folder, "R_objects/aggregated/HV/by_organ_age_condition/aggr_filtered_coding_highrp-organ-Lung-age-12w-condition--Wnt_without_oesophagus_mesenchymal_low_ncount_cluster.rds"))
organoid_ca <- readRDS(file.path(organoid_object_folder, "R_objects/clustassess/by_organ_age_condition/organ-Lung-age-12w-condition--Wnt_without_oesophagus_mesenchymal_low_ncount_cluster.rds"))
organoid_hv_genes <- organoid_ca[["Highly_Variable"]]$feature_list
organoid_ma_genes <- organoid_ca[["Most_Abundant"]]$feature_list
organoid <- RunPCA(organoid)
organoid@reductions$pca@cell.embeddings <- organoid_ca$Highly_Variable$"4500"$pca
organoid$target_clusters <- organoid_ca$Highly_Variable$"4500"$clustering_stability$"split_by_k"$SLM$"6"$partitions[[1]]$mb

rawlings <- readRDS(file.path(rawlings_object_folder, "data/R_objects/8w_epithelial.rds"))
rawlings_ca <- readRDS(file.path(rawlings_object_folder, "data/R_objects/clustassess/8w_epithelial.rds"))
rawlings_hv_genes <- rawlings_ca[["Highly_Variable"]]$feature_list
rawlings_ma_genes <- rawlings_ca[["Most_Abundant"]]$feature_list
rawlings <- RunPCA(rawlings)
rawlings@reductions$pca@cell.embeddings <- rawlings_ca$Highly_Variable$"3000"$pca
rawlings$target_clusters <- rawlings_ca$Highly_Variable$"3000"$clustering_stability$"split_by_k"$SLM$"6"$partitions[[1]]$mb
rawlings <- RunUMAP(rawlings, reduction = "pca", dims = 1:30)
rawlings@reductions$umap@cell.embeddings <- rawlings_ca$Highly_Variable$"3000"$umap

psocket <- parallel::makePSOCKcluster(50)
doParallel::registerDoParallel(psocket)

k_values <- c(seq(from = 5, to = 50, by = 5), seq(from = 60, to = 100, by = 10))
integration_quality_assessment <- assess_integration_quality(
    dataset_list = list(r8w = rawlings, o12w = organoid),
    k = k_values,
    target_clusters = 20,
    nfeatures = 2000,
    id = "r8w_vs_o12w",
    features_list = list(r8w = rawlings_hv_genes[seq_len(3000)], o12w = organoid_hv_genes[seq_len(4500)]),
    umap_arguments = list(seed = 42, min_dist = 0.3, verbose = FALSE),
    output_path = output_path,
    n_reps = 50,
    n_neigh_seq = seq(from = 5, to = 50, by = 5),
    resolution_seq = seq(from = 0.01, to = 1, by = 0.01),
    threshold = 5
)

# stop clusters
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()

saveRDS(integration_quality_assessment, file = file.path(output_path, "integration_quality", "r8w_vs_o12w", "dataset_specific_cells.rds"))
integration_quality_assessment <- readRDS(file = file.path(output_path, "integration_quality", "r8w_vs_o12w", "dataset_specific_cells.rds"))

display_matrix <- as.matrix(integration_quality_assessment[, 3:4])
colnames(display_matrix) <- c("r8w", "o12w")

integration_quality_assessment$display_r <- glue::glue("{sprintf('%.04f', integration_quality_assessment$r8w_percentage)}\n({integration_quality_assessment$r8w})")
integration_quality_assessment$display_o <- glue::glue("{sprintf('%.04f', integration_quality_assessment$o12w_percentage)}\n({integration_quality_assessment$o12w})")
# create heatmap of the integration qualtiy and insert the values in the cells
pheatmap::pheatmap(
    display_matrix,
    color = colorRampPalette(c("white", "orange"))(n = 299),
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 10,
    fontsize_number = 10,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    border_color = "black",
    cellwidth = 50,
    cellheight = 40,
    display_numbers = as.matrix(integration_quality_assessment)[, 5:6],
    number_format = "%.04f",
    main = glue::glue("Integration quality\nassessment - r8w vs o12w\n"),
    angle_col = 0,
    file = file.path(output_path, "integration_quality", "r8w_vs_o12w", "integration_quality_assessment.pdf"),
    width = 10,
    height = 10
)
