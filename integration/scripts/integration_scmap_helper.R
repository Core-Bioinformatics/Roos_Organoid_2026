library(scmap)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(ClustAssess)

ranking_functions <- list(
    "median" = median,
    "max" = max,
    "top_qt" = function(arr) {
        fivenum(arr)[4]
    },
    "top_qt_max" = function(arr) {
        arr_summary <- fivenum(arr)
        arr_summary[4] + arr_summary[5]
    },
    "iqr" = function(arr) {
        arr_summary <- fivenum(arr)
        arr_summary[2] - arr_summary[4]
    },
    "iqr_median" = function(arr) {
        arr_summary <- fivenum(arr)
        arr_summary[3] - (arr_summary[4] - arr_summary[2])
    },
    "iqr_median_coeff" = function(arr) {
        arr_summary <- fivenum(arr)
        arr_summary[3] / (arr_summary[4] - arr_summary[2])
    },
    "mean" = mean
)

rank_configs <- function(ecc_list, rank_by = "top_qt_max", return_type = "order") {
    if (!(rank_by %in% names(ranking_functions))) {
        rank_by <- "top_qt_max"
    }

    ranking_fun <- ranking_functions[[rank_by]]
    scores <- sapply(ecc_list, ranking_fun)
    # print(scores)

    if (return_type == "order") {
        return(order(scores, decreasing = TRUE))
    }

    return(length(ecc_list) + 1 - rank(scores, ties.method = "first"))
}

nn_rank_configurations <- function(nn_stability_result, ranking_criterion = "iqr") {
    config_names <- names(nn_stability_result$n_neigh_ec_consistency)

    nn_ecc_list <- list()
    nn_config_list <- c()
    nn_list <- c()
    index <- 1

    for (config_name in config_names) {
        n_neigh_values <- names(nn_stability_result$n_neigh_ec_consistency[[config_name]])
        for (n_neigh in n_neigh_values) {
            nn_ecc_list[[index]] <- nn_stability_result$n_neigh_ec_consistency[[config_name]][[n_neigh]]
            index <- index + 1
            nn_config_list <- c(nn_config_list, config_name)
            nn_list <- c(nn_list, n_neigh)
        }
    }

    ecc_medians <- sapply(nn_ecc_list, median)
    eligible_configs <- which(ecc_medians >= fivenum(ecc_medians)[4])
    best_config_index <- rank_configs(nn_ecc_list[eligible_configs], rank_by = ranking_criterion, return_type = "order")[1]
    best_config <- nn_config_list[eligible_configs][best_config_index]
    best_nn <- as.integer(nn_list[eligible_configs][best_config_index])

    return(list(best_config, best_nn))
}

build_adj_matrix <- function(embedding, n_neigh, graph_type = c("nn", "snn")) {
    if (length(graph_type) > 1) {
        graph_type <- graph_type[1]
    }

    n_neigh <- as.integer(n_neigh)

    adj_matrix <- Seurat::FindNeighbors(
        object = embedding,
        k.param = n_neigh,
        compute.SNN = FALSE,
        nn.method = "rann",
        verbose = FALSE
    )$nn

    if (graph_type == "snn") {
        adj_matrix_snn <- ClustAssess::get_highest_prune_param(adj_matrix, n_neigh)
        rownames(adj_matrix_snn$adj_matrix) <- rownames(embedding)
        colnames(adj_matrix_snn$adj_matrix) <- rownames(embedding)

        adj_matrix_snn$adj_matrix <- Seurat::as.Graph(adj_matrix_snn$adj_matrix)

        return(adj_matrix_snn)
    }

    rownames(adj_matrix) <- rownames(embedding)
    colnames(adj_matrix) <- rownames(embedding)

    return(list(adj_matrix = adj_matrix, prune_value = 0))
}

count_dataset_specific_cells <- function(contingency_table, threshold = 5) {
    nclusters <- ncol(contingency_table)
    nclasses <- nrow(contingency_table)
    ncells <- rep(0, nclasses)
    names(ncells) <- rownames(contingency_table)

    for (i in seq_len(nclusters)) {
        dominant_class <- which.max(contingency_table[, i])
        all_small_influence <- TRUE
        for (j in seq_len(nclasses)) {
            if (j == dominant_class) {
                next
            }

            if (contingency_table[j, i] >= threshold) {
                all_small_influence <- FALSE
                break
            }
        }

        if (all_small_influence) {
            ncells[dominant_class] <- ncells[dominant_class] + contingency_table[dominant_class, i]
        }
    }

    return(ncells)
}

get_nn_from_umap <- function(so, umap_arguments = list()) {
    pca_emb <- so@reductions$pca@cell.embeddings

    do.call(
        uwot::umap,
        c(
            list(X = pca_emb, ret_model = TRUE, ret_nn = TRUE),
            umap_arguments
        )
    )
}

get_scmap_similarities <- function(so1,
                                   so2,
                                   k = 30,
                                   nfeatures = 2000,
                                   features1 = NULL,
                                   features2 = NULL,
                                   seed = 42) {
    sce1 <- Seurat::as.SingleCellExperiment(so1)
    SummarizedExperiment::rowData(sce1)$feature_symbol <- rownames(sce1)
    sce2 <- Seurat::as.SingleCellExperiment(so2)
    SummarizedExperiment::rowData(sce2)$feature_symbol <- rownames(sce2)

    sce1 <- scmap::selectFeatures(sce1, n_features = nfeatures, suppress_plot = TRUE)
    sce2 <- scmap::selectFeatures(sce2, n_features = nfeatures, suppress_plot = TRUE)

    if (!is.null(features1)) {
        SummarizedExperiment::rowData(sce1)$scmap_features <- rownames(sce1) %in% features1
    }

    if (!is.null(features2)) {
        SummarizedExperiment::rowData(sce2)$scmap_features <- rownames(sce2) %in% features2
    }

    set.seed(seed)
    sce1 <- scmap::indexCell(sce1, M = NULL, k = NULL)
    set.seed(seed)
    sce2 <- scmap::indexCell(sce2, M = NULL, k = NULL)

    set.seed(seed)
    proj12 <- scmap::scmapCell(
        projection = sce1,
        index_list = list(
            sce2 = S4Vectors::metadata(sce2)$scmap_cell_index
        ),
        w = k
    )$sce2

    # convert to distance
    proj12$similarities <- 1 - proj12$similarities
    proj12$cells <- t(proj12$cells)
    proj12$similarities <- t(proj12$similarities)

    set.seed(seed)
    proj21 <- scmap::scmapCell(
        projection = sce2,
        index_list = list(
            sce1 = S4Vectors::metadata(sce1)$scmap_cell_index
        ),
        w = k
    )$sce1

    proj21$similarities <- 1 - proj21$similarities
    proj21$cells <- t(proj21$cells)
    proj21$similarities <- t(proj21$similarities)

    return(list(proj12, proj21))
}

get_integrated_nn <- function(so1,
                              so2,
                              k = 30,
                              nfeatures = 2000,
                              features1 = NULL,
                              features2 = NULL,
                              umap_arguments = list()) {
    n1 <- ncol(so1)

    nn1 <- get_nn_from_umap(so1, umap_arguments = umap_arguments)$nn[[1]]
    nn2 <- get_nn_from_umap(so2, umap_arguments = umap_arguments)$nn[[1]]

    # add the number of rows of the first dataset to the index of the second dataset
    nn2$idx <- nn2$idx + n1

    scmap_results <- get_scmap_similarities(
        so1,
        so2,
        k = k,
        nfeatures = nfeatures,
        features1 = features1,
        features2 = features2
    )

    # add the number of rows of the first dataset to the index of the second dataset
    scmap_results[[1]]$cells <- scmap_results[[1]]$cells + n1

    # merge the nn matrices inside and between the datasets
    final_nn <- list(
        idx = cbind(
            rbind(nn1$idx, nn2$idx),
            rbind(scmap_results[[1]]$cells, scmap_results[[2]]$cells)
        ),
        dist = cbind(
            rbind(nn1$dist, nn2$dist),
            rbind(scmap_results[[1]]$similarities, scmap_results[[2]]$similarities)
        )
    )

    # order the indices by the distance
    for (i in seq_len(nrow(final_nn$idx))) {
        actual_order <- order(final_nn$dist[i, ])
        final_nn$idx[i, ] <- final_nn$idx[i, actual_order]
        final_nn$dist[i, ] <- final_nn$dist[i, actual_order]
    }

    return(final_nn)
}

trim_nn <- function(nn, k = 30) {
    nn$idx <- nn$idx[, seq_len(k)]
    nn$dist <- nn$dist[, seq_len(k)]
    return(nn)
}

get_umap_df <- function(umap,
                        name_first,
                        name_second,
                        clusters_first,
                        clusters_second) {
    df <- data.frame(umap)
    df$dataset <- name_second
    df$dataset[seq_along(clusters_first)] <- name_first

    df$clusters <- c(
        paste0(name_first, "_", clusters_first),
        paste0(name_second, "_", clusters_second)
    )
    df$dataset <- factor(df$dataset)
    df$clusters <- factor(df$clusters)

    return(df)
}


assess_integration_quality <- function(dataset_list,
                                       k,
                                       target_clusters,
                                       nfeatures = 2000,
                                       id = "",
                                       features_list = NULL,
                                       umap_arguments = list(seed = 42, min_dist = 0.3, verbose = FALSE),
                                       output_path = ".",
                                       n_reps = 50,
                                       n_neigh_seq = seq(from = 5, to = 50, by = 5),
                                       resolution_seq = seq(from = 0.1, to = 2, by = 0.1),
                                       threshold = 5) {
    dataset_names <- names(dataset_list)
    if (id == "") {
        id <- paste(dataset_names, collapse = "_")
    }

    if (length(k) > 1) {
        results_list <- lapply(
            k,
            function(k_val) {
                assess_integration_quality(
                    dataset_list = dataset_list,
                    k = k_val,
                    target_clusters = target_clusters,
                    nfeatures = nfeatures,
                    id = id,
                    features_list = features_list,
                    umap_arguments = umap_arguments,
                    output_path = output_path,
                    n_reps = n_reps,
                    n_neigh_seq = n_neigh_seq,
                    resolution_seq = resolution_seq,
                    threshold = threshold
                )
            }
        )
        # convert results_list to data.frame
        names(results_list) <- k
        results_list <- as.data.frame(do.call(rbind, results_list))

        for (dts in dataset_names) {
            results_list[[glue::glue("{dts}_percentage")]] <- results_list[[dts]] / ncol(dataset_list[[dts]])
        }

        return(results_list)
    }

    folder_path <- file.path(output_path, glue::glue("integration_quality/{id}"))
    if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
    }

    dir.create(file.path(folder_path, "contingency_heatmaps"), recursive = TRUE)
    dir.create(file.path(folder_path, "umaps"), recursive = TRUE)
    dir.create(file.path(folder_path, "clustassess"), recursive = TRUE)

    if (is.null(features_list)) {
        features_list <- list()
        for (dataset_name in dataset_names) {
            features_list[[dataset_name]] <- NULL
        }
    }

    umap_arguments["n_neighbors"] <- k
    umap_arguments["metric"] <- "cosine"

    merged_nn <- get_integrated_nn(
        dataset_list[[1]],
        dataset_list[[2]],
        k = k,
        nfeatures = nfeatures,
        features1 = features_list[[1]],
        features2 = features_list[[2]],
        umap_arguments = umap_arguments
    )

    umap_arguments["n_neighbors"] <- 2 * k
    umap_df <- get_umap_df(
        do.call(
            uwot::umap,
            c(
                list(X = NULL, nn_method = merged_nn),
                umap_arguments
            )
        ),
        dataset_names[[1]],
        dataset_names[[2]],
        dataset_list[[1]]$target_clusters,
        dataset_list[[2]]$target_clusters
    )

    saveRDS(umap_df, file.path(folder_path, glue::glue("umaps/k_{k}_umap_df.rds")))

    graph_constr_assessment <- assess_nn_stability(
        embedding = as.matrix(umap_df[, c("X1", "X2")]),
        n_neigh_sequence = n_neigh_seq,
        n_repetitions = n_reps,
        graph_reduction_type = "PCA"
    )

    saveRDS(graph_constr_assessment, file.path(folder_path, glue::glue("clustassess/k_{k}_graph_constr_assessment.rds")))

    nn_ranking <- nn_rank_configurations(graph_constr_assessment)
    split_configs <- strsplit(nn_ranking[[1]], "_")[[1]]
    adj_matrix_list <- build_adj_matrix(as.matrix(umap_df[, c("X1", "X2")]), nn_ranking[[2]], graph_type = split_configs[2])

    pdf(file.path(folder_path, glue::glue("clustassess/k_{k}_graph_constr_assessment.pdf")), width = 12, height = 6)
    print(plot_n_neigh_ecs(graph_constr_assessment))
    dev.off()

    graph_clustering_assessment <- assess_clustering_stability(
        graph_adjacency_matrix = adj_matrix_list$adj_matrix,
        resolution = resolution_seq,
        n_repetitions = n_reps
    )

    saveRDS(graph_clustering_assessment, file.path(folder_path, glue::glue("clustassess/k_{k}_graph_clustering_assessment.rds")))

    pdf(file.path(folder_path, glue::glue("clustassess/k_{k}_graph_clustering_assessment.pdf")), width = 12, height = 6)
    print(plot_clustering_overall_stability(graph_clustering_assessment))
    dev.off()

    graph_clustering_assessment$split_by_k$Louvain <- NULL
    graph_clustering_assessment$split_by_k$Louvain.refined <- NULL
    graph_clustering_assessment$split_by_resolution$Louvain <- NULL
    graph_clustering_assessment$split_by_resolution$Louvain.refined <- NULL

    pdf(file.path(folder_path, glue::glue("clustassess/k_{k}_nclusters_assessment_k.pdf")), width = 25, height = 10)
    print(plot_k_n_partitions(graph_clustering_assessment))
    dev.off()

    while (TRUE) {
        if (as.character(target_clusters) %in% names(graph_clustering_assessment$split_by_k$SLM)) {
            break
        }

        target_clusters <- target_clusters + 1
    }

    mb <- graph_clustering_assessment$split_by_k$SLM[[as.character(target_clusters)]]$partitions[[1]]$mb
    umap_df$mb <- factor(mb)

    pdf(file.path(folder_path, glue::glue("umaps/k_{k}_clusters.pdf")), width = 10, height = 10)
    print(ggplot(umap_df, aes(x = .data$X1, y = .data$X2, color = .data$mb)) +
        geom_point(size = 0.5) +
        ggtitle(glue::glue("k = {k} all neighbours\nSLM {target_clusters} clusters")))
    print(ggplot(umap_df, aes(x = .data$X1, y = .data$X2, color = .data$dataset)) +
        geom_point(size = 0.5) +
        ggtitle(glue::glue("k = {k} all neighbours\nDatasets")))
    print(ggplot(umap_df, aes(x = .data$X1, y = .data$X2, color = .data$clusters)) +
        geom_point(size = 0.5) +
        ggtitle(glue::glue("k = {k} all neighbours\nOriginal lusters")))
    dev.off()

    contingency_table <- table(umap_df$clusters, umap_df$mb)
    pdf(file.path(folder_path, glue::glue("contingency_heatmaps/k_{k}_cont_table_clusters.pdf")), width = 12, height = 6)
    pheatmap::pheatmap(
        contingency_table,
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
        cellwidth = 30,
        cellheight = 20,
        display_numbers = TRUE,
        number_format = "%.0f",
        main = glue::glue("Contingency table - k = {k}")
    )
    dev.off()

    contingency_table <- table(umap_df$dataset, umap_df$mb)
    pdf(file.path(folder_path, glue::glue("contingency_heatmaps/k_{k}_cont_table_datasets.pdf")), width = 12, height = 6)
    pheatmap::pheatmap(
        contingency_table,
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
        cellwidth = 30,
        cellheight = 20,
        display_numbers = TRUE,
        number_format = "%.0f",
        main = glue::glue("Contingency table - k = {k}")
    )
    dev.off()

    return(count_dataset_specific_cells(contingency_table, threshold))
}
