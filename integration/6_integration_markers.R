library(Seurat)
library(dplyr)


find_markers <- function(so,
                         umap_df,
                         n_chosen_clusters,
                         id,
                         marker_folder,
                         threshold_cells = 3) {
    nclusters_dataset <- length(unique(umap_df[umap_df$dataset == id, "clusters"]))
    marker_folder <- file.path(marker_folder, glue::glue("k_{chosen_k}"), id)
    dir.create(marker_folder, recursive = TRUE, showWarnings = FALSE)

    for (j in seq_len(n_chosen_clusters)) {
        for (i in seq_len(nclusters_dataset)) {
            group_1 <- rownames(umap_df)[which(umap_df$clusters == glue::glue("{id}_{i}") & umap_df$mb == glue::glue("{j}"))]

            if (length(group_1) < threshold_cells) {
                next
            }

            Idents(so) <- factor(ifelse(
                colnames(so) %in% group_1,
                glue::glue("{id}_{i}_integrated_{j}"),
                "others"
            ))

            markers <- FindMarkers(
                so,
                ident.1 = glue::glue("{id}_{i}_integrated_{j}"),
                ident.2 = "others",
                min.pct = 0.25,
                logfc.threshold = 0.25,
                verbose = FALSE,
                only.pos = TRUE
            ) %>%
                filter(.data$p_val_adj < 0.01) %>%
                arrange(desc(.data$avg_log2FC))

            write.csv(
                markers,
                file.path(marker_folder, glue::glue("{id}_{i}_integrated_{j}_markers.csv"))
            )
        }

        group_1 <- rownames(umap_df)[which(umap_df$dataset == id & umap_df$mb == glue::glue("{j}"))]

        if (length(group_1) < 3) {
            next
        }

        Idents(so) <- factor(ifelse(
            colnames(so) %in% group_1,
            glue::glue("{id}_integrated_{j}"),
            "others"
        ))

        markers <- FindMarkers(
            so,
            ident.1 = glue::glue("{id}_integrated_{j}"),
            ident.2 = "others",
            min.pct = 0.25,
            logfc.threshold = 0.25,
            verbose = FALSE,
            only.pos = TRUE
        ) %>%
            filter(.data$p_val_adj < 0.01) %>%
            arrange(desc(.data$avg_log2FC))

        write.csv(
            markers,
            file.path(marker_folder, glue::glue("{id}_integrated_{j}_markers.csv"))
        )
    }
}


# considering that the parent directory is the root of the repository
output_path <- "output/integration"

# ============== CASE 1 ==============
# ORGANOID 8w vs RAWLINGS 6w

organoid <- readRDS(file.path(organoid_object_folder, "aggregated/HV/by_organ_age_condition/aggr_filtered_coding_highrp-organ-Lung-age-08w-condition--Wnt_without_oesophagus_mesenchymal.rds"))
rawlings <- readRDS(file.path(rawlings_object_folder, "data/R_objects/6w_epithelial.rds"))
# r6: 6, 8, 15
chosen_k_list <- c(6, 8, 14, 15)
chosen_k <- 50
id <- "r6w_vs_o8w"

# ORGANOID 12w vs RAWLINGS 8w

organoid <- readRDS(file.path(organod_object_folder, "aggregated/HV/by_organ_age_condition/aggr_filtered_coding_highrp-organ-Lung-age-12w-condition--Wnt_without_oesophagus_mesenchymal_low_ncount_cluster.rds"))
rawlings <- readRDS(file.path(rawlings_object_folder, "data/R_objects/8w_epithelial.rds"))
# r8: 3, 5, 7, 9, 10, 11, 12, 13, 16, 17, 18, 19, 21, 22, 23, 24, 25 (debatable from 22)
chosen_k_list <- c(3, 5, 7, 9, 10, 11, 12, 13, 15)
chosen_k <- 70
id <- "r8w_vs_o12w"

################

dts_names <- stringr::str_split(id, "_vs_")[[1]]
umap_df <- readRDS(file.path(output_path, "integration_quality", id, "umaps", glue::glue("k_{chosen_k}_umap_df.rds")))
ca_clustering <- readRDS(file.path(output_path, "integration_quality", id, "clustassess", glue::glue("k_{chosen_k}_graph_clustering_assessment.rds")))
ca_clustering$split_by_k$Louvain <- NULL
ca_clustering$split_by_k$Louvain.refined <- NULL
ca_clustering$split_by_resolution$Louvain <- NULL
ca_clustering$split_by_resolution$Louvain.refined <- NULL

# the following helped us in determining the optimal number of clusters
nclusters <- names(ca_clustering$split_by_k$SLM)

freqs <- sapply(
    nclusters,
    function(x) {
        sum(sapply(ca_clustering$split_by_k$SLM[[x]]$partitions, function(y) {
            y$freq
        }))
    }
)

print(freqs[freqs >= 30])

avg_ecc <- sapply(
    nclusters,
    function(x) {
        mean(ca_clustering$split_by_k$SLM[[x]]$ecc)
    }
)
print(avg_ecc[freqs >= 30])
ClustAssess::plot_k_n_partitions(ca_clustering)

# find the markers
for (chosen_k in as.character(chosen_k_list)) {
    print(paste(id, chosen_k))
    umap_df$mb <- factor(ca_clustering$split_by_k$SLM[[chosen_k]]$partitions[[1]]$mb)
    find_markers(rawlings, umap_df, chosen_k, dts_names[1], file.path(output_path, "markers", id))
    find_markers(organoid, umap_df, chosen_k, dts_names[2], file.path(output_path, "markers", id))
}
