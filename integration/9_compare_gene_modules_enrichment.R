library(Seurat)
library(ggplot2)
library(gprofiler2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(grid)

sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")

get_jsi <- function(clustering1, clustering2) {
    n1 <- length(clustering1)
    n2 <- length(clustering2)

    jaccard <- matrix(0, nrow = n1, ncol = n2)
    for (i in seq_len(n1)) {
        for (j in seq_len(n2)) {
            jaccard[i, j] <- length(intersect(clustering1[[i]], clustering2[[j]])) / length(union(clustering1[[i]], clustering2[[j]]))
        }
    }
    rownames(jaccard) <- names(clustering1)
    colnames(jaccard) <- names(clustering2)

    jaccard
}

generate_jsi_heatmap <- function(clustering1_list, clustering2_list, name_first, name_second, htmp_title = "") {
    jsi <- get_jsi(clustering1_list, clustering2_list)
    rownames(jsi) <- paste0(name_first, "_", rownames(jsi))
    colnames(jsi) <- paste0(name_second, "_", colnames(jsi))

    Heatmap(
        jsi,
        name = "JSI",
        col = viridis::viridis(nrow(jsi) * ncol(jsi) * 2), #colorRamp2(c(0, 1), c("white", "blue")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", jsi[i, j]), x, y, just = "center", gp = gpar(fontsize = 20))
        },
        column_title = htmp_title,
        column_names_rot = 0
    )
}

get_clustering_list <- function(enrichment_folder, source_filter) {
    enr_files <- list.files(enrichment_folder)

    clustering_list <- list()

    for (enr_file in enr_files) {
        enr <- read.csv(file.path(enrichment_folder, enr_file), sep = ",", header = TRUE, row.names = 1)
        enr <- enr %>% filter(source == source_filter)
        enr <- unique(enr$term_id)

        index_cluster <- strsplit(enr_file, "_")[[1]][2]

        if (length(enr) > 0) {
            clustering_list[[index_cluster]] <- enr
        }

    }
    clustering_list
}


######## Organoid 08w vs Rawlings 06w ########
rawling_enr_path <- file.path(rawlings_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "6w_epithelial", "no_loop", "Most_Abundant", "1500", "integration_trajectory_start")
organoid_enr_path <- file.path(sc_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "by_organ_age_condition", "organ-Lung-age-08w-condition--Wnt_without_oesophagus_mesenchymal", "no_loop", "Most_Abundant", "4000", "integration_trajectory_start")

output_path <- file.path(rawlings_project, "output", "comparative_analysis", "gene_modules")
enrich_path <- file.path(output_path, "enrichment_all", "r6w_vs_o8w")
if (!dir.exists(enrich_path)) {
    dir.create(enrich_path, recursive = TRUE)
}

rawling_enr_files <- list.files(file.path(rawling_enr_path, "clusters_moran_01", "enrichment_all"))
organoid_enr_files <- list.files(file.path(organoid_enr_path, "clusters_moran_01", "enrichment_all"))

for (source in sources) {
    print(source)
    pdf(file.path(enrich_path, paste0(stringr::str_replace_all(source, ":", ""), ".pdf")), width = 15, height = 15)
    for (rclust in rawling_enr_files) {
        rclustering_list <- get_clustering_list(file.path(rawling_enr_path, "clusters_moran_01", "enrichment_all", rclust), source)
        if (length(rclustering_list) == 0) {
            next
        }

        for (oclust in organoid_enr_files) {
            oclustering_list <- get_clustering_list(file.path(organoid_enr_path, "clusters_moran_01", "enrichment_all", oclust), source)

            if (length(oclustering_list) == 0) {
                next
            }

            print(generate_jsi_heatmap(
                oclustering_list,
                rclustering_list,
                "o8w",
                "r6w",
                glue::glue("Enriched terms {source} JSI: O8W - {oclust} vs R6W - {rclust}")))
        }
    }
    dev.off()
}


######## Organoid 12w vs Rawlings 08w ########
rawling_enr_path <- file.path(rawlings_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "8w_epithelial", "no_loop", "Highly_Variable", "3000", "integration_trajectory_start")
organoid_enr_path <- file.path(sc_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "by_organ_age_condition", "organ-Lung-age-12w-condition--Wnt_without_oesophagus_mesenchymal_low_ncount_cluster", "no_loop", "Highly_Variable", "4500", "integration_trajectory_start")

enrich_path <- file.path(output_path, "enrichment_all", "r8w_vs_o12w")
if (!dir.exists(enrich_path)) {
    dir.create(enrich_path, recursive = TRUE)
}

rawling_enr_files <- list.files(file.path(rawling_enr_path, "clusters_moran_01", "enrichment_all"))
organoid_enr_files <- list.files(file.path(organoid_enr_path, "clusters_moran_01", "enrichment_all"))

for (source in sources) {
    print(source)
    pdf(file.path(enrich_path, paste0(stringr::str_replace_all(source, ":", ""), ".pdf")), width = 15, height = 15)
    for (rclust in rawling_enr_files) {
        rclustering_list <- get_clustering_list(file.path(rawling_enr_path, "clusters_moran_01", "enrichment_all", rclust), source)
        if (length(rclustering_list) == 0) {
            next
        }

        for (oclust in organoid_enr_files) {
            oclustering_list <- get_clustering_list(file.path(organoid_enr_path, "clusters_moran_01", "enrichment_all", oclust), source)

            if (length(oclustering_list) == 0) {
                next
            }

            print(generate_jsi_heatmap(
                oclustering_list,
                rclustering_list,
                "o12w",
                "r8w",
                glue::glue("Enriched terms {source} JSI: O12W - {oclust} vs R8W - {rclust}")))
        }
    }
    dev.off()
}

