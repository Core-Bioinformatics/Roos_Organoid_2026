qlibrary(Seurat)
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
    two_distinc_values <- min(jsi) != max(jsi)
    if (two_distinc_values) {
        col_values <- viridis::viridis(nrow(jsi) * ncol(jsi) * 2)
    } else {
        col_values <- "white"
    }

    Heatmap(
        jsi,
        name = "JSI",
        col = col_values, #colorRamp2(c(0, 1), c("white", "blue")),
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


primary_psd_app <- "pseudotime_shiny_apps/second_part_analysis/primary_liver_fetal_hepatoblasts"
organoid_psd_app <- "pseudotime_shiny_apps/second_part_analysis/organoid_liver_hepatoblasts"

comparison_pairs <- list(
    c("02_6", "02_4")
)

for (comp_pair in comparison_pairs) {
    output_path <- file.path("output", "integration_primary", glue::glue("comparison_analysis_organoid_{comp_pair[2]}_vs_primary_{comp_pair[1]}"), "enrichment")
    if (!dir.exists(output_path)) {
        dir.create(output_path, recursive = TRUE)
    }

    primary_enr_files <- list.files(file.path(primary_psd_app, "enrichment_all", comp_pair[1]))
    organoid_enr_files <- list.files(file.path(organoid_psd_app, "enrichment_all", comp_pair[2]))
    print(paste("primary ", comp_pair[1], " vs organoid ", comp_pair[2]))

    for (source in sources) {
        print(source)
        pdf(file.path(output_path, paste0(stringr::str_replace_all(source, ":", ""), ".pdf")), width = 15, height = 15)
        organoid_clustering_list <- get_clustering_list(file.path(organoid_psd_app, "enrichment_all", comp_pair[2]), source)
        if (length(organoid_clustering_list) == 0) {
            next
        }

        primary_clustering_list <- get_clustering_list(file.path(primary_psd_app, "enrichment_all", comp_pair[1]), source)

        if (length(primary_clustering_list) == 0) {
            next
        }

        print(generate_jsi_heatmap(
            primary_clustering_list,
            organoid_clustering_list,
            "primary",
            "organoid",
            glue::glue("Enriched terms {source} JSI: primary - {length(primary_clustering_list)} vs organoid - {length(organoid_clustering_list)}")))
        dev.off()
    }
}

