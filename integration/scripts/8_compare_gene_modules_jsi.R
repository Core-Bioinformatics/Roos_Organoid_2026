library(Seurat)
library(ggplot2)
library(gprofiler2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(grid)

get_jsi <- function(clustering1, clustering2) {
    n1 <- length(clustering1)
    n2 <- length(clustering2)

    jaccard <- matrix(0, nrow = n1, ncol = n2)
    for (i in seq_len(n1)) {
        for (j in seq_len(n2)) {
            jaccard[i, j] <- length(intersect(clustering1[[i]], clustering2[[j]])) / length(union(clustering1[[i]], clustering2[[j]]))
        }
    }

    jaccard
}

generate_jsi_heatmap <- function(df1, df2, name_first, name_second, htmp_title = "", gene_pool = NULL) {
    if (!is.null(gene_pool)) {
        df1 <- df1[gene_pool, ]
        df2 <- df2[gene_pool, ]
    }

    df1$module <- as.factor(df1$module)
    clustering1_list <- split(rownames(df1), df1$module)

    df2$module <- as.factor(df2$module)
    clustering2_list <- split(rownames(df2), df2$module)

    jsi <- get_jsi(clustering1_list, clustering2_list)
    rownames(jsi) <- paste0(name_first, "_", seq_len(nrow(jsi)))
    colnames(jsi) <- paste0(name_second, "_", seq_len(ncol(jsi)))

    Heatmap(
        jsi,
        name = "JSI",
        col = viridis::viridis(nrow(jsi) * ncol(jsi) * 2), 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", jsi[i, j]), x, y, just = "center", gp = gpar(fontsize = 20))
        },
        column_title = htmp_title,
        column_names_rot = 0
    )
}


######## Organoid 08w vs Rawlings 06w ########
rawling_psd_path <- file.path(rawlings_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "6w_epithelial", "no_loop", "Most_Abundant", "1500", "integration_trajectory_start")
organoid_psd_path <- file.path(sc_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "by_organ_age_condition", "organ-Lung-age-08w-condition--Wnt_without_oesophagus_mesenchymal", "no_loop", "Most_Abundant", "4000", "integration_trajectory_start")

output_path <- file.path(rawlings_project, "output", "comparative_analysis", "gene_modules")
jsi_path <- file.path(output_path, "JSI", "r6w_vs_o8w")
if (!dir.exists(jsi_path)) {
    dir.create(jsi_path, recursive = TRUE)
}

rawling_psd_files <- list.files(file.path(rawling_psd_path, "clusters_moran_01"))
organoid_psd_files <- list.files(file.path(organoid_psd_path, "clusters_moran_01"))

options <- list(
    all_genes = NULL,
    common_genes = intersect(
        rownames(read.csv(file.path(rawling_psd_path, "clusters_moran_01", rawling_psd_files[[1]]), sep = ",", header = TRUE, row.names = 1)),
        rownames(read.csv(file.path(organoid_psd_path, "clusters_moran_01", organoid_psd_files[[1]]), sep = ",", header = TRUE, row.names = 1))
    )
)

for (op in names(options)) {
    pdf(file.path(jsi_path, paste0(op, ".pdf")), width = 15, height = 15)
    for (rclust in list.files(file.path(rawling_psd_path, "clusters_moran_01"))) {
        r_nclust <- strsplit(strsplit(rclust, "\\.")[[1]][1], "_")[[1]][2]

        r_modules <- read.csv(file.path(rawling_psd_path, "clusters_moran_01", rclust), sep = ",", header = TRUE, row.names = 1)

        for (oclust in list.files(file.path(organoid_psd_path, "clusters_moran_01"))) {
            o_nclust <- strsplit(strsplit(oclust, "\\.")[[1]][1], "_")[[1]][2]

            o_modules <- read.csv(file.path(organoid_psd_path, "clusters_moran_01", oclust), sep = ",", header = TRUE, row.names = 1)

            print(generate_jsi_heatmap(
                o_modules,
                r_modules,
                "o8w",
                "r6w",
                glue::glue("{op} JSI: O8W - {o_nclust} vs R6W - {r_nclust}"),
                gene_pool = options[[op]]))
        }
    }
    dev.off()
}


######## Organoid 12w vs Rawlings 08w ########
rawling_psd_path <- file.path(rawlings_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "8w_epithelial", "no_loop", "Highly_Variable", "3000", "integration_trajectory_start")
organoid_psd_path <- file.path(sc_project, "output", "pseudotime_shiny_apps", "second_part_analysis", "by_organ_age_condition", "organ-Lung-age-12w-condition--Wnt_without_oesophagus_mesenchymal_low_ncount_cluster", "no_loop", "Highly_Variable", "4500", "integration_trajectory_start")

jsi_path <- file.path(output_path, "JSI", "r8w_vs_o12w")
if (!dir.exists(jsi_path)) {
    dir.create(jsi_path, recursive = TRUE)
}

rawling_psd_files <- list.files(file.path(rawling_psd_path, "clusters_moran_01"))
organoid_psd_files <- list.files(file.path(organoid_psd_path, "clusters_moran_01"))

options <- list(
    all_genes = NULL,
    common_genes = intersect(
        rownames(read.csv(file.path(rawling_psd_path, "clusters_moran_01", rawling_psd_files[[1]]), sep = ",", header = TRUE, row.names = 1)),
        rownames(read.csv(file.path(organoid_psd_path, "clusters_moran_01", organoid_psd_files[[1]]), sep = ",", header = TRUE, row.names = 1))
    )
)

for (op in names(options)) {
    pdf(file.path(jsi_path, paste0(op, ".pdf")), width = 15, height = 15)
    for (rclust in list.files(file.path(rawling_psd_path, "clusters_moran_01"))) {
        r_nclust <- strsplit(strsplit(rclust, "\\.")[[1]][1], "_")[[1]][2]

        r_modules <- read.csv(file.path(rawling_psd_path, "clusters_moran_01", rclust), sep = ",", header = TRUE, row.names = 1)

        for (oclust in list.files(file.path(organoid_psd_path, "clusters_moran_01"))) {
            o_nclust <- strsplit(strsplit(oclust, "\\.")[[1]][1], "_")[[1]][2]

            o_modules <- read.csv(file.path(organoid_psd_path, "clusters_moran_01", oclust), sep = ",", header = TRUE, row.names = 1)

            print(generate_jsi_heatmap(
                o_modules,
                r_modules,
                "o12w",
                "r8w",
                glue::glue("{op} JSI: O12W - {o_nclust} vs R8W - {r_nclust}"),
                gene_pool = options[[op]]))
        }
    }
    dev.off()
}
