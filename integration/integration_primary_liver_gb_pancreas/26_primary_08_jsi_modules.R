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
        col = circlize::colorRamp2(c(0, 0.5, 1), viridis::viridis(3)),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", jsi[i, j]), x, y, just = "center", gp = gpar(fontsize = 20))
        },
        column_title = htmp_title,
        column_names_rot = 0
    )
}

primary_app_path <- "pseudotime_shiny_apps/second_part_analysis/primary_liver_fetal_hepatoblasts"
organoid_app_path <- "pseudotime_shiny_apps/second_part_analysis/organoid_liver_hepatoblasts"

comparison_pairs <- list(
    c("02_6", "02_4")
)

for (comp_pair in comparison_pairs) {
    primary_module_path <- file.path(primary_app_path, paste0("primary_moran_", comp_pair[1], "_modules.csv"))
    organoid_module_path <- file.path(organoid_app_path, paste0("organoid_moran_", comp_pair[2], "_modules.csv"))

    jsi_path <- file.path("/servers/sutherland-scratch/andi/projects/0_2208_organoid/output", "integration_primary", glue::glue("comparison_analysis_organoid_{comp_pair[2]}_vs_primary_{comp_pair[1]}"), "gene_modules")
    if (!dir.exists(jsi_path)) {
        dir.create(jsi_path, recursive = TRUE)
    }

    options <- list(
        all_genes = NULL,
        common_genes = intersect(
            rownames(read.csv(primary_module_path, sep = ",", header = TRUE, row.names = 1)),
            rownames(read.csv(organoid_module_path, sep = ",", header = TRUE, row.names = 1))
        )
    )

    for (op in names(options)) {
        pdf(file.path(jsi_path, paste0(op, ".pdf")), width = 15, height = 15)

        primary_modules <- read.csv(primary_module_path, sep = ",", header = TRUE, row.names = 1)
        primary_nclust <- length(unique(primary_modules$module))

        organoid_modules <- read.csv(organoid_module_path, sep = ",", header = TRUE, row.names = 1)
        organoid_nclust <- length(unique(organoid_modules$module))


        print(generate_jsi_heatmap(
            primary_modules,
            organoid_modules,
            "primary",
            "organoid",
            glue::glue("{op} JSI: primary - {primary_nclust} vs organoid - {organoid_nclust}"),
            gene_pool = options[[op]])
        )
        dev.off()
    }
}

