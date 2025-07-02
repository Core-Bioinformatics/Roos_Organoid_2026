library(monocle3)
library(dplyr)
library(reshape2)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(gplots)
library(qualpalr)
library(ComplexHeatmap)
library(circlize)

# the way the clusters of genes are picked
# on the psd second shiny app, filter the genes by moran's I test (0.2, 0.15, .. 0)
# the n cluster the genes using 50 reps with res 0.1 -> 5
# filter out the clusters with freq < 30 or ecc < 0.9
# pick the highest #clusters
# find an ordering based on that

# pseudotime_app_path <- file.path(project_dir, "output/pseudotime_apps/second_part_analysis/SITTF7_ALI/no_loop/Highly_Variable/1000/basal_cells_to_secretory_cells")
min_avg_expr <- 0
k_weights <- 30 
# apply_conv <- TRUE
apply_conv <- FALSE 
cap_value <- 3
percentage_cutoff <- 1.01
value_cutoff <- 0.25

mon_obj <- readRDS(file.path(pseudotime_app_path, "monocle_object.rds"))
ordering_df <- read.csv(file.path(pseudotime_app_path, "cluster_orders.csv"), sep = ";", comment.char = "#")

for (apply_conv in c(TRUE, FALSE)) {
    for (i in seq_len(nrow(ordering_df))) {
        if (ordering_df[["Cluster_order"]][i] == "") {
            next
        }
        print(ordering_df[i, ])

        gene_clusters_order <- ordering_df[["Cluster_order"]][i]

        if (!is.integer(gene_clusters_order)) {
            gene_clusters_order <- as.integer(strsplit(gene_clusters_order, ",")[[1]])
        } else {
            gene_clusters_order <- c(gene_clusters_order)
        }

        prefix <- glue::glue("moran_{ordering_df[['Moran_thresh']][i]}_cluster_{ordering_df[['Nclusters']][i]}")

        gene_modules <- read.csv(file.path(pseudotime_app_path, glue::glue("{prefix}.csv"))) %>% filter(.data$avg_expression >= min_avg_expr)

        genes_in_order <- c()
        gene_clusters <- c()

        for (cl in gene_clusters_order) {
            gene_to_add <- (gene_modules %>% filter(.data$module == cl))$X
            genes_in_order <- c(genes_in_order, gene_to_add)
            gene_clusters <- c(gene_clusters, rep(cl, length(gene_to_add)))
        }

        cells_in_order <- colnames(mon_obj)[order(pseudotime(mon_obj))]

        submatrix <- as.matrix(exprs(mon_obj)[genes_in_order, cells_in_order])
        submatrix[submatrix > cap_value] <- cap_value
        perc_expressed <- rowSums(submatrix > value_cutoff) / ncol(submatrix)
        submatrix <- submatrix[perc_expressed < percentage_cutoff, ]
        genes_in_order <- genes_in_order[perc_expressed < percentage_cutoff]
        gene_clusters <- gene_clusters[perc_expressed < percentage_cutoff]

        gene_clusters <- factor(gene_clusters)
        nclusters <- length(gene_clusters_order)

        weights <- rep(1 / k_weights, k_weights)

        if (apply_conv) {

            submatrix_convolved <- matrix(0, nrow = nrow(submatrix), ncol = ncol(submatrix) + k_weights - 1)
            for (i in seq_len(nrow(submatrix))) {
                submatrix_convolved[i, ] <- convolve(submatrix[i, ], weights, type = "open")
            }
            prefix <- glue::glue("{prefix}_conv_k{k_weights}")
            submatrix_convolved <- submatrix_convolved[, -seq_len(k_weights / 2)]
            submatrix_convolved <- submatrix_convolved[, seq_len(ncol(mon_obj))]
        } else {
            submatrix_convolved <- submatrix
            prefix <- glue::glue("{prefix}_no_conv")
        }

        rownames(submatrix_convolved) <- genes_in_order
        
        if (nclusters > 1) {
            colors_rows <- qualpalr::qualpal(nclusters)$hex
        } else {
            colors_rows <- "#072c07"
        }

        # NEW METHOD
        st_target <- c(min(submatrix_convolved ), (max(submatrix_convolved) + min(submatrix_convolved)) / 2, max(submatrix_convolved))

        cols_gene_clusters <- qualpal(n = nlevels(gene_clusters), colorspace = "pretty_dark")$hex
        names(cols_gene_clusters) <- as.character(levels(gene_clusters))
        nsamples <- length(unique(colData(mon_obj)$barcode))
        if (nsamples > 1) {
            cols_barcode <- qualpal(n = nsamples, colorspace = "pretty")$hex
        } else {
            cols_barcode <- "#072c07"
        }
        names(cols_barcode) <- unique(colData(mon_obj)$barcode)

        bottom_ha <- HeatmapAnnotation(
            samples = colData(mon_obj)$barcode,
            pseudotime = pseudotime(mon_obj)[cells_in_order],
            col = list(
                samples = cols_barcode,
                pseudotime = colorRamp2(c(0, max(pseudotime(mon_obj)[cells_in_order])), c("black", "#288deb"))
            ),
            show_annotation_name = TRUE,
            annotation_name_side = "left"
        )

        left_ha <- rowAnnotation(
            gene_clusters = gene_clusters,
            col = list(gene_clusters = cols_gene_clusters),
            show_annotation_name = FALSE,
            show_legend = TRUE
        )

        unique_clusters <- unique(gene_clusters)
        
        remapped <- as.character(sapply(seq_along(unique_clusters), function(i) { sprintf("%.2d", i) }))
        names(remapped) <- as.character(unique_clusters)
        args <- c(list(as.character(gene_clusters)), remapped)
        remapped <- do.call(recode, args)

        htmp <- Heatmap(
            submatrix_convolved,
            name = "expression level",
            row_order = seq_len(nrow(submatrix_convolved)),
            column_order = seq_len(ncol(submatrix_convolved)),
            show_column_names = FALSE,
            show_row_names = FALSE,
            use_raster = TRUE,
            raster_by_magick = TRUE,
            col = colorRamp2(st_target, viridis(3)),
            row_title = unique_clusters,
            row_split = remapped,# gene_clusters,
            # split = as.character(gene_clusters),
            # cluster_rows = FALSE,
            # cluster_row_slices = FALSE,
            # raster_magick_filter = "Gaussian",
            raster_quality = 1,
            bottom_annotation = bottom_ha,
            left_annotation = left_ha,
            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(direction = "horizontal", legend_width = unit(5, "cm"))

        )

        pdf(file.path(pseudotime_app_path, glue::glue("{prefix}_heatmap.pdf")), width = 15, height = 15) #0.14 * length(genes_in_order))
        draw(
            htmp,
            heatmap_legend_side = "top",
            annotation_legend_side = "right",
            legend_grouping = "original"
        )
        dev.off()
    }
}
