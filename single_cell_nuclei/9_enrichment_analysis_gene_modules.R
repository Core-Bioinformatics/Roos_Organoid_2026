library(gprofiler2)

sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")
organism <- "hsapiens"

psd_folder <- "output/pseudotime_apps/second_part_analysis"

gene_paths <- file.path(psd_folder, c(
    file.path("aggr_filtered_without_mesenchymal/no_loop/Highly_Variable/1250", c(
        "basal_cells_to_ciliated_cells",
        "organoid_to_ciliated_cells",
        "organoid_to_ciliated_cells_2"
    )),
    file.path("SITTF7_ALI/no_loop/Highly_Variable/1000", c(
        "basal_cells_to_secretory_cells",
        "basal_cells_to_secretory_cells_2"
    ))
), "pseudotime_heatmaps")

for (gene_path in gene_paths) {
    psd_files <- list.files(gene_path)
    psd_files <- psd_files[grepl("csv", psd_files)]

    enrichment_path <- file.path(gene_path, "enrichment")
    if (!dir.exists(enrichment_path)) {
        dir.create(enrichment_path, recursive = TRUE)
    }

    gene_pool <- as.character(rhdf5::h5read(file.path(dirname(gene_path), "expression.h5"), "genes"))

    for (psd_file in psd_files) {
        if (psd_file == "cluster_orders.csv") {
            next
        }
        print(paste(gene_path, psd_file))
        target_folder <- file.path(enrichment_path, strsplit(psd_file, "\\.csv")[[1]][1])
        modules <- read.csv(file.path(gene_path, psd_file), sep = ",", header = TRUE, row.names = 1)
        # gene_pool <- rownames(modules)
        modules$module <- as.factor(modules$module)
        genes_from_modules <- rownames(modules)
        clustering_list <- split(genes_from_modules, modules$module)

        if (!dir.exists(file.path(target_folder))) {
            dir.create(file.path(target_folder), recursive = TRUE)
        }

        for (i in seq_along(clustering_list)) {
            genes <- clustering_list[[i]]
            print(paste(length(gene_pool), length(genes)))
            gprof_res <- gprofiler2::gost(
                query = genes,
                sources = sources,
                organism = organism,
                evcodes = TRUE,
                domain_scope = "custom",
                custom_bg = gene_pool
            )

            if (!is.null(gprof_res)) {
                gprof_res$result$parents <- sapply(gprof_res$result$parents, toString)
                gprof_res <- gprof_res$result
                # uncomment this if you want to remove the evcodes and the genes participating to the intersection
                # gprof_res <- gprof_res[ , seq_len(ncol(gprof_res) - 2)]

                write.csv(gprof_res, file.path(target_folder, paste0("module_", sprintf("%02d", i), "_enrichment.csv")))
            }
        }
    }
}
