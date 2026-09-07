library(gprofiler2)

sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")
organism <- "hsapiens"

app_paths <- c(
    "pseudotime_shiny_apps/second_part_analysis/primary_liver_fetal_hepatoblasts"
)

gene_paths <- unlist(lapply(app_paths, function(app_path) {
    list.files(app_path, pattern = "modules.csv", full.names = TRUE)
}))


for (gene_path in gene_paths) {
    app_dir <- dirname(gene_path)
    cluster_details <- strsplit(basename(gene_path), "moran_")[[1]][2]
    cluster_details <- strsplit(cluster_details, "_modules.csv")[[1]][1]
    print(paste("Processing", cluster_details))

    enrichment_path <- file.path(app_dir, "enrichment_all", cluster_details)
    if (!dir.exists(enrichment_path)) {
        dir.create(enrichment_path, recursive = TRUE)
    }

    gene_pool <- as.character(rhdf5::h5read(file.path(app_dir, "expression.h5"), "genes"))

    
    modules <- read.csv(gene_path, sep = ",", header = TRUE, row.names = 1)
    # gene_pool <- rownames(modules)
    modules$module <- as.factor(modules$module)
    nclust <- length(unique(modules$module))
    genes_from_modules <- rownames(modules)
    clustering_list <- split(genes_from_modules, modules$module)

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

            write.csv(gprof_res, file.path(enrichment_path, paste0("module_", sprintf("%02d", i), "_enrichment.csv")))
        }
    }
}
