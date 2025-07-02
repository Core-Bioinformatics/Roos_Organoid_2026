library(gprofiler2)

sources <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")
organism <- "hsapiens"

gene_paths <- c(
    file.path(rawlings_psd_folder, "6w_epithelial/no_loop/Most_Abundant/1500/integration_trajectory_start/clusters_moran_01"),
    file.path(rawlings_psd_folder, "8w_epithelial/no_loop/Highly_Variable/3000/integration_trajectory_start/clusters_moran_01"),
    file.path(organoid_psd_folder, "by_organ_age_condition/organ-Lung-age-12w-condition--Wnt_without_oesophagus_mesenchymal_low_ncount_cluster/no_loop/Highly_Variable/4500/integration_trajectory_start/clusters_moran_01"),
    file.path(organoid_psd_folder, "by_organ_age_condition/organ-Lung-age-08w-condition--Wnt_without_oesophagus_mesenchymal/no_loop/Most_Abundant/4000/integration_trajectory_start/clusters_moran_01")
)

for (gene_path in gene_paths) {
    psd_files <- list.files(gene_path)
    psd_files <- psd_files[grepl("csv", psd_files)]

    enrichment_path <- file.path(gene_path, "enrichment_all")
    if (!dir.exists(enrichment_path)) {
        dir.create(enrichment_path, recursive = TRUE)
    }

    gene_pool <- as.character(rhdf5::h5read(file.path(dirname(gene_path), "expression.h5"), "genes"))

    for (psd_file in psd_files) {
        print(paste(gene_path, psd_file))
        nclust <- strsplit(strsplit(psd_file, "\\.")[[1]][1], "_")[[1]][2]
        modules <- read.csv(file.path(gene_path, psd_file), sep = ",", header = TRUE, row.names = 1)
        # gene_pool <- rownames(modules)
        modules$module <- as.factor(modules$module)
        genes_from_modules <- rownames(modules)
        clustering_list <- split(genes_from_modules, modules$module)

        if (!dir.exists(file.path(enrichment_path, nclust))) {
            dir.create(file.path(enrichment_path, nclust), recursive = TRUE)
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

                write.csv(gprof_res, file.path(enrichment_path, nclust, paste0("module_", sprintf("%02d", i), "_enrichment.csv")))
            }
        }
    }
}
