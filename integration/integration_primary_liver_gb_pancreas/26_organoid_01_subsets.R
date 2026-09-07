library(Seurat)

project_folder <- ""
objects_folder <- file.path(project_folder, "R_objects")

organoid_used_clusters <- list(
    "Liver" = list(
        ftype = "Highly_Variable",
        fsize = 4500,
        clmethod = "SLM",
        nclusters = 13,
        selected_clusters = c(4, 6, 11),
        new_name = "organoid_liver_hepatoblasts"
    )
)

for (organ in names(organoid_used_clusters)) {
    
    so <- readRDS(file.path(objects_folder, "aggregated", "HV", "by_organ", paste0("aggr_filtered_coding_highrp-organ-", organ, ".rds")))
    if (!is.null(so) && length(organoid_used_clusters[[organ]]) > 0) {
        ca_object <- readRDS(file.path(objects_folder, "clustassess", "by_organ", paste0("organ-", organ, ".rds")))
        used_partition <- ClustAssess::get_clusters_from_clustassess_object(
            clustassess_object = ca_object,
            feature_type = organoid_used_clusters[[organ]]$ftype,
            feature_size = organoid_used_clusters[[organ]]$fsize,
            clustering_method = organoid_used_clusters[[organ]]$clmethod,
            nclusters = organoid_used_clusters[[organ]]$nclusters
        )[[1]]$partitions[[1]]$mb
        so <- subset(so, cells = colnames(so)[used_partition %in% organoid_used_clusters[[organ]]$selected_clusters])
    }

    DefaultAssay(so) <- "RNA"
    so@assays$SCT <- NULL
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 4000)
    so <- ScaleData(so, features = rownames(so), verbose = FALSE)
    so <- CellCycleScoring(so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
    so <- RunPCA(so)
    so <- RunUMAP(so, dims = 1:30, reduction = "pca")

    qs::qsave(so, file.path(objects_folder, paste0(organoid_used_clusters[[organ]]$new_name, "_filtered_normalized.qs")), nthreads = 30)
}


