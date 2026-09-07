library(Seurat)
library(harmony)
library(UpSetR)

project_folder <- ""
objects_folder <- file.path(project_folder, "R_objects")

target_organs <- list(
    # "Intestine" = "",
    "Pancreas" = c("Ductal cells"),
    "Liver" = c("Hepatoblasts")
)

organoid_used_clusters <- list(
    "Liver" = list(
        ftype = "Highly_Variable",
        fsize = 4500,
        clmethod = "SLM",
        nclusters = 13,
        selected_clusters = c(4, 6, 11)
    ),
    "Pancreas" = list()
)
cao_objects <- list()
organoid_objects <- list()

for (organ in names(target_organs)) {
    target_subtype <- target_organs[[organ]]
    new_organ <- organ
    if (target_subtype != "" && !is.null(target_subtype)) {
        new_organ <- paste(organ, target_subtype, sep = "_")
    }
    print(paste(organ, new_organ))
    cao_objects[[organ]] <- qs::qread(file.path(objects_folder, paste0("cao_", new_organ, "_filtered_normalized.qs")), nthreads = 30)
    if (length(target_subtype) > 0 && target_subtype != "") {
        cao_objects[[organ]] <- subset(cao_objects[[organ]], subset = Main_cluster_name %in% target_subtype)
    }
    organoid_objects[[organ]] <- readRDS(file.path(objects_folder, "aggregated", "HV", "by_organ", paste0("aggr_filtered_coding_highrp-organ-", organ, ".rds")))
    if (!is.null(organoid_used_clusters[[organ]]) && length(organoid_used_clusters[[organ]]) > 0) {
        ca_object <- readRDS(file.path(objects_folder, "clustassess", "by_organ", paste0("organ-", organ, ".rds")))
        used_partition <- ClustAssess::get_clusters_from_clustassess_object(
            clustassess_object = ca_object,
            feature_type = organoid_used_clusters[[organ]]$ftype,
            feature_size = organoid_used_clusters[[organ]]$fsize,
            clustering_method = organoid_used_clusters[[organ]]$clmethod,
            nclusters = organoid_used_clusters[[organ]]$nclusters
        )[[1]]$partitions[[1]]$mb
        organoid_objects[[organ]] <- subset(organoid_objects[[organ]], cells = colnames(organoid_objects[[organ]])[used_partition %in% organoid_used_clusters[[organ]]$selected_clusters])
    }

    genes_expressed_in_cao <- rownames(cao_objects[[organ]])[which(Matrix::rowSums(GetAssayData(cao_objects[[organ]], assay = "RNA", layer = "counts") > 0) > 0)]
    genes_expressed_in_organoid <- rownames(organoid_objects[[organ]])[which(Matrix::rowSums(GetAssayData(organoid_objects[[organ]], assay = "RNA", layer = "counts") > 0) > 0)]

    # upset plot between the two sets of genes
    gene_list <- list(
        "Cao" = genes_expressed_in_cao,
        "Organoid" = genes_expressed_in_organoid
    )
    print(upset(fromList(gene_list), order.by = "freq", nsets = 2, nintersects = NA, sets.bar.color = c("#E41A1C", "#377EB8"), main.bar.color = "#4DAF4A", text.scale = c(2, 2, 2, 1.5, 2, 1.5)))

    mtd <- cao_objects[[organ]]@meta.data
    mtd <- mtd[, c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp", "Organ", "Development_day", "Fetus_id")]
    mtd$project <- "Cao"
    colnames(mtd) <- c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp", "organ", "age", "sample_name", "project")
    cao_objects[[organ]]@meta.data <- mtd
    cao_objects[[organ]]@assays$RNA@layers$scale.data <- NULL



    mtd <- organoid_objects[[organ]]@meta.data
    mtd <- mtd[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "organ", "age", "sample_name")]
    mtd$project <- "Organoid"
    colnames(mtd) <- c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp", "organ", "age", "sample_name", "project")
    organoid_objects[[organ]]@meta.data <- mtd
    # remove the SCT assay
    DefaultAssay(organoid_objects[[organ]]) <- "RNA"
    organoid_objects[[organ]][["SCT"]] <- NULL

    integrated_so <- merge(cao_objects[[organ]], organoid_objects[[organ]], add.cell.ids = c("Cao", "Organoid"))
    # integrated_so <- JoinLayers(integrated_so)
    integrated_so <- NormalizeData(integrated_so, normalization.method = "LogNormalize", scale.factor = 10000)
    integrated_so <- FindVariableFeatures(integrated_so, selection.method = "vst", nfeatures = 2000)
    integrated_so <- ScaleData(integrated_so, features = rownames(integrated_so), verbose = FALSE)
    integrated_so <- CellCycleScoring(integrated_so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
    integrated_so <- RunPCA(integrated_so)
    integrated_so <- RunUMAP(integrated_so, dims = 1:30, reduction = "pca")
    print(DimPlot(integrated_so, group.by = c("sample_name", "project")))
    integrated_so <- IntegrateLayers(
        object = integrated_so,
        method = CCAIntegration,
        orig.reduction = "pca",
        new.reduction = "cca",
        # scale.layer = "scale.data",
        # features = rownames(integrated_so)
    )
    integrated_so <- RunUMAP(integrated_so, dims = 1:30, reduction = "cca")
    print(DimPlot(integrated_so, group.by = c("sample_name", "project")))

    integrated_so <- RunHarmony(integrated_so, group.by.vars = "project", theta = 10) # try with theta = 30 as well
    integrated_so <- RunUMAP(integrated_so, dims = 1:30, reduction = "harmony")
    print(DimPlot(integrated_so, group.by = c("sample_name", "project")))

    # qs::qsave(integrated_so, file.path(objects_folder, paste0("cao_", new_organ, "_organoid_integrated_filtered_normalized.qs")), nthreads = 30)
}

