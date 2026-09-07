library(Seurat)
library(harmony)
library(ggplot2)

project_folder <- ""
objects_folder <- file.path(project_folder, "R_objects")


# Liver
so_liver <- qs::qread(file.path(objects_folder, "primary_integrated_liver_so.qs"), nthreads = 30)
ca_liver <- qs::qread(file.path(objects_folder, "clustassess", "primary_integrated_liver_cca.qs"), nthreads = 30)
used_partition <- ClustAssess::get_clusters_from_clustassess_object(
    clustassess_object = ca_liver,
    feature_type = "Highly_Variable",
    feature_size = 3000,
    clustering_method = "SLM",
    nclusters = 18
)[[1]]$partitions[[1]]$mb
target_clusters <- c(3, 8, 9, 5, 12)
mask_target <- used_partition %in% target_clusters
# mask_project <- so_liver$project == "primary Fetal Development"
removed_cells <- colnames(so_liver)[mask_target]
so_liver <- subset(so_liver, cells = removed_cells)
so_liver <- NormalizeData(so_liver, normalization.method = "LogNormalize", scale.factor = 10000)
so_liver <- FindVariableFeatures(so_liver, selection.method = "vst", nfeatures = 4000)
so_liver <- ScaleData(so_liver, features = rownames(so_liver), verbose = FALSE)

qs::qsave(so_liver, file.path(objects_folder, "primary_integrated_liver_filtered_so.qs"), nthreads = 30)

# GB
so_gb <- qs::qread(file.path(objects_folder, "primary_integrated_gb_so.qs"), nthreads = 30)
ca_gb <- qs::qread(file.path(objects_folder, "clustassess", "primary_integrated_gb_harmony_10.qs"), nthreads = 30)

norm_matrix <- GetAssayData(JoinLayers(so_gb), assay = "RNA", layer = "data")
used_partition <- ClustAssess::get_clusters_from_clustassess_object(
    clustassess_object = ca_gb,
    feature_type = "Highly_Variable",
    feature_size = 2000,
    clustering_method = "SLM",
    nclusters = 16
)[[1]]$partitions[[1]]$mb
target_clusters <- c(1, 11)
mask_target <- used_partition %in% target_clusters
mask_project <- so_gb$project != "organoid"
mask_expression <- norm_matrix["KRT19", ] > 0
kept_cells <- colnames(so_gb)[mask_target & mask_project & mask_expression]
removed_cells <- setdiff(colnames(so_gb)[so_gb$project != "organoid"], kept_cells)
length(removed_cells)
so_gb <- subset(so_gb, cells = setdiff(colnames(so_gb), removed_cells))
so_gb <- NormalizeData(so_gb, normalization.method = "LogNormalize", scale.factor = 10000)
so_gb <- FindVariableFeatures(so_gb, selection.method = "vst", nfeatures = 4000)
so_gb <- ScaleData(so_gb, features = rownames(so_gb), verbose = FALSE)

qs::qsave(so_gb, file.path(objects_folder, "primary_integrated_gb_filtered_so.qs"), nthreads = 30)
