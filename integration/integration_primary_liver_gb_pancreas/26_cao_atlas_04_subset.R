library(Seurat)

project_folder <- ""
objects_folder <- file.path(project_folder, "R_objects")

target_subsets <- list(
    "Intestine" = "",
    "Pancreas" = c("Ductal cells"),
    "Liver" = c("Hepatoblasts")
)

for (organ in names(target_subsets)) {
    target_subset <- target_subsets[[organ]]
    print(organ)
    if (target_subset == "" || is.null(target_subset)) {
        next
    }
    so <- qs::qread(file.path(objects_folder, paste0("cao_", organ, "_filtered_normalized.qs")), nthreads = 30)
    so <- subset(so, subset = Main_cluster_name %in% target_subset)

    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
    so <- ScaleData(so, features = rownames(so), verbose = FALSE)
    so <- RunPCA(so)
    so <- RunUMAP(so, dims = 1:30, reduction = "pca")

    qs::qsave(so, file.path(objects_folder, paste0("cao_", organ, "_", paste(target_subset, collapse = "_"), "_filtered_normalized_subset.qs")), nthreads = 30)
}
 e