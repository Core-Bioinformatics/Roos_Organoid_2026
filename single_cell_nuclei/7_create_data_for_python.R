library(Seurat)
library(readr)
library(dplyr)
library(tidyr)
library(ggrastr)
library(ggplot2)
library(ggplotify)
source("pseudotime_roots_generate_shiny.R")

assay_name <- "SCT"

psd_configs <- read.csv(file.path(project_folder, "metadata", "sample_combinations.csv"), comment.char = "#")
print(psd_configs)
project_folder <- file.path(project_folder, "output/data_for_python")
combinations <- c()
seurat_ca_metadata_path <- "metadata/seurat_paths_to_clustassess.csv"
seurat_ca_metadata <- read.csv(seurat_ca_metadata_path, comment.char = "#")

common_seurat_dir <- find_common_directory(seurat_ca_metadata$seurat)

se_prefix <- "aggr_filtered_coding_highrp-"

for (i in seq_len(nrow(psd_configs))) {
    sampleid <- psd_configs$Id[i]
    ftype <- as.character(psd_configs$Ftype[i])
    fsize <- as.character(psd_configs$Fsize[i])
    clm <- psd_configs$Cl_method[i]
    k <- as.character(psd_configs$k[i])

    if (length(combinations) > 0 && sampleid %in% combinations) {
        next
    }

    corresp_index <- find_corresp_index(sampleid, seurat_ca_metadata, se_prefix)
    print(corresp_index)

    if (corresp_index == -1) {
        next
    }
    
    seurat_path <- dirname(seurat_ca_metadata$seurat[corresp_index])
    clustassess_path <- seurat_ca_metadata$clustassess[corresp_index]

    to_add_folders <- c(sampleid, ftype, fsize)

    print(to_add_folders)

    while (basename(seurat_path) != common_seurat_dir) {
        to_add_folders <- c(basename(seurat_path), to_add_folders)
        seurat_path <- dirname(seurat_path)
    }
    seurat_path <- seurat_ca_metadata$seurat[corresp_index]

    to_add_folders <- as.list(to_add_folders)
    target_folder <- project_folder

    for (j in seq_along(to_add_folders)) {
        target_folder <- file.path(target_folder, to_add_folders[[j]])
    }

    if (file.exists(file.path(target_folder, "metadata.csv"))) {
        next
    }

    if (file.exists(file.path(target_folder, "genes.txt"))) {
        next
    }

    print(paste(target_folder, sampleid, "\n"))
    dir.create(target_folder, showWarnings = FALSE, recursive = TRUE)



    ca_object <- readRDS(clustassess_path)
    se_object <- readRDS(seurat_path)
    
    meta_data <- se_object@meta.data
    meta_data[["library"]] <- meta_data[["barcode"]]
    meta_data[["UMAP1"]] <- ca_object[[ftype]][[fsize]]$umap[, 1]
    meta_data[["UMAP2"]] <- ca_object[[ftype]][[fsize]]$umap[, 2]

    npcas <- ncol(ca_object[[ftype]][[fsize]]$pca)
    for (j in seq_len(npcas)) {
        meta_data[[paste0("PCA", j)]] <- ca_object[[ftype]][[fsize]]$pca[, j]
    }

    meta_data[["stable_clusters"]] <- factor(as.integer(ca_object[[ftype]][[fsize]]$clustering_stability$split_by_k[[clm]][[k]]$partitions[[1]]$mb))
    meta_data[["barcode"]] <- rownames(meta_data)

    write.csv(meta_data, file.path(target_folder, "metadata.csv"))
    write(
        paste(ca_object[[ftype]]$feature_list[seq_len(as.integer(fsize))], collapse = "\n"),
        file.path(target_folder, "genes.txt")
    )
}


