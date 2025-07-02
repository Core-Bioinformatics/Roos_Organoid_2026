library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(gridExtra)
library(patchwork)
library(ggrastr)
theme_set(theme_bw())
set.seed(777)

today_date <- function() {
    format(Sys.Date(), "%y-%m-%d")
}

create.so <- function(h5_path, use_pc_genes = TRUE) {
    counts <- Read10X_h5(h5_path)

    pc_genes <- read.table("~/protein_coding_gene_names.txt")[[1]]
    npc_genes <- read.table("~/non_protein_coding_gene_names.txt")[[1]]
    pc_genes <- stringr::str_replace_all(string = pc_genes, pattern = "_", replacement = "-")
    npc_genes <- stringr::str_replace_all(string = npc_genes, pattern = "_", replacement = "-")

    first_filter_genes <- intersect(rownames(counts), ifelse(use_pc_genes, pc_genes, npc_genes))
    other_genes <- setdiff(setdiff(rownames(counts), ifelse(use_pc_genes, pc_genes, npc_genes)), ifelse(use_pc_genes, npc_genes, pc_genes))
    basic_name <- substr(other_genes, 1, regexpr(".", other_genes, fixed = T) - 1)
    second_filter_genes <- other_genes[which(basic_name %in% ifelse(use_pc_genes, pc_genes, npc_genes))]
    protein_coding_genes <- union(first_filter_genes, second_filter_genes)

    so <- CreateSeuratObject(
        counts = counts[intersect(rownames(counts), pc_genes), ],
        assay = "RNA"
    )
    # mt.genes=grep("^MT-", rownames(so), value=FALSE)
    # rp.genes=grep("^RP[SL]", rownames(so), value=FALSE)
    # mt.rp.genes=grep("^MRP[SL]", rownames(so), value=FALSE)
    # so[['percent.mt']] = PercentageFeatureSet(so, features=rownames(so)[mt.genes])
    # so[['percent.rp']] = PercentageFeatureSet(so, features=rownames(so)[c(mt.rp.genes, rp.genes)])
    # so = subset(so, features=-c(mt.genes, rp.genes, mt.rp.genes))

    return(so)
}

min.cells <- 0
min.features <- 0

# load up date
current.files <- paste0(floris_counts_files, list.files(floris_counts_files))

so_list <- lapply(current.files, function(x) {
    create.so(x, use_pc_genes = TRUE)
})
short_name <- sapply(list.files(floris_counts_files), function(x) {
    strsplit(x, "[.]")[[1]][1]
})
names(so_list) <- short_name

for (name in short_name) {
    if (file.exists(glue::glue("{project_folder}/data/R_objects/raw_{name}.rds"))) {
        so_list[[name]] <- readRDS(glue::glue("{project_folder}/data/R_objects/raw_{name}.rds"))
    } else {
        saveRDS(so_list[[name]], glue::glue("{project_folder}/data/R_objects/raw_{name}.rds"))
    }
}

print("Merge the objects")
so_merged <- merge(
    x = so_list[[1]],
    y = so_list[2:length(so_list)],
    add.cell.ids = short_name
)

rm(so_list)
gc()

sample_names <- sapply(so_merged@assays[["RNA"]]@data@Dimnames[[2]], function(x) {
    substr(x, 0, nchar(x) - 19)
})
so_merged@meta.data[["library"]] <- factor(sample_names)

samples_metadata <- read.csv(file.path(project_folder, "metadata", "metadata.csv"))

age <- rep("", ncol(so_merged))
for (i in 1:nrow(samples_metadata)) {
    age[so_merged@meta.data$library == samples_metadata$Barcode[i]] <- samples_metadata$After_fertilization[i]
}
so_merged@meta.data[["age"]] <- factor(age)

saveRDS(so_merged, glue::glue("{project_folder}/data/R_objects/raw_merged_coding.rds"))

so_merged2 <- subset(so_merged, nFeature_RNA > 1000 & nCount_RNA > 100 & nCount_RNA < 2e4)

mt.genes <- grep("^MT-", rownames(so_merged2), value = FALSE)
rp.genes <- grep("^RP[SL]", rownames(so_merged2), value = FALSE)
mt.rp.genes <- grep("^MRP[SL]", rownames(so_merged2), value = FALSE)
so_merged2[["percent.mt"]] <- PercentageFeatureSet(so_merged2, features = rownames(so_merged2)[mt.genes])
so_merged2[["percent.rp"]] <- PercentageFeatureSet(so_merged2, features = rownames(so_merged2)[c(rp.genes, mt.rp.genes)])
so_merged2 <- subset(so_merged2, features = -c(mt.rp.genes, mt.genes, rp.genes))
so_merged2 <- subset(so_merged2, percent.mt < 15 & percent.rp < 40)

for (libr in unique(so_merged2$library)) {
    subso <- subset(so_merged2, library == libr)
    subso <- SCTransform(subso, return.only.var.genes = F, verbose = F, variable.features.n = 5000)

    saveRDS(subso, file.path(project_folder, "data", "R_objects", glue::glue("filtered_{libr}.rds")))
}



pdf(glue::glue("{project_folder}/preprocessing/qc/qc.pdf"), width = 15, height = 6)

p <- VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0, log = TRUE) + plot_annotation(title = "After filtering - log scale")
print(rasterize(p))
p <- VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0) + plot_annotation(title = "After filtering")
print(rasterize(p))
p <- VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0, group.by = "library", log = TRUE) + plot_annotation(title = "After filtering - log scale")
print(rasterize(p))
p <- VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0, group.by = "library") + plot_annotation(title = "After filtering")
print(rasterize(p))
p <- VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0, group.by = "age", log = TRUE) + plot_annotation(title = "After filtering - log scale")
print(rasterize(p))
p <- VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4, pt.size = 0, group.by = "age") + plot_annotation(title = "After filtering")
print(rasterize(p))

p <- ggplot(so_merged2@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(alpha = 0.03) +
    ggtitle("after filtering - aggregate")
print(rasterize(p))
p <- ggplot(so_merged2@meta.data, aes(x = percent.mt, y = percent.rp)) +
    geom_point(alpha = 0.03) +
    ggtitle("after filtering - aggregate")
print(rasterize(p))
p <- ggplot(so_merged2@meta.data, aes(x = percent.mt, y = nCount_RNA)) +
    geom_point(alpha = 0.03) +
    ggtitle("after filtering - aggregate")
print(rasterize(p))
p <- ggplot(so_merged2@meta.data, aes(x = percent.rp, y = nCount_RNA)) +
    geom_point(alpha = 0.03) +
    ggtitle("after filtering - aggregate")
print(rasterize(p))

#
for (metadata_names in c("library", "age")) {
    metadata_values <- levels(so_merged2@meta.data[[metadata_names]])
    print(metadata_values)

    for (v in metadata_values) {
        p <- ggplot(so_merged2[, so_merged2[[metadata_names]] == v]@meta.data, aes(x = nCount_RNA, y = nFeature_RNA)) +
            geom_point(alpha = 0.03) +
            ggtitle(paste("after filtering", metadata_names, ":", v))
        print(rasterize(p))
        p <- ggplot(so_merged2[, so_merged2[[metadata_names]] == v]@meta.data, aes(x = percent.mt, y = percent.rp)) +
            geom_point(alpha = 0.03) +
            ggtitle(paste("after filtering", metadata_names, ":", v))
        print(rasterize(p))
        p <- ggplot(so_merged2[, so_merged2[[metadata_names]] == v]@meta.data, aes(x = percent.mt, y = nCount_RNA)) +
            geom_point(alpha = 0.03) +
            ggtitle(paste("after filtering", metadata_names, ":", v))
        print(rasterize(p))
        p <- ggplot(so_merged2[, so_merged2[[metadata_names]] == v]@meta.data, aes(x = percent.rp, y = nCount_RNA)) +
            geom_point(alpha = 0.03) +
            ggtitle(paste("after filtering", metadata_names, ":", v))
        print(rasterize(p))
    }
}



dev.off()

so_merged2 <- SCTransform(so_merged2, return.only.var.genes = F, verbose = F, variable.features.n = 5000)
so_merged2 <- RunPCA(so_merged2, npcs = 50, verbose = F)
so_merged2 <- RunUMAP(so_merged2, dims = 1:50, reduction = "pca", verbose = F)
saveRDS(so_merged2, glue::glue("{project_folder}/data/R_objects/merged_filtered_sct.rds"))
