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

create.so = function(h5_path, use_pc_genes = TRUE){
    counts <- Read10X_h5(h5_path)
    
    so <-  CreateSeuratObject(
        counts = counts,
        assay = "RNA"
    )
    
    return(so)
}

# load up date
current.files = file.path(counts_files, list.files(counts_files))

so_list <- lapply(current.files, function(x) {create.so(x, use_pc_genes = TRUE)})
short_name <- sapply(list.files(counts_files), function(x) { strsplit(x, "[.]")[[1]][1]})
names(so_list) <- short_name

so_merged = merge(
  x = so_list[[1]],
  y = so_list[2:3],
  add.cell.ids = short_name)

mt.genes=grep("^MT-", rownames(so_merged), value=FALSE)
rp.genes=grep("^RP[SL]", rownames(so_merged), value=FALSE)
mt.rp.genes=grep("^MRP[SL]", rownames(so_merged), value=FALSE)
so_merged[['percent.mt']] = PercentageFeatureSet(so_merged, features=rownames(so_merged)[mt.genes])
so_merged[['percent.rp']] = PercentageFeatureSet(so_merged, features=rownames(so_merged)[c(rp.genes, mt.rp.genes)])
so_merged = subset(so_merged, features=-c(mt.rp.genes,mt.genes, rp.genes))
so_merged$percent.mt[is.na(so_merged$percent.mt)] <- 100
so_merged$percent.rp[is.na(so_merged$percent.rp)] <- 100
so_merged$barcode <- sapply(colnames(so_merged), function(x) { strsplit(x, "_")[[1]][1]})


so_merged2 <- subset(so_merged, percent.rp < 5 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA > 1000 & nCount_RNA < 1e+04)

pdf("preprocessing/qc/qc.pdf", width = 15, height = 6)
print(VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol=4, pt.size=0, log = TRUE, raster = TRUE))
print(VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol=4, pt.size=0, log = FALSE, raster = TRUE))
print(VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol=4, pt.size=0, log = TRUE, group.by = "barcode", raster = TRUE))
print(VlnPlot(so_merged2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol=4, pt.size=0, log = FALSE, group.by = "barcode", raster = TRUE))
print(rasterize(ggplot(so_merged2@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(alpha=0.03)+ ggtitle("after filtering - aggregate")))
print(rasterize(ggplot(so_merged2@meta.data, aes(x=percent.mt, y=percent.rp)) + geom_point(alpha=0.03)+ ggtitle("after filtering - aggregate")))
print(rasterize(ggplot(so_merged2@meta.data, aes(x=percent.mt, y=nCount_RNA)) + geom_point(alpha=0.03)+ ggtitle("after filtering - aggregate")))
print(rasterize(ggplot(so_merged2@meta.data, aes(x=percent.rp, y=nCount_RNA)) + geom_point(alpha=0.03)+ ggtitle("after filtering - aggregate")))

for (mtd_val in unique(so_merged2$barcode)) {
  print(rasterize(ggplot(so_merged2[, so_merged2$barcode == mtd_val]@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(alpha=0.03)+ ggtitle(paste0("after filtering ", mtd_val))))
  print(rasterize(ggplot(so_merged2[, so_merged2$barcode == mtd_val]@meta.data, aes(x=percent.mt, y=percent.rp)) + geom_point(alpha=0.03)+   ggtitle(paste0("after filtering ", mtd_val))))
  print(rasterize(ggplot(so_merged2[, so_merged2$barcode == mtd_val]@meta.data, aes(x=percent.mt, y=nCount_RNA)) + geom_point(alpha=0.03)+   ggtitle(paste0("after filtering ", mtd_val))))
  print(rasterize(ggplot(so_merged2[, so_merged2$barcode == mtd_val]@meta.data, aes(x=percent.rp, y=nCount_RNA)) + geom_point(alpha=0.03)+   ggtitle(paste0("after filtering ", mtd_val))))
}
dev.off()


so_merged2 <- SCTransform(so_merged2, variable.features.n = 3000, return.only.var.genes = FALSE)
so_merged2 <- RunPCA(so_merged2, features = so_merged2@assays$SCT@var.features[1:2000])
so_merged2 <- RunUMAP(so_merged2, reduction = "pca", dims = 1:30)

saveRDS(so_merged2, "R_objects/aggr_filtered.rds")

DimPlot(so_merged2, group.by = "barcode")

scConf <- ShinyCell::createConfig(so_merged2)
ShinyCell::makeShinyApp(so_merged2, scConf,
                      gene.mapping = TRUE,
                      gex.assay='SCT', 
                      shiny.dir = paste0(project_folder, "/output/shiny_apps/", "dummy_shiny_app"), 
                      shiny.title=paste0('Floris snRNAseq - initial analysis'))


## Remove mesenchymal

so_merged2 <- readRDS("R_objects/aggr_filtered.rds")
DimPlot(so_merged2)

mask_cells <- (so_merged2@reductions$umap@cell.embeddings[, 1] > 0 & so_merged2@reductions$umap@cell.embeddings[, 2] < -4)
so_merged2 <- so_merged2[ , !mask_cells]

so_merged2 <- SCTransform(so_merged2, variable.features.n = 3000, return.only.var.genes = FALSE)
so_merged2 <- RunPCA(so_merged2, features = so_merged2@assays$SCT@var.features[1:2000])
so_merged2 <- RunUMAP(so_merged2, reduction = "pca", dims = 1:30)

saveRDS(so_merged2, "R_objects/aggr_filtered_without_mesenchymal.rds")
DimPlot(so_merged2)

## Only ALI
so_merged <- readRDS("R_objects/aggr_filtered_without_mesenchymal.rds")
so_merged <- subset(so_merged, barcode == "SITTF7")
so_merged <- SCTransform(so_merged, variable.features.n = 3000, return.only.var.genes = FALSE)
so_merged <- RunPCA(so_merged, features = so_merged@assays$SCT@var.features[1:2000])
so_merged <- RunUMAP(so_merged, reduction = "pca", dims = 1:30)

saveRDS(so_merged, "R_objects/SITTF7_ALI.rds")
