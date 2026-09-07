library(Seurat)
library(harmony)
library(ggplot2)

project_folder <- ""
fetal_dev_project_folder <- ""
hnf_dev_project_folder <- ""
primary_dev_project_folder <- ""

# GB
## organoid
gb_organoid <- readRDS(file.path(project_folder, "R_objects", "aggregated", "HV", "by_organ", "aggr_filtered_coding_highrp-organ-GB.rds"))
DefaultAssay(gb_organoid) <- "RNA"
gb_organoid@assays$SCT <- NULL
gb_organoid@meta.data <- gb_organoid@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "organ", "age", "sample_name", "condition")]
gb_organoid$project <- "organoid"
gb_organoid$organ <- "GB"
gb_organoid$condition <- "N/A"
gb_organoid[["RNA"]] <- CreateAssay5Object(
    counts = GetAssayData(gb_organoid, assay = "RNA", slot = "counts"),
    data = GetAssayData(gb_organoid, assay = "RNA", slot = "data")
)
Key(gb_organoid[["RNA"]]) <- "rna_"
gb_organoid@reductions$pca <- NULL
gb_organoid@reductions$umap <- NULL

## Fetal dev
gb_fetal <- readRDS(file.path(fetal_dev_project_folder, "R_objects", "samples", "filter-100k-counts", "not_used", "so-gb-2274-sctransformed.rds"))
DefaultAssay(gb_fetal) <- "RNA"
gb_fetal@assays$SCT <- NULL
gb_fetal@assays$ATAC <- NULL
gb_fetal <- gb_fetal[intersect(rownames(gb_fetal), rownames(gb_organoid)), ]
gb_fetal@meta.data <- gb_fetal@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp")]
gb_fetal$organ <- "GB"
gb_fetal$age <- "11w"
gb_fetal$sample_name <- "2274 GB"
gb_fetal$project <- "primary Fetal Development"
gb_fetal[["RNA"]] <- CreateAssay5Object(
    counts = GetAssayData(gb_fetal, assay = "RNA", slot = "counts"),
    data = GetAssayData(gb_fetal, assay = "RNA", slot = "data")
)
Key(gb_fetal[["RNA"]]) <- "rna_"
gb_fetal@reductions$pca <- NULL
gb_fetal@reductions$umap <- NULL
gb_fetal@reductions$lsi <- NULL
gb_fetal@reductions$umap.atac <- NULL

## Hnf1b
gb_hnf <- readRDS(file.path(hnf_dev_project_folder, "R_objects", "aggregate", "so_sct_aggregate.rds"))
colnames(gb_hnf@meta.data)
gb_hnf <- subset(gb_hnf, sample_name == "2373 GB (nuclei)")
DefaultAssay(gb_hnf) <- "RNA"
gb_hnf@assays$SCT <- NULL
gb_hnf <- gb_hnf[intersect(rownames(gb_hnf), rownames(gb_organoid)), ]
gb_hnf@meta.data <- gb_hnf@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "age", "sample_name", "organ" )]
gb_hnf$project <- "primary HNF1b"
gb_hnf$condition <- "N/A"
gb_hnf[["RNA"]] <- CreateAssay5Object(
    counts = GetAssayData(gb_hnf, assay = "RNA", slot = "counts"),
    data = GetAssayData(gb_hnf, assay = "RNA", slot = "data")
)
Key(gb_hnf[["RNA"]]) <- "rna_"
gb_hnf@reductions$pca <- NULL
gb_hnf@reductions$umap <- NULL



integrated_gb_so <- merge(gb_organoid, y = list(gb_fetal, gb_hnf), add.cell.ids = c("organoid", "primary Fetal Development", "primary HNF1b"))
integrated_gb_so <- NormalizeData(integrated_gb_so, normalization.method = "LogNormalize", scale.factor = 10000)
integrated_gb_so <- FindVariableFeatures(integrated_gb_so, selection.method = "vst", nfeatures = 4000)
integrated_gb_so <- ScaleData(integrated_gb_so, features = rownames(integrated_gb_so), verbose = FALSE)
integrated_gb_so <- CellCycleScoring(integrated_gb_so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
integrated_gb_so <- RunPCA(integrated_gb_so)
integrated_gb_so <- RunUMAP(integrated_gb_so, dims = 1:30, reduction = "pca")
integrated_gb_so <- RunHarmony(object = integrated_gb_so, group.by.vars = "project", reduction.use = "pca", theta = 2)
integrated_gb_so <- RunUMAP(integrated_gb_so, dims = 1:30, reduction = "harmony")
integrated_gb_so <- IntegrateLayers(
    object = integrated_gb_so,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "cca"
)
integrated_gb_so <- RunUMAP(integrated_gb_so, dims = 1:30, reduction = "cca")
DimPlot(integrated_gb_so, group.by = c("project", "sample_name"))
qs::qsave(integrated_gb_so, file.path(project_folder, "R_objects", "primary_integrated_gb_so.qs"), nthreads = 30)


# LIVER
liver_fetal <- qs::qread(file.path(primary_dev_project_folder, "objects", "R", "seurat", "aggregate_strong_filtering_without_cc_without_14w_sct.qs"), nthreads = 30)
liver_fetal <- subset(liver_fetal, age %in% c("07w", "09w", "11w", "13w"))
DefaultAssay(liver_fetal) <- "RNA"
liver_fetal@assays$SCT <- NULL
liver_fetal@meta.data <- liver_fetal@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "age", "sample_name")]
liver_fetal$project <- "primary Fetal Development"
liver_fetal$organ <- "Liver"
liver_fetal$condition <- "N/A"
liver_fetal[["RNA"]] <- CreateAssay5Object(
    counts = GetAssayData(liver_fetal, assay = "RNA", slot = "counts"),
    data = GetAssayData(liver_fetal, assay = "RNA", slot = "data")
)
Key(liver_fetal[["RNA"]]) <- "rna_"
liver_fetal@reductions$pca <- NULL
liver_fetal@reductions$umap <- NULL

liver_organoid <- readRDS(file.path(project_folder, "R_objects", "aggregated", "HV", "by_organ", "aggr_filtered_coding_highrp-organ-Liver.rds"))
DefaultAssay(liver_organoid) <- "RNA"
liver_organoid@assays$SCT <- NULL
liver_organoid@meta.data <- liver_organoid@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp", "organ", "age", "sample_name", "condition")]
liver_organoid$project <- "organoid"
liver_organoid[["RNA"]] <- CreateAssay5Object(
    counts = GetAssayData(liver_organoid, assay = "RNA", slot = "counts"),
    data = GetAssayData(liver_organoid, assay = "RNA", slot = "data")
)
Key(liver_organoid[["RNA"]]) <- "rna_"
liver_organoid@reductions$pca <- NULL
liver_organoid@reductions$umap <- NULL


integrated_liver_so <- merge(liver_fetal, liver_organoid, add.cell.ids = c("primary Fetal Development", "organoid"))
# integrated_liver_so <- JoinLayers(integrated_liver_so)
integrated_liver_so <- NormalizeData(integrated_liver_so, normalization.method = "LogNormalize", scale.factor = 10000)
integrated_liver_so <- FindVariableFeatures(integrated_liver_so, selection.method = "vst", nfeatures = 4000)
integrated_liver_so <- ScaleData(integrated_liver_so, features = rownames(integrated_liver_so), verbose = FALSE)

# integrated_liver_so <- SCTransform(integrated_liver_so, assay = "RNA", verbose = FALSE)
integrated_liver_so <- CellCycleScoring(integrated_liver_so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
integrated_liver_so <- RunPCA(integrated_liver_so)
integrated_liver_so <- RunUMAP(integrated_liver_so, dims = 1:30, reduction = "pca")

integrated_liver_so <- RunHarmony(object = integrated_liver_so, group.by.vars = "sample_name", reduction.use = "pca", theta = 2)
integrated_liver_so <- RunUMAP(integrated_liver_so, dims = 1:30, reduction = "harmony")



integrated_liver_so <- IntegrateLayers(
    object = integrated_liver_so,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "cca"
)
integrated_liver_so <- RunUMAP(integrated_liver_so, dims = 1:30, reduction = "cca")

DimPlot(integrated_liver_so, group.by = c("project", "sample_name"))
qs::qsave(integrated_liver_so, file.path(project_folder, "R_objects", "primary_integrated_liver_so.qs"), nthreads = 30)

