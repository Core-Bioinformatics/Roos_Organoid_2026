library(Seurat)

project_folder <- ""
objects_folder <- file.path(project_folder, "R_objects")

# Liver
so_liver <- qs::qread(file.path(objects_folder, "primary_integrated_liver_filtered_so.qs"), nthreads = 30)
so_liver <- subset(so_liver, cells = colnames(so_liver)[so_liver$project == "primary Fetal Development"])

so_liver <- NormalizeData(so_liver, normalization.method = "LogNormalize", scale.factor = 10000)
so_liver <- FindVariableFeatures(so_liver, selection.method = "vst", nfeatures = 4000)
so_liver <- ScaleData(so_liver, features = rownames(so_liver), verbose = FALSE)
so_liver <- CellCycleScoring(so_liver, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
so_liver <- RunPCA(so_liver)
so_liver <- RunUMAP(so_liver, dims = 1:30, reduction = "pca")

qs::qsave(so_liver, file.path(objects_folder, "primary_liver_fetal_hepatoblasts_so.qs"), nthreads = 30)

