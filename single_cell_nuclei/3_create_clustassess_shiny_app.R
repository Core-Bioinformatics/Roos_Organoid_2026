library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)

prefix <- "aggr_filtered" # or ""
high_res_sufix <- "_res_3"
so <- readRDS(file.path(project_folder, "R_objects", paste0(prefix, ".rds"))) # _without_mesenchymal

clustassess_object <- readRDS(file.path(project_folder, "R_objects", "clustassess", glue::glue("{prefix}{high_res_sufix}.rds")))

sct_matrix <- as.matrix(so@assays$SCT@data)
sct_matrix <- rbind(sct_matrix, so@assays$RNA@data["FOXI1", ])
rownames(sct_matrix)[nrow(sct_matrix)] <- "FOXI1"

so@assays$SCT@data <- as(sct_matrix, "sparseMatrix")


write_shiny_app(
  object = so,
  assay_name = "SCT",
  clustassess_object = clustassess_object,
  project_folder = glue::glue("{project_folder}/output/clustassess_shiny/{prefix}{high_res_sufix}"),
  shiny_app_title = prefix
)
