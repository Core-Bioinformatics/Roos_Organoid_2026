library(Seurat)
library(ClustAssess)

combinations <- list(
    library = c("6w_epithelial", "8w_epithelial")
)

for (lbr in combinations$library) {
    print(lbr)
    so <- readRDS(glue::glue("{project_folder}/data/R_objects/{lbr}.rds"))
    assay_name <- "SCT"
    clustassess_object <- readRDS(glue::glue("{project_folder}/data/R_objects/clustassess/{lbr}_clustassess.rds"))

    write_shiny_app(
        seurat_object = so,
        assay_name = "SCT",
        clustassess_object = clustassess_object,
        project_folder = glue::glue("{project_folder}/shiny/clustassess/{lbr}")
    )
}
