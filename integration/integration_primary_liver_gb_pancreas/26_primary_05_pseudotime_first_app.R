
source("scripts/pseudotime_roots_generate_shiny.R")

library(Seurat)
library(monocle3)
so <- qs::qread("R_objects/primary_liver_fetal_hepatoblasts_so.qs")
so <- JoinLayers(so)
ca <- qs::qread("R_objects/clustassess/primary_liver_fetal_hepatoblasts.qs")

x <- write_object.Seurat(
    seurat_object = so,
    clustassess_object = ca,
    assay_used = "RNA",
    app_name = "primary Liver Fetal Hepatoblasts",
    stable_config = list(
        ftype = "Highly_Variable",
        fsize = 4000,
        clmethod = "SLM",
        k = 6
    ),
    output_dir = "pseudotime_shiny_apps/first_part_analysis/primary_liver_fetal_hepatoblasts",
    use_closed_loops = FALSE,
    learn_graph_controls = list(
        eps = 1e-5,
        maxiter = 100
    ),
    nodes_per_log10_cells = 75
)


