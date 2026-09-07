
source("scripts/pseudotime_roots_generate_shiny.R") # check Core-Bioinformatics/Starlng for an up-to-date version

library(Seurat)
library(monocle3)
so <- qs::qread("R_objects/organoid_liver_hepatoblasts_filtered_normalized.qs")
ca <- qs::qread("R_objects/clustassess/organoid_liver_hepatoblasts.qs")

x <- write_object.Seurat(
    seurat_object = so,
    clustassess_object = ca,
    assay_used = "RNA",
    app_name = "Organoid Liver Hepatoblasts",
    stable_config = list(
        ftype = "Highly_Variable",
        fsize = 3500,
        clmethod = "SLM",
        k = 6
    ),
    output_dir = "output/pseudotime_shiny_apps/first_part_analysis/organoid_liver_hepatoblasts",
    use_closed_loops = FALSE,
    learn_graph_controls = list(
        eps = 1e-5,
        maxiter = 100
    ),
    nodes_per_log10_cells = 75
)

