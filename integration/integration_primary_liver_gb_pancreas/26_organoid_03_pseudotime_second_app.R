source("scripts/pseudotime_subset_generate_shiny.R") # check Core-Bioinformatics/Starlng for an up-to-date version

library(Seurat)
library(monocle3)
so <- qs::qread("R_objects/organoid_liver_hepatoblasts_filtered_normalized.qs")
ca <- qs::qread("R_objects/clustassess/organoid_liver_hepatoblasts.qs")


first_dir <- "pseudotime_shiny_apps/first_part_analysis/organoid_liver_hepatoblasts"
second_dir <- "pseudotime_shiny_apps/second_part_analysis/organoid_liver_hepatoblasts"

if (!dir.exists(second_dir)) {
    dir.create(second_dir, recursive = TRUE)
}

write_object(
    mon_obj = readRDS(file.path(first_dir, "monocle_object.rds")),
    trajectory_id = "organoid_hepatoblasts_subset",
    start_genes = "AFP",
    start_expression_thresh = 6.18,
    start_relax_ngenes = 0,
    end_genes = c("CES1", "CYP2E1"),
    end_expression_thresh = 0.87,
    end_relax_ngenes = 0,
    ncores = 30,
    use_closed_loops = FALSE,
    use_partitions = FALSE,
    nodes_per_log10_cells = 75,
    learn_graph_controls = list(
        eps = 1e-5,
        maxiter = 100
    ),
    output_dir = second_dir
)
