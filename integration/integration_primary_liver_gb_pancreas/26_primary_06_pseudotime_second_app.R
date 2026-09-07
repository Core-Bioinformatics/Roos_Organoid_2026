source("scripts/pseudotime_subset_generate_shiny.R")

library(Seurat)
library(monocle3)
so <- qs::qread("R_objects/primary_liver_fetal_hepatoblasts_so.qs")
so <- JoinLayers(so)
ca <- qs::qread("R_objects/clustassess/primary_liver_fetal_hepatoblasts.qs")


first_dir <- "pseudotime_shiny_apps/first_part_analysis/primary_liver_fetal_hepatoblasts"
second_dir <- "pseudotime_shiny_apps/second_part_analysis/primary_liver_fetal_hepatoblasts"

if (!dir.exists(second_dir)) {
    dir.create(second_dir, recursive = TRUE)
}

write_object(
    mon_obj = readRDS(file.path(first_dir, "monocle_object.rds")),
    trajectory_id = "primary_liver_fetal_hepatoblasts_subset",
    start_genes = c("GSTP1", "FST", "DLK1"),
    start_expression_thresh = 2.7,
    start_relax_ngenes = 1,
    end_genes = c("LIPC", "CYP3A7", "ANGPTL3", "LEPR", "GC"),
    end_expression_thresh = 3.23,
    end_relax_ngenes = 2,
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

