
# this is a precursor of the Starlng app
# we recommend using the Starlng framework instead https://github.com/Core-Bioinformatics/Starlng
source("pseudotime_subset_generate_shiny.R")


write_object.data_frame(
    ids_df = read.csv("metadata/pseudotime_genes.csv", comment.char = "#"),
    seurat_ca_corresp = read.csv("metadata/seurat_paths_to_clustassess.csv", comment.char = "#"),
    target_app_dir_first_part = "output/pseudotime_apps/first_part_analysis",
    target_app_dir_second_part = "output/pseudotime_apps/second_part_analysis",
    use_closed_loops = FALSE,
    prefix = "",
    ncores = 75,
    learn_graph_controls = list(
        eps = 1e-5,
        maxiter = 100
    ),
    nodes_per_log10_cells = 30 
)