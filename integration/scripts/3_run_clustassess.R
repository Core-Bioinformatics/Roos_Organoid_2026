library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)

combinations <- list(
     library = c("6w_epithelial", "8w_epithelial")
)

RhpcBLASctl::blas_set_num_threads(1)
ncores <- 50
my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
)

doParallel::registerDoParallel(cl = my_cluster)
for (lbr in combinations$library) {
    print(lbr)
    so <- readRDS(glue::glue("{project_folder}/data/R_objects/{lbr}.rds"))
    assay_name <- "SCT"

    features <- dimnames(so@assays[[assay_name]])[[1]]
    var_features <- so@assays[[assay_name]]@var.features
    n_abundant <- 4000
    most_abundant_genes <- rownames(so@assays[[assay_name]])[order(Matrix::rowSums(so@assays[[assay_name]]),
        decreasing = TRUE
    )][1:n_abundant]
    gene_list <- list(
        "Highly_Variable" = so@assays[[assay_name]]@var.features,
        "Most_Abundant" = most_abundant_genes
    )

    steps_list <- list(
        "Highly_Variable" = seq(from = 500, by = 500, to = 4000),
        "Most_Abundant" = seq(from = 500, by = 500, to = 4000)
    )

    expr_matrix <- so@assays$SCT@scale.data
    rm(so)
    gc()

    test_automm <- automatic_stability_assessment(
        expression_matrix = expr_matrix,
        n_repetitions = 50,
        temp_file = "clustassess_temp.rds",
        n_neigh_sequence = seq(from = 5, to = 50, by = 5),
        resolution_sequence = seq(from = 0.1, to = 2, by = 0.1),
        features_sets = gene_list,
        steps = steps_list,
        n_top_configs = 2,
        npcs = 30,
        umap_arguments = list(
            min_dist = 0.3,
            n_neighbors = 30,
            metric = "cosine",
            init = "spca"
        )
    )

    saveRDS(test_automm, glue::glue("{project_folder}/data/R_objects/clustassess/{lbr}_clustassess.rds"))
}
parallel::stopCluster(cl = my_cluster)
foreach::registerDoSEQ()
