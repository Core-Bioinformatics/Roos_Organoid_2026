library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)

ncores <- 8
my_cluster <- parallel::makeCluster(
  ncores,
  type = "PSOCK"
)

read.table()

RhpcBLASctl::blas_set_num_threads(1)
doParallel::registerDoParallel(cl = my_cluster)
so <- readRDS(file.path(project_folder, "R_objects", "SITTF7_ALI.rds"))

assay_name <- "SCT"

expr_matrix <- so@assays[[assay_name]]@scale.data
features <- dimnames(so@assays[[assay_name]])[[1]]
var_features <- so@assays[[assay_name]]@var.features
n_abundant <- 3000
most_abundant_genes <- rownames(expr_matrix)[order(Matrix::rowSums(expr_matrix), decreasing = TRUE)][1:n_abundant]

gene_list = list(
  "Most_Abundant" = most_abundant_genes,
  "Highly_Variable" = so@assays[[assay_name]]@var.features
)

steps_list = list(
  "Most_Abundant" = seq(from = 500, by = 500, to = 3000),
  "Highly_Variable" = seq(from = 500, by = 500, to = 3000)
)

rm(so)
gc()

test_automm = automatic_stability_assessment(expression_matrix = expr_matrix,
                                             n_cores = ncores,
                                             n_repetitions = 50,
                                             temp_file = "clustassess_temp.rds",
                                             n_neigh_sequence = seq(from = 5, to = 50, by = 5),
                                             resolution_sequence = seq(from = 0.1, to = 2, by = 0.1), #c(0.5,1, 1.5, 2, 2.5, 3),
                                             features_sets = gene_list,
                                             steps = steps_list,
                                             n_top_configs = 2,
                                             npcs = 30,
                                             min_dist = 0.3,
                                             n_neighbors = 30,
                                             metric = "cosine")


saveRDS(test_automm, file.path(project_folder, "R_objects", "SITTF7_ALI_clustassess_object.rds"))
on.exit(parallel::stopCluster(cl = my_cluster))
