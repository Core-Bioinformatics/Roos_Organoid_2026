library(Seurat)

sample_names <- c(
    "WSSS8012016",
    "5478STDY7698210",
    "5891STDY9030809"
)

sample_so <- lapply(sample_names, function(x) {
    readRDS(file.path(project_folder,
        "data/R_objects",
        glue::glue("filtered_{x}.rds")
    ))
})
names(sample_so) <- sample_names

clustassess_so <- lapply(sample_names, function(x) {
    readRDS(file.path(project_folder,
        "data/R_objects/clustassess",
        glue::glue("{x}.rds")
    ))
})
names(clustassess_so) <- sample_names

# keep only epithelial cells - 6w
sample_so$"WSSS8012016"$stable_clusters <- factor(as.numeric(clustassess_so$WSSS8012016$"Highly_Variable"$"1500"$clustering_stability$split_by_k$"SLM"$"9"$partitions[[1]]$mb))
sample_so$"5478STDY7698210"$stable_clusters <- factor(as.numeric(clustassess_so$"5478STDY7698210"$"Highly_Variable"$"2000"$clustering_stability$split_by_k$"SLM"$"4"$partitions[[1]]$mb))
# sample_so$"5891STDY9030809"$stable_clusters <- factor(as.numeric(clustassess_so$"5891STDY9030809"$"Most_Abundant"$"3000"$clustering_stability$split_by_k$"SLM"$"9"$partitions[[1]]$mb))

so <- merge(
    subset(sample_so$"WSSS8012016", stable_clusters == 3),
    subset(sample_so$"5478STDY7698210", stable_clusters == 2)
)
so <- SCTransform(so, variable.features.n = 5000, return.only.var.genes = F, verbose = F)

so$stable_clusters <- NULL
saveRDS(so, file.path("/servers/sutherland-scratch/andi/projects/0_2306_Floris_Rawlings/data/R_objects", "6w_epithelial.rds"))

so <- subset(sample_so$"5891STDY9030809", stable_clusters %in% c(1,2,4,6,9))
so <- SCTransform(so, variable.features.n = 5000, return.only.var.genes = F, verbose = F)
so$stable_clusters <- NULL
saveRDS(so, file.path("/servers/sutherland-scratch/andi/projects/0_2306_Floris_Rawlings/data/R_objects", "8w_epithelial.rds"))
