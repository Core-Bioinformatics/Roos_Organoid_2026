library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)
library(qs)
library(harmony)

project_folder <- ""
ncores <- 30
objects_folder <- file.path(project_folder, "R_objects")
ca_folder <- file.path(objects_folder, "clustassess")
mtd_folder <- file.path(project_folder, "metadata")
if (!dir.exists(ca_folder)) {
    dir.create(ca_folder, recursive = TRUE)
}
shiny_ca_folder <- file.path(project_folder, "output", "clustassess_shiny_apps")
if (!dir.exists(shiny_ca_folder)) {
    dir.create(shiny_ca_folder, recursive = TRUE)
}

chosen_batch_correction <- "cca"
chosen_theta <- 10
choose_processing_function <- function(proc_type = "default", dt_mtx, npcs = 30, categ = NULL, ...) {
    if (proc_type == "default") {
        return(
            function(dt_mtx, actual_npcs = 30) {
                actual_npcs <- min(actual_npcs, ncol(dt_mtx)%/%2)
                
                RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
                embedding <- stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x
                
                RhpcBLASctl::blas_set_num_threads(1)
                rownames(embedding) <- rownames(dt_mtx)
            
                colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))
                
                return(embedding)
            }
        )
    }

    if (proc_type == "harmony") {
        return(
            function(dt_mtx, actual_npcs = 30) {
                actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)

                RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
                embedding <- stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x

                RhpcBLASctl::blas_set_num_threads(1)
                rownames(embedding) <- rownames(dt_mtx)
                colnames(embedding) <- paste0("PC_", seq_len(actual_npcs))

                embedding <- harmony::RunHarmony(embedding, categ, verbose = FALSE, ...)

                return(embedding)
            }
        )
    }

    if (proc_type == "cca") {
        return(
            function(dt_mtx, actual_npcs = 30) {
                actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)
                dt_mtx <- t(dt_mtx)
                # convert it to a normalised matrix from scaled
                print(dim(dt_mtx))
                sample_names <- unique(categ)
                print(length(categ))

                so_objects <- lapply(sample_names, function(sample_name) {
                    
                    ClustAssess::create_seurat_object_default(
                        normalized_expression_matrix = expm1(dt_mtx[, categ == sample_name])
                    )
                })
                names(so_objects) <- sample_names

                unified_so <- merge(
                    x = so_objects[[1]],
                    y = so_objects[-1],
                    add.cell.ids = sample_names
                )

                unified_so <- NormalizeData(unified_so, verbose = FALSE)
                unified_so <- FindVariableFeatures(unified_so, selection.method = "vst", nfeatures = nrow(dt_mtx), verbose = FALSE)
                unified_so <- ScaleData(unified_so, features = rownames(dt_mtx), verbose = FALSE)
                unified_so@assays$RNA@layers$scale.data <- dt_mtx

                RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
                unified_so <- RunPCA(unified_so, features = rownames(dt_mtx), npcs = actual_npcs, verbose = FALSE, approx = FALSE)

                unified_so <- IntegrateLayers(
                    object = unified_so,
                    method = CCAIntegration,
                    orig.reduction = "pca",
                    new.reduction = "cca",
                    scale.layer = "scale.data",
                    features = rownames(dt_mtx),
                    verbose = FALSE,
                    k.weight = min(100, min(sapply(so_objects, ncol))),
                    ...
                )
                RhpcBLASctl::blas_set_num_threads(1)

                return(unified_so@reductions$cca@cell.embeddings)
            }
        )
    }
}



mtd_configs <- read.csv(file.path(mtd_folder, "26_primary_clustassess_configs.csv"), header = TRUE, comment.char = "#")
nreps <- 100
neigh_seq <- seq(from = 5, to = 50, by = 5)
res_seq <- seq(from = 0.1, to = 2, by = 0.1)
assay_name <- "RNA"

for (i in seq_len(nrow(mtd_configs))) {
    id <- mtd_configs$id[i]
    so_path <- file.path(objects_folder, paste0(id, "_so.qs"))
    ca_path <- file.path(ca_folder, paste0(id, ".qs"))
    used_id <- id

    print(id)

    so <- qread(so_path, nthreads = ncores)
    DefaultAssay(so) <- assay_name

    # initialise the parallel context
    RhpcBLASctl::blas_set_num_threads(1)
    my_cluster <- parallel::makeCluster(
        ncores,
        type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = my_cluster)

    # extract the matrix and the features
    expr_matrix <- GetAssayData(so, assay = assay_name, layer = "scale.data")
    features <- rownames(expr_matrix)
    var_features <- VariableFeatures(so)
    max_ngenes <- length(var_features)
    if (max_ngenes > 2000) {
        fstep <- 500
    } else {
        fstep <- 250
    }

    most_abundant_genes <- rownames(expr_matrix)[order(Matrix::rowSums(expr_matrix), decreasing = TRUE)]

    gene_list <- list(
        "Most_Abundant" = most_abundant_genes[seq_len(max_ngenes)],
        "Highly_Variable" = var_features[seq_len(max_ngenes)]
    )

    steps_list <- list(
        "Most_Abundant" = seq(from = 500, by = fstep, to = max_ngenes),
        "Highly_Variable" = seq(from = 500, by = fstep, to = max_ngenes)
    )


    mtd_categ <- NULL
    mtd_harmony <- mtd_configs$mtd_harmony[i]
    if (is.na(mtd_harmony)) {
        mtd_harmony <- ""
    }
    if (mtd_harmony != "") {
        mtd_categ <- so@meta.data[, mtd_harmony]
        args_list <- list(
            proc_type = chosen_batch_correction,
            dt_mtx = expr_matrix,
            categ = mtd_categ
        )

        if (chosen_batch_correction == "harmony") {
            args_list$theta <- chosen_theta
        }

        processing_function <- do.call(choose_processing_function, args_list)
        used_id <- paste0(id, "_", chosen_batch_correction)
        if (chosen_batch_correction == "harmony") {
            used_id <- paste0(used_id, "_", gsub("\\.", "", chosen_theta))
        }
        ca_path <- file.path(ca_folder, paste0(used_id, ".qs"))
    } else {
        processing_function <- choose_processing_function(
            proc_type = "default",
            dt_mtx = expr_matrix
        )
    }

    if (file.exists(ca_path) && file.size(ca_path) > 0) {
        test_automm <- qs::qread(ca_path, nthreads = ncores)
    } else {
        rm(so)
        gc()
        test_automm <- automatic_stability_assessment(
            expression_matrix = expr_matrix,
            n_repetitions = nreps,
            temp_file = "clustassess_temp.rds",
            save_temp = FALSE,
            n_neigh_sequence = neigh_seq,
            resolution_sequence = res_seq,
            features_sets = gene_list,
            steps = steps_list,
            n_top_configs = 2,
            umap_arguments = list(
                min_dist = 0.3,
                n_neighbors = 30,
                metric = "cosine"
            ),
            matrix_processing = processing_function,
            verbose = TRUE
        )
        parallel::stopCluster(cl = my_cluster)
        qsave(test_automm, ca_path, nthreads = ncores)

        so <- qread(so_path, nthreads = ncores)
        DefaultAssay(so) <- assay_name
    }

    mtd_cols <- setdiff(colnames(so@meta.data), "cell_id")
    so@meta.data <- so@meta.data[, mtd_cols]
    app_path <- file.path(shiny_ca_folder, used_id)
    
    if (file.exists(file.path(app_path, "expression.h5")) && file.size(file.path(app_path, "expression.h5")) > 0) {
        next
    }

    so <- JoinLayers(so)
    write_shiny_app(
        object = so,
        assay_name = assay_name,
        clustassess_object = test_automm,
        shiny_app_title = used_id,
        project_folder = app_path,
        prompt_feature_choice = FALSE
    )
    gc()
}



