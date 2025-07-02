# Roos_Organoid_2026
Data analysis scripts for the Organoid project.

Manuscript link: [TO BE ADDED]


Short description of scripts
- bulk
    - `1_noisyr.Rmd` - Performs denoising of the bulk datasets using the [noisyR package](https://github.com/Core-Bioinformatics/noisyR). The denoising is done per sample at transcript level.
    - `2_bulkanalyser.Rmd` - Performs modifications of the metadata of the samples and creates, based on the metadata and the denoised expression matrix, a [bulkAnalyseR shiny app](https://github.com/Core-Bioinformatics/bulkAnalyseR).
- single_cell_nuclei
    - `1_create_seurat_object.R` - Create the Seurat object for each sample from the raw h5 CellRanger output. The samples are further merged together, where QC metrics are assessed.
    - `2_run_clustassess.R` - Run the automatic stability assessment pipeline as described in the [ClustAssess package](https://github.com/Core-Bioinformatics/ClustAssess).
    - `3_create_clustassess_shiny_app.R` - Creates the ClustAssess shiny app that is built upon the normalised matrix and the output of the stability assessment.
    - `4/5_create_pseudotime_first/second_app.R` - Creates the pseudotime apps that are used to infer the pseudotime trajectory by selecting the root and end points (first app) and identify the driver genes that define the transition between states by showing the distribution of modules (second app). We recommend using the [Starlng framework](https://github.com/Core-Bioinformatics/Starlng) as it contains up-to-date analysis and functionalities for these tasks.
    - `6_create_pseudotime_gene_heatmaps.R` - Visually describes the functionality of the gene modules in transitioning between cell states by creating cell-resolution heatmaps. This functionality can be recapitulated in the Starlng package.
    - `7_create_data_for_python.R` - Translate the information about datasets: normalised expression matrices, metadata dataframe, the set of highly variable genes, the dimensionality reduction embeddings into python-readable formats. These will be used in the RNA velocity analysis.
    - `8_run_velocity.sh` - Runs a Python script taylored based on the documentation of the [scvelo framework](https://github.com/theislab/scvelo).
    - `9_enrichment_analysis_gene_modules.R` - Enables enrichment analysis of the modules identified at steps 4,5,6 to describe functionally the cell groups. This feature is found in the Starlng framework as well.
- integration
    - `1_create_seurat_object.R` - Create the Seurat object for each sample used from the Rawlins dataset using the h5 output from CellRanger. Perform filtering and merge in one single Seurat object. Apply some basic preprocessing steps from the Seurat pipeline.
    - `2_split_objects_by_metadata.R` - Split the merged Seurat object in two independent datasets: 6w epithelial and 8w epithelial.
    - `3_run_clustassess.R` - Run the automatic stability assessment pipeline as described in the [ClustAssess package](https://github.com/Core-Bioinformatics/ClustAssess)
    - `4_write_shiny_ca.R` - Use the ClustAssess output to generate a ClustAssess shiny app for the two datasets.
    - `integration_scmap_helper.R` - Helper functions to run the integration on the Rawlins dataset with the Organoid. The integration is on UMAP level and consists of combining the nearest neigbhours identified at dataset level with the neighbours that can be identified using the [scmap package](https://github.com/hemberg-lab/scmap). The quality is assessed using stability functions similar with the ones from the ClustAssess pipeline.
    - `5_integration_scmap_assessment.R` - Run the integration quality assessment between r6w - o8w, and r8w - o12w respectively. The assessment is performed by varying the number of nearest neighbours.
    - `6_integration_markers.R` - Perform intersections between the dataset-specific stable clusters with the stable clusters identified on the integrated UMAP. For the intersection we find the Differentially expressed genes between the intersected set of cells and the rest of the cells. The DEG analysis is performed inside each dataset, since the integration is not performed on the expression level as well.
    - `7_enrichment_analysis_on_modules.R` - Performs enrichment analysis on each cluster from the gene modules identified in the pseudotime app. The pseudotime app is generated using a pipeline that is yet to be published, therefore it is not included in this repository.
    - `8_compare_gene_modules_jsi.R` - Compare the gene modules using the Jaccard Similarity Index (JSI) between Organoid and Rawlins dataset: r6w - o8w, r8w - o12w.
    - `9_compare_gene_modules_enrichment.R` - Compare the datasets by calculating the JSI between the enriched terms identified in the gene modules from the Organoid and Rawlins dataset: r6w - o8w, r8w - o12w.






*Note*: The scripts contain some variables that were not defined, such as `project_folder`. This is to mantain the confidentiality of the environment where the analysis was performed. The user should attempt to replicate the folder structure, as defined in the scripts.