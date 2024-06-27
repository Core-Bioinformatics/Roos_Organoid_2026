# Roos_Organoid_2024
Data analysis scripts for the Organoid project.


Short description of scripts
- single_cell
- single_nuclei
- integration
    - `1_create_seurat_object.R` - Create the Seurat object for each sample used from the Rawlings dataset using the h5 output from CellRanger. Perform filtering and merge in one single Seurat object. Apply some basic preprocessing steps from the Seurat pipeline.
    - `2_split_objects_by_metadata.R` - Split the merged Seurat object in two independent datasets: 6w epithelial and 8w epithelial.
    - `3_run_clustassess.R` - Run the automatic stability assessment pipeline as described in the [ClustAssess package](https://github.com/Core-Bioinformatics/ClustAssess)
    - `4_write_shiny_ca.R` - Use the ClustAssess output to generate a ClustAssess shiny app for the two datasets.
    - `integration_scmap_helper.R` - Helper functions to run the integration on the Rawlings dataset with the Organoid. The integration is on UMAP level and consists of combining the nearest neigbhours identified at dataset level with the neighbours that can be identified using the [scmap package](https://github.com/hemberg-lab/scmap). The quality is assessed using stability functions similar with the ones from the ClustAssess pipeline.
    - `5_integration_scmap_assessment.R` - Run the integration quality assessment between r6w - o8w, and r8w - o12w respectively. The assessment is performed by varying the number of nearest neighbours.
    - `6_integration_markers.R` - Perform intersections between the dataset-specific stable clusters with the stable clusters identified on the integrated UMAP. For the intersection we find the Differentially expressed genes between the intersected set of cells and the rest of the cells. The DEG analysis is performed inside each dataset, since the integration is not performed on the expression level as well.
    - `7_enrichment_analysis_on_modules.R` - Performs enrichment analysis on each cluster from the gene modules identified in the pseudotime app. The pseudotime app is generated using a pipeline that is yet to be published, therefore it is not included in this repository.
    - `8_compare_gene_modules_jsi.R` - Compare the gene modules using the Jaccard Similarity Index (JSI) between Organoid and Rawlings dataset: r6w - o8w, r8w - o12w.
    - `9_compare_gene_modules_enrichment.R` - Compare the datasets by calculating the JSI between the enriched terms identified in the gene modules from the Organoid and Rawlings dataset: r6w - o8w, r8w - o12w.






*Note*: The scripts contain some variables that were not defined, such as `project_folder`. This is to mantain the confidentiality of the environment where the analysis was performed. The user should attempt to replicate the folder structure, as defined in the scripts.