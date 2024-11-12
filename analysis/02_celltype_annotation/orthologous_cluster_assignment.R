# Identify orthologous clusters across species
setwd("/data/share/htp/EBgrant/analysis_scRNA_allRuns/publication_scripts")

source("EB-analyses/analysis/02_celltype_annotation/cluster_matching_functions.R")

output_directory <- "zenodo/cell_type_assignment/intermediate_files"

# specify the input: one seurat object per species
species <- c("human","cyno","rhesus","orang")
seu_file_list <- lapply(setNames(species,species), function(spec){
  paste0("zenodo/processing_per_species/seu_",spec,".RDS")
})

# 1) High resolution clustering per species
# Clustering resolution set based on total cell number per species
# Scanorama integrated embedding as input to correct batch effects from different experiments
run_initial_clustering(seu_file_list = seu_file_list,
                       out_dir = output_directory,
                       assay = "RNA",
                       reduction = "pca_scanorama",
                       resolution = "CellNumber")

# 2) Reciprocal classification
# Calculate within and across species similarity scores between clusters
run_reciprocal_classification(seu_file_list = seu_file_list, #lapply(setNames(species,species), function(spec){paste0(output_directory, "/seu_",spec,".RDS")}),
                              init_clust_dir = output_directory,
                              outfile = paste0(output_directory, "/reciprocal_classification.RDS"),
                              label_name = "spec_clust")

# 3) Create distance matrix
# Convert classification based cluster similarity scores into a distance matrix
convert_RC_to_dist(RC_output = paste0(output_directory, "/reciprocal_classification.RDS"),
                   summarize_by = "mean",
                   outfile = paste0(output_directory,"/distmat_mean.RDS"))

# 4) Identify orthologous clusters
# Hierarchical clustering on the cluster distance matrix to find groups of clusters ("orthologous clusters") that are most similar to eachother
# Different parameters for this step can be explored in the shiny app and the choice of k and clustering method is based on this
cluster_on_distmat(distmat = paste0(output_directory,"/distmat_mean.RDS"), 
                   method = "fixed_k",
                   k = 26, 
                   outfile = paste0(output_directory,"/orthologous_clusters.RDS"))

# 5) Merge similar orthologous clusters
# Combine orthologous clusters with a very similar pseudobulk expression profile
# Different parameter choices can be explored in the shiny app
run_cluster_merging(seu_file_list = seu_file_list,
                    init_clust_dir = output_directory,
                    orth_clusters_out = paste0(output_directory,"/orthologous_clusters.RDS"), 
                    cut_height = 0.1, 
                    outfile = paste0(output_directory,"/orthologous_clusters_merged.RDS"))


plot_final_clusters(clustered_seu_list, 
                    distmat = paste0(output_directory,"/distmat_mean.RDS"),
                    merged_consensus = merged_consensus_out, 
                    outfile = plotting_out)