library(tidyverse)
library(Seurat)
library(MetaNeighbor)
library(SummarizedExperiment)

seu_list <- readRDS("zenodo/cell_type_assignment/seu_list.RDS")


# Reciprocal classification ####
source("analysis/02_celltype_annotation/cluster_matching_functions.R")

# save separate objects per species
species <- c("human","cyno","rhesus","orang")
for(spec in species){
  saveRDS(seu_list[[spec]], paste0("input/seu_",spec,".RDS"))
}


seu_file_list <- lapply(setNames(species,species), function(spec){
  paste0("input/seu_",spec,".RDS")
})

run_reciprocal_classification(seu_file_list = seu_file_list,
                              init_clust_dir = output_directory,
                              outfile = "output/reciprocal_classication_manual_annotation.RDS",
                              label_name = "manual_annotation",
                              import_clustering_results = F)


# import results
convert_RC_to_dist <- function(RC_output, summarize_by = "mean"){
  mat <- readRDS(RC_output)$similarity %>% 
    ungroup() %>% 
    filter(cluster.x != "unassigned", cluster.y != "unassigned") %>% 
    mutate(X = paste(species.x,cluster.x, sep = "__"),
           Y = paste(species.y,cluster.y, sep = "__")) %>% 
    select(X,Y,similarity_value) %>% 
    pivot_wider(names_from = X, values_from = similarity_value) %>% 
    column_to_rownames("Y")
  
  mat <- mat[colnames(mat),]
  mat[is.na(mat)] <- 0
  
  if(summarize_by == "mean"){
    mat_sym <- (as.matrix(mat) + t(mat)) / 2
  }
  if(summarize_by == "max"){
    mat_sym <- pmax(as.matrix(mat), t(mat))
  }
  if(summarize_by == "min"){
    mat_sym <- pmin(as.matrix(mat), t(mat))
  }
  if(summarize_by == "directed"){
    mat_sym <- mat
  }
  
  #d <- as.dist(1-mat_sym)
  #return(d)
  
  return(mat_sym)
}

sim_RC_mean <- convert_RC_to_dist(RC_output = "output/reciprocal_classification_manual_annotation.RDS", summarize_by = "mean")
saveRDS(sim_RC_mean, "zenodo/cell_type_assignment/correspondence_validation/reciprocal_classification_out.RDS")

# MetaNeighbor ####

# Convert to SummarizedExperiment 
count_matrices <- list()
col_data_list <- list()

for (species_name in names(seu_list)) {
  seu_obj <- seu_list[[species_name]]
  
  # Extract the count matrix
  count_matrix <- GetAssayData(seu_obj, slot = "counts")
  
  # Extract the metadata for colData
  metadata <- seu_obj@meta.data
  col_data <- data.frame(
    sample_id = rownames(metadata),        # Cell barcodes
    study_id = species_name,               # Species (study_id)
    cell_type = metadata$manual_annotation # Cell type
  )
  
  # Add to the lists
  count_matrices[[species_name]] <- count_matrix
  col_data_list[[species_name]] <- col_data
}

# Combine the count matrices and colData from all species
combined_counts <- do.call(cbind, count_matrices)
combined_col_data <- do.call(rbind, col_data_list)
rownames(combined_col_data) <- combined_col_data$sample_id

# Create the SummarizedExperiment object
se_obj <- SummarizedExperiment(
  assays = list(counts = combined_counts),
  colData = combined_col_data
)

se_obj <- se_obj[,se_obj$cell_type != "unassigned"]

# Run MetaNeighbor
var_genes = variableGenes(dat = se_obj, exp_labels = se_obj$study_id)

celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = se_obj,
                             study_id = se_obj$study_id,
                             cell_type = se_obj$cell_type,
                             fast_version = T)

rownames(celltype_NV) <- gsub("\\|","__",rownames(celltype_NV))
colnames(celltype_NV) <- gsub("\\|","__",colnames(celltype_NV))

saveRDS(celltype_NV, "zenodo/cell_type_assignment/correspondence_validation/MetaNeighbor_out.RDS")