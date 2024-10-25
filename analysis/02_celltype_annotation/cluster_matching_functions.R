suppressPackageStartupMessages({
library(tidyverse)
library(tidyseurat)
library(foreach)
library(doParallel)
library(SingleR)
library(parallel) 
library(cowplot)
library(ggdendro)
  library(Seurat)
})

### 1) initial clustering ###################

# input: list of Seurat objects
# parameters: resolution (make cell number dependent?)
# output: 
#     - seurat objects
#     - plot: UMAP with cluster labels

run_initial_clustering <- function(seu_file_list, out_dir, assay = "RNA", reduction = "pca", resolution = "CellNumber", force_rerun = F){
  if (dir.exists(out_dir) & force_rerun == F) {
    message(paste0("Directory already exists: ", out_dir, " - clustering is not executed"))
  } else {
    message(paste0("Saving clustering output in ", out_dir))
    dir.create(out_dir, recursive = T)
  
    registerDoParallel(cores = 12)
    foreach(name = names(seu_file_list), .combine = "c" ) %dopar% {
        print(paste0("Load ",name," input"))
        seu <- readRDS(seu_file_list[[name]])
        
        if(class(seu[["RNA"]]) == "Assay5"){
          print("Changing Assay type")
          seu[["RNA"]] <- as(seu[["RNA"]], Class = "Assay")
        }
        
        DefaultAssay(seu) <- assay
        
        # determine cluster resolution
        if(resolution == "CellNumber"){
          ncell <- ncol(seu)
          reso <- case_when(ncell < 1000 ~ 0.7, 
                                (1000 <= ncell & ncell < 2000) ~ 1, 
                                (2000 <= ncell & ncell < 5000) ~ 1.3, 
                                (5000 <= ncell & ncell < 10000) ~ 1.7, 
                                ncell >= 10000 ~ 2)
        } else {
          reso <- resolution
        }
        
        print(paste0("Performing initial clustering for ", name, " with resolution ", reso))
        seu <- ScaleData(seu, verbose = F) %>% 
          FindVariableFeatures(nfeatures = 2000, verbose = F) %>% 
          RunPCA(verbose = F) %>%
          RunUMAP(dims = 1:20,verbose = F, reduction = reduction) %>%
          FindNeighbors(dims = 1:20,verbose = F, reduction = reduction) %>%
          Seurat::FindClusters(resolution = reso,verbose = F)
        
        seu$spec_clust <- paste(name, seu$seurat_clusters, sep = "__")
        
        print(paste0("Finished for ",name," - saving result and UMAP"))
        p <- DimPlot(seu, group.by = "spec_clust", reduction = "umap")+theme(legend.position = "bottom")+guides(color = guide_legend(ncol = 3, override.aes = list(size = 2)))+ggtitle(name)
        
        ggsave(paste0(out_dir,"/UMAP_",name,".pdf"), plot = p, width = 5, height = 6)
        #saveRDS(seu, paste0(out_dir,"/seu_",name,".RDS"))
        
        init_clust <- data.frame(spec_clust = seu$spec_clust, row.names = colnames(seu))
        saveRDS(init_clust, paste0(out_dir,"/initial_clusters_",name,".RDS"))
        
    }
    stopImplicitCluster()
  }
}
  
  
### 2) reciprocal classification ################

# input: list of Seurat objects with seurat clusters
# output: 
#     - dataframe of pairwise cluster similarity

## helper functions ###
cross_validation <- function(seurat_obj, label_name){
  print("Split 80:20 and perform cross validation")
  # split species 80:20 into a training and test set
  training_cells <- seurat_obj@meta.data %>% 
    rownames_to_column("cellID") %>% 
    group_split(!!label_name) %>% map(~ sample(.x$cellID, size = floor(0.8 * nrow(.x)), replace = FALSE)) %>%
    flatten() %>%
    unlist() 
  
  train_obj <- seurat_obj[,training_cells]
  test_obj <- seurat_obj[,!(colnames(seurat_obj) %in% training_cells)]
  
  model <- trainSingleR(ref = train_obj@assays$RNA@data, labels = train_obj@meta.data[[label_name]], aggr.ref=TRUE)
  pred <- classifySingleR(test_obj@assays$RNA@data, model, assay.type=1)
  
  df <- data.frame(species.x = test_obj$species,
                   cluster.x = test_obj@meta.data[[label_name]],
                   species.y = test_obj$species,
                   cluster.y = pred$labels) %>% 
    group_by(species.x, cluster.x,species.y, cluster.y) %>% 
    summarize(n = n()) %>% 
    group_by(species.x, cluster.x) %>% 
    mutate(similarity_value = n/sum(n),
           metric = "classification")
  
  return(df)
}
perform_classification <- function(seu_list, model_list, spec1, spec2, label_name){
  print(paste0("Classification of ",spec1," with ",spec2," as reference"))
  pred <- classifySingleR(seu_list[[spec1]]@assays$RNA@data, 
                          model_list[[spec2]], 
                          assay.type=1)
  
  df <- data.frame(species.x = spec1,
                   cluster.x = seu_list[[spec1]]@meta.data[[label_name]],
                   species.y = spec2,
                   cluster.y = pred$labels) %>% 
    group_by(species.x, cluster.x,species.y, cluster.y) %>% 
    summarize(n = n()) %>% 
    group_by(species.x, cluster.x) %>% 
    mutate(similarity_value = n/sum(n),
           metric = "classification")
  
  return(df)
}

## full function ###
run_reciprocal_classification <- function(seu_file_list,init_clust_dir, outfile, label_name = "spec_clust", force_rerun = F){
  if (file.exists(outfile) & force_rerun == F) {
    message("Output file already exists, skipping reciprocal classification")
    return(NULL)
  }
  
  print("Read in seurat object")
  seu_list <- mclapply(seu_file_list, readRDS, mc.cores = 4)
  
  print("Add initial clustering results")
  for(spec in names(seu_list)){
    init_clusters <- readRDS(paste0(init_clust_dir,"/initial_clusters_",spec,".RDS"))
    seu_list[[spec]] <- AddMetaData(seu_list[[spec]], init_clusters)
  }
  
  
  # Change Assay to RNA and convert seurat v5 to earlier version if necessary
  seu_list <- lapply(seu_list, function(seu){
    DefaultAssay(seu) <- "RNA"
    if(class(seu[["RNA"]]) == "Assay5"){seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay")}
    return(seu)
    })
  
  # Find intersection of genes
  common_genes <- Reduce(intersect, lapply(seu_list, rownames))
  seu_list <- lapply(seu_list, function(seu){seu[common_genes,]}) 

  # Perform cross validation for each species to get within-species similarity
  self_class <- mclapply(seu_list, function(seu){cross_validation(seu,"spec_clust")}, mc.cores = length(seu_list))
  print(head(bind_rows(self_class)))
  
  # Train model for each species
  model_list <- mclapply(seu_list, function(ref){
    print("Training model")
    trainSingleR(ref = ref@assays$RNA@data, labels =ref@meta.data[[label_name]], aggr.ref=TRUE)
  }, mc.cores = length(seu_list))
  
  species_pairs <- combn(names(seu_list), 2, simplify = F)
  reciprocal_class <-  mclapply(species_pairs, function(x, seu_list, model_list, label_name){
    rbind(perform_classification(seu_list, model_list, x[1],x[2], label_name),
          perform_classification(seu_list, model_list, x[2],x[1], label_name))    # reverse classification
  }, seu_list, model_list, label_name, mc.cores = 12)
  
  
  # combine everything
  print("Combine within and across species classification")
  res <- rbind(bind_rows(reciprocal_class), bind_rows(self_class))
  
  
  saveRDS(list(similarity = res,
               seu_metadata = lapply(seu_list, function(x)x@meta.data)),
          outfile)
}


### 3) create distance matrix ###################
# input: dataframe of pairwise cluster similarity
# parameters: summarization method
# output: distance matrix

convert_RC_to_dist <- function(RC_output, summarize_by = "mean", outfile, force_rerun = F){
  if (file.exists(outfile) & force_rerun == F) {
    message("Distance matrix output file already exists, not creating a new one")
    return(NULL)
  }
  
  mat <- readRDS(RC_output)$similarity %>% 
    ungroup() %>% 
    #mutate(X = paste(species.x,cluster.x, sep = "__"),
    #       Y = paste(species.y,cluster.y, sep = "__")) %>% 
    select(cluster.x,cluster.y,similarity_value) %>% 
    pivot_wider(names_from = cluster.x, values_from = similarity_value) %>% 
    column_to_rownames("cluster.y")
  
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
  
  d <- as.dist(1-mat_sym)
  saveRDS(d, outfile)
}


### 4) Hierarchical clustering and tree cutting #########
cluster_on_distmat <- function(distmat, method = "fixed_k", k = NULL, resolution = NULL, outfile, force_rerun = F){
  if (file.exists(outfile) & force_rerun == F) {
    message("Initial consensus output file already exists, not creating a new one")
    return(NULL)
  }
  
  d <- readRDS(distmat)
  
  if(method == "fixed_k"){
    hc <- hclust(d, method = "average")
    cut_k <- cutree(hc, k = k)
    cut_join <- data.frame(spec_clust = names(cut_k), consensus = as.numeric(cut_k))
  }
  if(method == "dynamic"){
    hc <- hclust(d, method = "average")
    ctD <- cutreeDynamic(hc, minClusterSize = 1, distM = as.matrix(d))
    cut_join <- data.frame(spec_clust = hc$labels, consensus = as.numeric(ctD))
  }
  if(method == "louvain"){
    sim_mat <- 1 - as.matrix(d)
    graph <- graph_from_adjacency_matrix(sim_mat, mode = "undirected", diag = F, weighted = T)
    clusters <- cluster_louvain(graph, resolution = resolution)
    cluster_membership <- membership(clusters)
    cut_join <- data.frame(spec_clust = names(cluster_membership), consensus = as.vector(cluster_membership))
  }
  
  res <- cut_join %>% 
    #mutate(species = gsub("__.*","",spec_clust)) %>% 
    group_by(consensus) %>% 
    #mutate(nSpecies = length(unique(species))) %>% #filter(nSpecies > 2) %>% 
    ungroup()
  
  saveRDS(res, outfile)
}

### 5) Merge similar clusters ######
run_cluster_merging <- function(seu_file_list,init_clust_dir, orth_clusters_out, cut_height = 0.1, outfile, force_rerun = F){
  if (file.exists(outfile) & force_rerun == F) {
    message("Merged clusters output file already exists, not creating a new one")
    return(NULL)
  }
  
  print("-------------------------------------Running cluster merging -------------------------------------------")
  seu_list <- mclapply(seu_file_list, readRDS, mc.cores = 4)
  tree_cut <- readRDS(orth_clusters_out)
  
  print("Add consensus clusters")
  for(spec in names(seu_list)){
    init_clusters <- readRDS(paste0(init_clust_dir,"/initial_clusters_",spec,".RDS")) %>% 
      left_join(tree_cut)
    seu_list[[spec]] <- AddMetaData(seu_list[[spec]], init_clusters)
  }
  
  #seu_list <- lapply(seu_list, AddMetaData_byJoin,tree_cut)
  
  print("Filter for highly variable genes")
  seu_list <- mclapply(seu_list,  FindVariableFeatures, nfeatures = 2000, mc.cores = 4)
  HVG <- Reduce(union, lapply(seu_list, VariableFeatures))
  seu_list <- lapply(seu_list, function(x)x[HVG,])
  
  #print(lapply(seu_list, function(x){head(x@meta.data)}))
  
  print("Create pseudobulk per consensus cluster per species")
  pseudobulk_list <- mclapply(seu_list, function(seu){
    if(class(seu[["RNA"]]) == "Assay5"){seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay")}
    AverageExpression(seu, group.by = "consensus")$RNA
  }, mc.cores = 4)
  
  print(head(pseudobulk_list$human))
  
  print("Spearman correlation")
  cor_df <- mclapply(pseudobulk_list, function(pb){
    cor(as.matrix(pb), method = "spearman") %>% data.frame %>%
      rownames_to_column("consensus1") %>% 
      pivot_longer(cols = -consensus1, names_to = "consensus2", values_to = "spearman") %>% 
      mutate(consensus2 = gsub("X","",consensus2))
  }, mc.cores = 4) %>% 
    bind_rows(.id = "species")
  
  print(head(cor_df))
  
  cor_mat <- cor_df %>% 
    group_by(consensus1, consensus2) %>% 
    summarize(mean_spearman = mean(spearman)) %>% 
    pivot_wider(names_from = consensus2, values_from = mean_spearman) %>% 
    column_to_rownames("consensus1")
  cor_mat <- cor_mat[,rownames(cor_mat)]
  d_pb <- as.dist(cor_mat)
  
  print("done^2")
  
  print(paste0("Hierarchical clustering and tree cutting at height ",cut_height))
  d_pb[is.na(d_pb)] <- min(d_pb, na.rm = T)
  
  cut <- cutree(tree = hclust(1-d_pb), h = cut_height)
  cut_join_merged <- data.frame(consensus = as.numeric(names(cut)), consensus_merged = cut) %>% 
    right_join(select(tree_cut, c(consensus, spec_clust)))
  
  metadata <- bind_rows(lapply(seu_list, function(x)x@meta.data)) %>% 
    rownames_to_column("BC") %>% 
    select(BC, spec_clust) %>% 
    left_join(cut_join_merged) %>% 
    column_to_rownames("BC")
  
  saveRDS(metadata, outfile)
  
}

### 6) Plotting #######
# *** Distance matrix and merging ####
plot_distmap_and_hc <- function(d, merged_consensus, cluster_colors, species_colors = color_schemes$species){
  # hierarchical clustering
  hc <- hclust(d, method = "average")
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc, type = "rectangle")
  
  label_data <- ddata$labels %>% 
    mutate(species = word(label,1,1,"__")) %>% 
    left_join(distinct(merged_consensus), by = c("label" = "spec_clust"))
  
  plot_annotation <- function(label_data, color_by){
    ggplot(label_data, aes_string(x="x", y="1", fill = paste0("as.factor(",color_by,")")))+
      geom_bar(stat = "identity")+
      theme_nothing()+
      theme(legend.position = "none", axis.text = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 8))+
      labs(y = color_by)
  }
  
  p.dendro <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    theme_nothing()
  
  p.mat <- reshape2::melt(as.matrix(d)) %>% 
    mutate(Var1 = factor(Var1, levels = rownames(as.matrix(d))[hc$order]),
           Var2 = factor(Var2, levels = rownames(as.matrix(d))[hc$order]),
           value = ifelse(Var1 == Var2, NA,value)) %>% 
    ggplot(aes(x = as.numeric(Var1), y = Var2, fill = 1-value))+
    geom_tile(color = "black")+
    scale_fill_gradient2(midpoint = 0.01)+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          axis.text.y = element_text(size = 4),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  p <- plot_grid(p.dendro,
                 plot_annotation(label_data, "consensus"),      
                 plot_annotation(label_data, "consensus_merged")+scale_fill_manual(values = cluster_colors),
                 #plot_annotation(label_data2, "celltype_assigned")+scale_fill_manual(values = color_schemes$celltype),
                 plot_annotation(label_data, "species")+scale_fill_manual(values = species_colors),
                 p.mat,
                 ncol = 1, align = "v", axis = "lr", rel_heights = c(0.5,0.1,0.1,0.1,1))
  
  return(p)
}

# *** Consensus_cluster results #######
plot_consensus_clusters <- function(metadata_combined, cluster_colors, species_colors = color_schemes$species){
  p_inp <- metadata_combined %>% 
    group_by(species,consensus_merged) %>% 
    dplyr::summarise(nCells = n()) %>% 
   # distinct() %>% 
    group_by(consensus_merged) %>% 
    mutate(nSpecies = length(unique(species)),
           consensus_merged = as.factor(consensus_merged))
  
  p1 <-  ggplot(p_inp, aes(x = consensus_merged, y = species, size = nCells, fill = consensus_merged))+
    geom_point(shape = 21)+
    facet_grid(~nSpecies, space = "free", scales = "free")+
    scale_fill_manual(values = cluster_colors, limits = force)+
    theme_bw()+
    theme(axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
    )+
    labs(x = "consensus cluster")
  
  p2 <- p_inp %>% 
    group_by(species, nSpecies) %>% 
    summarise(total_cells = sum(nCells)) %>% 
    
    ggplot(aes(x = factor(nSpecies, levels = 1:4), y = total_cells, fill = species))+
    geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
    scale_fill_manual(values = species_colors, limits = force)+
    theme_bw()+
    labs(x = "species specificity level", y = "number of cells")
  
  p <- plot_grid(p1,p2, nrow = 1, align = "h", rel_widths = c(1,0.5), axis = "b")
  
  return(p)
}

plot_consensus_UMAPs <- function(seu_list, cluster_colors){
  plot_grid(plotlist = lapply(seu_list, function(seu){DimPlot(seu, group.by = "consensus_merged", label = T, reduction = "umap")+
    scale_color_manual(values = cluster_colors, limits = force)+ggtitle(unique(seu$species))}), 
    nrow = 2)
}

plot_final_clusters <- function(seu_file_list, distmat, merged_consensus, outfile){
  if (file.exists(outfile)) {
    message("Plotting output file already exists, not creating a new one")
    return(NULL)
  }
  
  print("Load output files for plotting")
  seu_list <- mclapply(seu_file_list, readRDS, mc.cores = 4)
  cut_join <- readRDS(merged_consensus)
  d <- readRDS(distmat)
  
  #seu_list <- lapply(seu_list, function(x)AddMetaData_byJoin(x, cut_join))
  seu_list <- lapply(seu_list, function(x)AddMetaData(x, cut_join))
  metadata_all <- lapply(seu_list, function(x)x@meta.data) %>% bind_rows(.id = "species")
  
  # define color scheme
  cluster_colors <- setNames(scales::hue_pal()(length(unique(metadata_all$consensus_merged))),
                             unique(metadata_all$consensus_merged))
  
  p <- plot_grid(
    plot_distmap_and_hc(d, cut_join, cluster_colors),
    plot_grid(
      plot_consensus_clusters(metadata_all, cluster_colors), 
      plot_consensus_UMAPs(seu_list, cluster_colors), 
      nrow = 2, rel_heights = c(0.3,1)
    ),
    nrow = 1)
  
  ggsave(outfile, plot = p, width = 20, height = 10)
}

#seu_integrated <- readRDS("/data/home/EBgrant/analysis_scRNA_allRuns/files/seu_ectoderm_integrated_Seurat5.RDS")
#merged_consensus <- "/data/home/EBgrant/analysis_scRNA_allRuns/files/cluster_matching/ectoderm/res_CellNumber/integration_consensus_clusters_harmony_res0.2.RDS"


plot_integrated_UMAP <- function(seu_integrated, reduction, color_by, color_palette){
  df <- seu_integrated@reductions[[reduction]]@cell.embeddings
  colnames(df) <- c("UMAP1", "UMAP2")
  p <- df %>% data.frame %>% 
    mutate(color = seu_integrated[[color_by]][[1]]) %>% 
    ggplot(aes(x = UMAP1, y = UMAP2, color = color))+
      geom_point(size = 0.3)+
      scale_color_manual(values = color_palette)+
      theme_bw()+
      theme(legend.title = element_blank())
  return(p)
}

plot_final_clusters_integrated_vs_individual <- function(seu_list, seu_integrated, merged_consensus, outfile){
  cut_join <- readRDS(merged_consensus) %>% mutate(consensus_merged = as.factor(consensus_merged))
  seu_list <- lapply(seu_list, function(x)AddMetaData(x, cut_join))
  metadata_all <- lapply(seu_list, function(x)x@meta.data) %>% bind_rows(.id = "species")
  seu_integrated <- AddMetaData(seu_integrated, cut_join)
  
  # define color scheme
  cluster_colors <- setNames(scales::hue_pal()(length(unique(metadata_all$consensus_merged))),
                             unique(metadata_all$consensus_merged))
  
  p <- plot_grid(
    plot_grid(
      plot_integrated_UMAP(seu_integrated, reduction = "umap.harmony", color_by = "consensus_merged", color_palette = cluster_colors), 
      plot_consensus_clusters(metadata_all, cluster_colors, species_colors = color_schemes$species), 
      nrow = 2, rel_heights = c(1,0.3)
    ),
    plot_consensus_UMAPs(seu_list, cluster_colors), 
    nrow = 1
  )
  
  ggsave(outfile, plot = p, width = 20, height = 10)
}


#### Run full pipeline #####
run_cluster_matching <- function(seu_file_list, out_dir, IC_resolution, summarize_RC_by, reduction = "pca", assay = "scanorama", tree_cut_k, merge_height){
  species <- names(seu_file_list)
  dir.create(out_dir)
  # Define interdependent output files of each step:
  initial_clustering_dir <- paste0(out_dir,"/res_", IC_resolution)
  RC_out <- paste0(initial_clustering_dir,"/reciprocal_classification.RDS")
  dist_out <- paste0(initial_clustering_dir,"/distmat_",summarize_RC_by,".RDS")
  initial_consensus_out <- paste0(initial_clustering_dir,"/initial_consensus_clusters_RC",summarize_RC_by,"_k",tree_cut_k,".RDS")
  merged_consensus_out <- paste0(initial_clustering_dir,"/merged_consensus_clusters_RC",summarize_RC_by,"_k",tree_cut_k,"_merge",merge_height,".RDS")
  plotting_out <- paste0(initial_clustering_dir,"/plots_RC",summarize_RC_by,"_k",tree_cut_k,"_merge",merge_height,".pdf")
  run_initial_clustering(seu_file_list,
                         out_dir = initial_clustering_dir,
                         reduction = reduction,
                         assay = assay,
                         resolution = IC_resolution)
  
  clustered_seu_list <- lapply(setNames(species,species), function(spec){
    paste0(initial_clustering_dir,"/seu_",spec,".RDS")
  })
  
  run_reciprocal_classification(clustered_seu_list, 
                                outfile = RC_out, 
                                label_name = "spec_clust")
  
  convert_RC_to_dist(RC_out = RC_out,
                     summarize_by = summarize_RC_by,
                     outfile = dist_out)
  
  cluster_on_distmat(distmat = dist_out, 
                     method = "fixed_k",
                     k = tree_cut_k, 
                     outfile = initial_consensus_out,
                     force_rerun = F)
  
  run_cluster_merging(clustered_seu_list, 
                      initial_consensus = initial_consensus_out, 
                      cut_height = merge_height, 
                      outfile = merged_consensus_out,
                      force_rerun = F)
  
  plot_final_clusters(clustered_seu_list, 
                      distmat = dist_out,
                      merged_consensus = merged_consensus_out, 
                      outfile = plotting_out)
}
