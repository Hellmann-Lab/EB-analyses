# color schemes ####
# species
spec_colors <-        c("#e9ecef","lightgrey", "#B0144F", "#92AEC8", "#3AA6A6", "#F2A518")
names(spec_colors) <- c("unassigned","doublet","human","rhesus","cynomolgus", "orang")

# individuals
ind_colors <-        c("#D81159", "#871744",  "#E59A10", "#FFBC42", "#218380","#A7BED3","blue" ,"lightgrey","#e9ecef")
names(ind_colors) <- c("human_29B5","human_63Ab2.2", "orang_69A1","orang_68A20","cyno_56A1", "rhesus_83D1", "rhesus_83Ab1.1", "doublet", "unassigned")

# cell types
ct_colors <- setNames(c("#a06cd5", "#c19ee0", "#dec9e9",
                        "#ddf5ff", "#90e0ef",  "#0096c7","#00b4d8" ,"#023e8a", "#023e8a","#03045e", "#03045e", 
                        "#68d8d6", "#07beb8",
                        "#F86363","#F52929","#BF0D0D","#830B0B", "#f1a7a9","#621708",
                        "#ffcad4", "#FF85A9","#E00043","#FF4782", "#ffe3e0", 
                        "#81c14b", "#1a7431", "lightgrey"),
                      nm= c("iPSCs_I", "iPSCs_II", "iPSCs_III", 
                            "early_ectoderm_I", "early_ectoderm_II", "early_ectoderm_III", "early_ectoderm_IV","neural_cells_I" ,"neural_cells_II" , "neurons_I", "neurons_II", 
                            "neural_crest_I", "neural_crest_II", 
                            "cardiac_I", "cardiac_II", "cardiac_III", "cardiac_IV", "cardiac_V", "cardiomyocytes",
                            "mesoderm_I","mesoderm_II","mesoderm_III","smooth_muscle_cells",  "early_mesoderm",
                            "endoderm","hepatocytes", "unassigned"))

rhodes_colors <- setNames(c("#c49be8", 
                            "#90e0ef", "#407ba7", "#023e8a", 
                            "#ffafcc","#c9184a",
                            "#1a7431"),
                          nm= c("Pluripotent_Cells",  
                                "Early_Ectoderm",  "Neural_Crest", "Neurons",
                                "Mesoderm", "Endothelial_Cells",
                                "Endoderm"))

color_schemes = list(species = spec_colors,
                     individuals = ind_colors,
                     celltype = ct_colors,
                     celltype_rhodes = rhodes_colors)


# functions ####
create_seurat_species <- function(mapped_cnts, species_name){
  filtered_assignment <- filter(species_assignment, species == species_name)
  seu_list <- list()
  for(experiment in names(mapped_cnts)){
    counts <- mapped_cnts[[experiment]]
    colnames(counts) <- paste(experiment,colnames(counts), sep = "_")  # add experiment ID to cell barcodes
    counts_species <- counts[, colnames(counts) %in% rownames(filtered_assignment)] # filter for species
    info <- filtered_assignment[colnames(counts_species),]
    
    if(nrow(info) > 0){   # not all species were used in all experiments - only return if there are cells
      seu <- CreateSeuratObject(counts = counts_species, meta.data = info)
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
      seu_list[[experiment]] <- seu
    }
    
  }
  return(seu_list)
}


plot.QC <- function(seu, mito_thresh = 8, gene_thresh=1000){
  #seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # plot QC metrics
  p.nUMI <- ggplot(seu, aes(x = experiment, y = nCount_RNA))+
    ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.2)+
    scale_color_manual(values = color_schemes$species)+
    geom_boxplot(alpha = 0.7, outlier.size = 0.4)+
    facet_grid(~species)+
    scale_y_log10()+
    theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(y = "nUMI/cell")
  
  p.nGene <- ggplot(seu, aes(x = experiment, y = nFeature_RNA))+
    ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.2)+
    scale_color_manual(values = color_schemes$species)+
    geom_boxplot(alpha = 0.7, outlier.size = 0.4)+
    facet_grid(~species)+
    scale_y_log10()+
    theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+
    geom_hline(yintercept = gene_thresh, color = "red")+
    labs(y = "nGene/cell")
  
  p.Mito <- ggplot(seu, aes(x = experiment, y = percent.mt))+
    ggbeeswarm::geom_quasirandom(aes(color = species), size = 0.2)+
    scale_color_manual(values = color_schemes$species)+
    geom_boxplot(alpha = 0.7, outlier.size = 0.4)+
    facet_grid(~species)+
    #scale_y_log10()+
    theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    geom_hline(yintercept = mito_thresh, color = "red")+
    labs(y = "% Mitochondrial / cell")
  
  plot_grid(p.nUMI, p.nGene, p.Mito, nrow = 3, align = "hv")
}

seurat_pipeline <- function(x, normalize = T, nVarGenes = 5000, nPC = 15, resolution = 0.5, umap = TRUE){
  if(normalize == T){x <- NormalizeData(x)}
  x <- FindVariableFeatures(x, nfeatures = nVarGenes)
  x <- ScaleData(x)
  # dim reduction
  x <- RunPCA(x)
  # clustering
  x <- FindNeighbors(x, dims = 1:nPC)
  x <- FindClusters(x, resolution = resolution)
  
  if(isTRUE(umap)){
    x <- RunUMAP(x, dims = 1:nPC)
  }
}

scDbl_pipeline <- function(seu, dbr_expected){
  sce <- as.SingleCellExperiment(seu)
  
  # Calculate doublet score for each cell
  print("--Calculating doublet score--")
  sce$doubletScore <- scDblFinder::computeDoubletDensity(sce)
  
  # Convert to doublet calls
  print("**Calling doulets**")
  dbl.calls <- scDblFinder::doubletThresholding(data.frame(score=sce$doubletScore),
                                                method="auto", returnType="call", dbr = dbr_expected)
  
  sce$doubletCall <- dbl.calls
  
  # Add doublet info  
  dbl_inf <- as.data.frame(colData(sce)) %>%
    select(c(doubletScore, doubletCall))
  seu <- AddMetaData(seu, metadata = dbl_inf)
  return(seu)
}

scran_norm <- function(seu){
  sce <- as.SingleCellExperiment(DietSeurat(seu))
  clust <- quickCluster(sce)
  sce <- computeSumFactors(sce, cluster = clust)
  sce <- logNormCounts(sce, log= FALSE, transform= "none")
  seu@assays$RNA@data <- as.matrix(log(x = assay(sce, "normcounts") + 1))
  return(seu)
}

scanorama_integration <- function(seu_list){
  dataset_list <- lapply(seu_list, function(x){
    t(GetAssayData(x, "data"))
  })
  
  genes_list <- lapply(dataset_list, colnames)
  
  # remove names from list
  names(dataset_list) <- NULL
  names(genes_list) <- NULL
  
  # perform integration
  scanorama <- import("scanorama")
  scanorama$integrate
  integrated <- scanorama$integrate(dataset_list, genes_list)
  corrected <- scanorama$correct(dataset_list, genes_list, return_dense = TRUE)
  
  return(list(integrated = integrated,
              corrected = corrected))
}

scanorama_to_seurat <- function(seu_list, scanorama_out, add_corrected_counts = F){
  # merge input Seurat object
  seu_merged <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
  
  # add integrated dimensional reductions to seurat object
  panorama <- do.call(rbind, scanorama_out$integrated[[1]])
  colnames(panorama) <- paste0("PC_", 1:100)
  rownames(panorama) <- colnames(seu_merged)
  
  seu_merged[["pca_scanorama"]] <- CreateDimReducObject(embeddings = panorama, key = "PC_")
  
  seu_merged <- RunUMAP(seu_merged, dims = 1:20, reduction = "pca_scanorama")
  
  
  if(add_corrected_counts){
    # merge scanorama output
    sc_list <- scanorama_out$corrected[[1]]
    names(sc_list) <- names(seu_list)
    scanorama_combined <- dplyr::bind_rows(lapply(sc_list, data.frame))
    # invert rows and columns
    scanorama_combined <- t(scanorama_combined)
    # add row and column names
    colnames(scanorama_combined) <- colnames(seu_merged)
    rownames(scanorama_combined) <- scanorama_out$corrected[[2]]
    
    # create seurat assay object and add scanorama output
    sc_assay <- CreateAssayObject(data = scanorama_combined)
    seu_merged[["scanorama"]] <- sc_assay
  }
  
  return(seu_merged)
}

HarmonyIntegrationMultiVar <- function (object, orig.reduction = "pca", groups, features = NULL, scale.layer = "scale.data", 
                                        layers = NULL, npcs = 50L, key = "harmony_", theta = NULL, 
                                        lambda = NULL, sigma = 0.1, nclust = NULL, tau = 0, block.size = 0.05, 
                                        max.iter.harmony = 10L, max.iter.cluster = 20L, epsilon.cluster = 1e-05, 
                                        epsilon.harmony = 1e-04, verbose = TRUE, ...) 
{
  orig <- object@reductions[[orig.reduction]]
  
  harmony.embed <- harmony::HarmonyMatrix(data_mat = Embeddings(object = orig), 
                                          meta_data = object@meta.data, vars_use = groups, do_pca = FALSE, 
                                          npcs = 0L, theta = theta, lambda = lambda, sigma = sigma, 
                                          nclust = nclust, tau = tau, block.size = block.size, 
                                          max.iter.harmony = max.iter.harmony, max.iter.cluster = max.iter.cluster, 
                                          epsilon.cluster = epsilon.cluster, epsilon.harmony = epsilon.harmony, 
                                          return_object = FALSE, verbose = verbose)
  rownames(x = harmony.embed) <- Cells(x = orig)
  dr <- suppressWarnings(expr = CreateDimReducObject(embeddings = harmony.embed, 
                                                     key = key, assay = DefaultAssay(object = orig)))
  
  object[["harmony"]] <- dr
  
  return(object)
}

HarmonyIntegration_seulist <- function(seu_list, groups){
  print("Combine into one object and separate out layers")
  obj <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
  DefaultAssay(obj) <- "RNA"
  
  print("Run basic pipeline")
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, features = rownames(obj), verbose = F)
  
  obj <- HarmonyIntegrationMultiVar(obj, groups = groups)
  
  return(obj)
}
