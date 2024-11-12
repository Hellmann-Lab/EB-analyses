library(tidyverse)
library(cowplot)
library(ggsankey)
library(Seurat)
library(Matrix)
library(ggh4x)
library(ggdendro)

# * Colour schemes #####
spec_colors <-        c("#B0144F", "#92AEC8", "#3AA6A6", "#F2A518")
names(spec_colors) <- c("human","rhesus","cyno", "orang")

ind_colors <-        c("#871744","#D81159",   "#E59A10", "#FFBC42", "#218380", "#52B3B7", "#83c5be", "#014f86", "#2c7da0", "#89c2d9")
names(ind_colors) <- c("human_29B5","human_63Ab2.2", "orang_69A1","orang_68A20","cyno_82A3","cyno_56B1", "cyno_56A1","rhesus_87B1", "rhesus_83D1", "rhesus_83Ab1.1")

ct_colors <- setNames(c("#a06cd5", 
                        "#ddf5ff", "#90e0ef", "#00b4d8", "#023e8a",
                        "#68d8d6", "#07beb8",
                        "#ffcad4", "#FF4782","#FF4242", "#cb0b0a","#8e0413", 
                        "#c5edac" ,"#81c14b", "#1a7431",
                        "lightgrey"),
                      nm= c("iPSCs", 
                            "early_ectoderm", "astrocyte_progenitor","granule_precursor_cells","neurons",  
                            "neural_crest_I", "neural_crest_II",
                            "fibroblasts","smooth_muscle_cells", "cardiac_progenitor_cells","cardiac_endothelial_cells","cardiac_fibroblasts",
                            "early_epithelial_cells" ,"epithelial_cells","hepatocytes",
                            "unassigned"))

cons_clust_colors <- setNames(c("#a06cd5",
                                "#ddf5ff", "#90e0ef","#023e8a","#68d8d6", "#07beb8",
                                "#FF4782","#cb0b0a","#8e0413", "coral3",
                                "#81c14b", "#1a7431", 
                                "lightgrey", "grey60", "grey40", "grey20","grey50"),
                              nm= c(4, 
                                    2,12,6,5,13,
                                    3,7,10,1,
                                    17,16,
                                    14,8,9,15,11))
names(cons_clust_colors) <- paste0("m", names(cons_clust_colors))

sankey_colors <- c(spec_colors, ct_colors,cons_clust_colors)

# * Load data ####
metadata_full <- readRDS("zenodo/cell_type_assignment/metadata_full.RDS")
seu_list <- readRDS("zenodo/cell_type_assignment/seu_list.RDS")

#.................................................................................................................. ####
# Figure 2 ####
# 2B) Sankey plot #####
species_order <- c("human","orang","cyno","rhesus")
ct_order <- names(ct_colors)
consensus_merged_order <- names(cons_clust_colors)
consensus_order <- paste0("c", c(6,12,
                                 5,10,15,23, 21,14,13,16, 24,
                                 11, 17, 2,3,4, 1,9,22,
                                 8, 7, 
                                 25, 18, 19, 26, 20))
spec_clust_order <- c(paste0("human__",1:100), paste0("orang__",1:100), paste0("cyno__",1:100), paste0("rhesus__",1:100))
spec_clust_order <- spec_clust_order[spec_clust_order %in% unique(metadata_full$spec_clust)]

sankey_order <- rev(c(species_order, ct_order, consensus_merged_order, consensus_order, spec_clust_order))


df_sankey <- metadata_full %>% 
  filter(!is.na(manual_annotation)) %>% 
  mutate(consensus = paste0("c",consensus),
         consensus_merged = paste0("m",consensus_merged)) %>% 
  select(species, spec_clust, consensus, consensus_merged, manual_annotation) %>% 
  distinct() %>% 
  make_long(species, spec_clust, consensus, consensus_merged,manual_annotation) %>% 
  mutate(fill = case_when(x %in% c("species","manual_annotation", "consensus_merged") ~ node,
                          x == "spec_clust" ~ gsub("__.*","",node),
                          x %in% c("consensus") ~ next_node,
                          T ~ "other"),
         node = factor(node, levels = sankey_order),
         next_node = factor(next_node, levels = sankey_order), 
         label = ifelse(x %in% c("species", "consensus", "consensus_merged", "manual_annotation"), as.character(node), NA)) %>% 
  drop_na(node)

fig2b <- ggplot(df_sankey, aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = fill,
                      label = label
)) +
  geom_sankey(type = "sankey", space = 0.3, node.color = 1)+
    scale_fill_manual(values = sankey_colors)+
    theme_void()+
    theme(legend.position = "none")

fig2b


# 2C) Stacked bar of cell types #####
fig2c <- metadata_full %>% 
  filter(manual_annotation != "unassigned") %>%
  mutate(species = factor(species, levels = c("rhesus", "cyno", "orang", "human")),
         manual_annotation = factor(manual_annotation, levels = rev(names(ct_colors)))) %>%
  ggplot(aes(x = species, fill = manual_annotation))+
  geom_bar()+
  scale_fill_manual(values = ct_colors)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 12),axis.text.x = element_text(size=12), axis.title.x = element_text(size=12)
        )+
  labs(y = "Number of cells")

fig2c

# celltype legend
ct_legend <- get_legend(ggplot(metadata_full, aes(x = species, fill = manual_annotation))+
  geom_bar()+
  scale_fill_manual(values = ct_colors)+
  theme(legend.position = "bottom", legend.title = element_blank())+
  guides(fill = guide_legend(nrow = 4)))


# 2D) UMAP per species #####
fig2d <- metadata_full %>% 
  filter(manual_annotation != "unassigned") %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = manual_annotation))+
  geom_point(size = 0.2)+
  scale_color_manual(values = ct_colors)+
  facet_grid(~species, scales = "free")+
  theme_bw()+
  labs(x = "UMAP1", y = "UMAP2")+
  theme(legend.position = "none",
        strip.text = element_text(size = 12), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

fig2d


# 2E) Heatmap #####

# remove unassigned cells from seurat objects
seu_list_filt <- lapply(seu_list, function(seu){
  tidyseurat::filter(seu, manual_annotation != "unassigned")
})

# List of representative marker genes
manual_markers <- Reduce(rbind, list(
  data.frame(celltype = "iPSCs", gene = c("POU5F1","NANOG","L1TD1")),               
  data.frame(celltype = "early_ectoderm", gene = c("SOX2","HES5","RFX4")),           
  data.frame(celltype = "astrocyte_progenitor", gene = c("SPON1","TFF3","SHH")),    
  data.frame(celltype = "granule_precursor_cells", gene = c("NFIA","ZIC1","ZIC4")), 
  data.frame(celltype = "neural_crest_I", gene = c("SOX10","S100B","FOXD3")),      
  data.frame(celltype = "neural_crest_II", gene = c("FXYD1","NPR3","SERPINA3")),   
  data.frame(celltype = "neurons", gene = c("STMN2","TAGLN3","DCX")),                
  data.frame(celltype = "fibroblasts", gene = c("PDGFRA","PCOLCE","PRRX1")),         
  data.frame(celltype = "smooth_muscle_cells", gene = c("COL8A1","ACTG2","ACTA2")), 
  data.frame(celltype = "cardiac_endothelial_cells", gene = c("NNMT","AQP1","CEBPD")),   
  data.frame(celltype = "cardiac_fibroblasts", gene = c("TNNT2","DCN","HAND2")), 
  data.frame(celltype = "cardiac_progenitor_cells", gene = c("HAS2","MESP1","APLNR")),
  #data.frame(celltype = "early_epithelial_cells", gene = c("EPCAM","")),
  data.frame(celltype = "epithelial_cells", gene = c("CDH1","EPCAM","CLDN7")),   
  data.frame(celltype = "hepatocytes", gene = c("TTR","APOA1","APOA2"))       
))


# Summarize expression fraction per gene & cell type
gene_expr <- lapply(seu_list_filt, function(seu){
  lapply(unique(seu$manual_annotation), function(ctype){
    cells <- seu@meta.data %>%
      rownames_to_column("cell_barcode") %>% 
      filter(manual_annotation == ctype) %>%
      pull(cell_barcode)
    
    frac <- Matrix::rowSums(seu@assays$RNA@counts[manual_markers$gene, cells] > 0) / ncol(seu@assays$RNA@counts[manual_markers$gene, cells])
    
    return(data.frame(manual_annotation = ctype,
                      gene = names(frac),
                      pct.expr = frac))
  }) %>% bind_rows()
}) %>% bind_rows(.id = "species")
    
# Filter for cell types present in at least 3 species: 
gene_expr_3plus <- gene_expr %>% 
  group_by(manual_annotation) %>% 
  filter(length(unique(species)) > 2) %>% 
  ungroup() %>% 
  left_join(manual_markers) %>% 
  filter(celltype %in% unique(manual_annotation)) %>%
  mutate(manual_annotation = factor(manual_annotation, levels = ct_order),
         celltype = factor(celltype, levels = ct_order),
         species = factor(species, levels = rev(c("human","orang","cyno","rhesus"))))


# ** Plot heatmap
fig2e <- plot_grid(
  ggplot(gene_expr_3plus, aes(x = "", y = species, fill = species))+
    geom_tile(color = "lightgrey")+
    facet_grid2(manual_annotation~., scales = "free", space = "free")+
    scale_fill_manual(values = spec_colors)+
    #scale_y_discrete(position = "right")+
    theme_minimal()+
    theme(strip.background.y = element_blank(),strip.text.y = element_blank(),
          axis.title = element_blank(),
          legend.position = "none"),
  
  
  ggplot(gene_expr_3plus, aes(x = gene, y = species, fill = pct.expr))+
    geom_tile(color = "lightgrey")+
    facet_grid2(manual_annotation~celltype, scales = "free", space = "free", 
                strip = strip_themed(background_y = elem_list_rect(fill = c("#dabfff", rep("lightblue",4), rep("lightpink1",2), rep("#c5edac",2)))))+
    scale_fill_gradient(low = "grey100", high = "darkblue")+
    scale_x_discrete(position = "top")+
    theme_bw()+
    theme(axis.text.x.top =  element_text(angle = 45, hjust = 0, vjust = 0, size = 8),
          #axis.text.y = element_text(size = 7),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
          strip.background.x = element_blank(),strip.text.x = element_blank(),
          axis.title = element_blank()),
  
  nrow = 1, rel_widths = c(0.1,1), align = "h", axis = "tb"
)
fig2e


# Assemble figure ####
fig2 <- plot_grid(
  plot_grid(
    plot_grid(
      NULL, # workflow scheme goes here
      fig2b,
      nrow = 2, rel_heights = c(0.3,1), labels = c("A","B"), scale = 0.95
    ),
    fig2c,
    nrow = 1, rel_widths = c(1,0.7), labels = c("","C"), scale = c(1,1,0.95)
  ),
  NULL,
  ct_legend,
  NULL,
  fig2d,
  fig2e,
  
  ncol = 1, rel_heights = c(0.8,0.1,0.25,0.05,0.6,1.3), labels = c("","","","","D","E"), scale = c(1,1,0.9,1,1,0.98)
)

ggsave("Figures/Figure2.pdf", plot = fig2, width = 8.27, height = 11.7)

#.................................................................................................................. ####
# Supplementary Figure S4 and S5 ####
RC_results <- readRDS("zenodo/cell_type_assignment/correspondence_validation/reciprocal_classification_out.RDS")
MetaNeighbor_results <- readRDS("zenodo/cell_type_assignment/correspondence_validation/MetaNeighbor_out.RDS")

# ** plotting functions ####
plot_distmap_and_hc <- function(mat_sym, cluster_colors, species_colors = species_colors, mid.point = 0.4, legend.title = ""){
  
  # hierarchical clustering
  d <- as.dist(1-mat_sym)
  hc <- hclust(d, method = "average")
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc, type = "rectangle")
  
  label_data <- ddata$labels %>% 
    mutate(species = word(label,1,1,"__"),
           celltype = word(label,2,2,"__")) 
  
  plot_annotation <- function(label_data, color_by){
    ggplot(label_data, aes_string(x="x", y="1", fill = paste0("as.factor(",color_by,")")))+
      geom_bar(stat = "identity")+
      theme_nothing()+
      theme(legend.position = "none", axis.text = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
            legend.text = element_text(size = 10))+
      labs(y = color_by)
  }
  
  p.dendro <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    theme_nothing()
  
  p.anno_ct <-  plot_annotation(label_data, "celltype")+scale_fill_manual(values = cluster_colors, breaks = names(cluster_colors))
  p.anno_species <- plot_annotation(label_data, "species")+scale_fill_manual(values = species_colors, breaks = names(species_colors))
  
  p.mat <- reshape2::melt(as.matrix(mat_sym)) %>% 
    mutate(Var1 = factor(Var1, levels = rownames(as.matrix(mat_sym))[hc$order]),
           Var2 = factor(Var2, levels = rownames(as.matrix(mat_sym))[hc$order])
    ) %>% 
    ggplot(aes(x = as.numeric(Var1), y = Var2, fill = value))+
    geom_tile(color = "black")+
    scale_fill_gradient2(low = "lightblue", mid = "white", high = "red", midpoint = mid.point)+
    labs(fill = legend.title)+
    theme_bw()+
    theme(axis.title = element_blank(), #legend.title = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none")
  
  ct_legend <- get_legend(p.anno_ct+theme(legend.position = "left",legend.justification = "left", legend.title = element_blank(), legend.key.size = unit(0.4, 'cm')))
  species_legend <- get_legend(p.anno_species+theme(legend.position = "left",legend.justification = "left", legend.title = element_blank(), legend.key.size = unit(0.4, 'cm')))
  heatmap_legend <- get_legend(p.mat+theme(legend.position = "left", legend.justification = "left"))
  
  p <- plot_grid(
    plot_grid(p.dendro,
              p.anno_ct,      
              p.anno_species,
              p.mat,
              ncol = 1, align = "v", axis = "lr", rel_heights = c(0.3,0.1,0.1,1), scale = c(0.985,1,1,1)),
    plot_grid(ct_legend, species_legend,NULL,heatmap_legend,NULL, ncol = 1, rel_heights = c(1,0.3,0.1,0.5,0.1)),
    nrow = 1, rel_widths = c(1,0.3)
  )
  
  return(p)
}
ct_dictionary = setNames(c("iPSCs",
                           "EE","NCI","NCII","AP","GPC","Neu",
                           "CEndo","SMC","Fib","CFib","MDnR","CPC",
                           "EEC","EC","Hepa",
                           "unass"),
                         c("iPSCs", 
                           "early_ectoderm","neural_crest_I","neural_crest_II","astrocyte_progenitor" ,"granule_precursor_cells","neurons",
                           "cardiac_endothelial_cells", "smooth_muscle_cells","fibroblasts",  "cardiac_fibroblasts", "mesoderm_noRibo", "cardiac_progenitor_cells",
                           "early_epithelial_cells", "epithelial_cells",   "hepatocytes",
                           "unassigned"))

create_pairwise_plot <- function(data, species_x, species_y) {
  sim_filt <- data %>% 
    data.frame() %>% 
    rownames_to_column("cluster_x") %>% 
    pivot_longer(-cluster_x, names_to = "cluster_y", values_to = "similarity") %>%
    
    mutate(
      c1 = pmin(cluster_x, cluster_y),  # Always the alphabetically first
      c2 = pmax(cluster_x, cluster_y)   # Always the alphabetically second
    ) %>%
    distinct(c1, c2, .keep_all = TRUE) %>%  # Keep unique pairs
    select(-c1, -c2) %>% 
    
    mutate(species.x = word(cluster_x,1,1, sep= "__"),
           species.y = word(cluster_y,1,1, sep= "__"),
           cluster.x = word(cluster_x,2,2, sep= "__"),
           cluster.y = word(cluster_y,2,2, sep= "__")) %>% 
    group_by(species.x, species.y) %>% 
    filter(cluster.y %in% cluster.x, 
           cluster.x %in% cluster.y,
           species.x == species_x & species.y == species_y) %>% 
    mutate(species.x = factor(species.x, levels = c("human","orang","cyno","rhesus")),
           species.y = factor(species.y, levels = rev(c("human","orang","cyno","rhesus"))),
           ct_short.x = ct_dictionary[cluster.x],
           ct_short.y = ct_dictionary[cluster.y],
           ct_short.x = factor(ct_short.x, levels = ct_dictionary),
           ct_short.y = factor(ct_short.y, levels = ct_dictionary))
  
  ggplot(sim_filt, aes(x = ct_short.x, y = ct_short.y, fill = similarity)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "lightblue", mid = "white", high = "red", midpoint = 0.4) +
    facet_grid(species.y ~ species.x, scales = "free") +
    theme_bw() +
    labs(x = "query", y = "reference", fill = "classification\nfraction") +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),
      axis.title = element_blank(),
      strip.text = element_text(size = 13),
      legend.position = "none")
}
plot_all_combinations <- function(sim){
  # Generate individual plots for each pairwise species combination
  plot1 <- create_pairwise_plot(sim, "human", "orang")
  plot2 <- create_pairwise_plot(sim, "human", "cyno")
  plot3 <- create_pairwise_plot(sim, "human", "rhesus")
  plot4 <- create_pairwise_plot(sim, "cyno", "orang")
  plot5 <- create_pairwise_plot(sim, "rhesus", "orang")
  plot6 <- create_pairwise_plot(sim, "cyno", "rhesus")
  
  # Combine the plots into a 3x2 matrix using patchwork
  p <- (plot1 + plot2 + plot3) / (plot4 + plot5 + plot6)
  
  return(p)
}

# ** Assemble S4 ####
figS4a <- plot_distmap_and_hc(RC_results, cluster_colors = ct_colors, species_colors = spec_colors, mid.point = 0.4, legend.title = "classification\nfraction")
figS4b <- plot_all_combinations(RC_results)
plot_grid(figS4a, figS4b,
          nrow = 2, rel_heights = c(1,0.9), labels = c("A","B"))

ggsave("Figures/suppl_figure_4_validation_RC.pdf", width = 8.3, height = 11.2)

# ** Assemble S5 ####
figS5a <- plot_distmap_and_hc(MetaNeighbor_results, cluster_colors = ct_colors, species_colors = spec_colors, mid.point = 0.4, legend.title = "classification\nfraction")
figS5b <- plot_all_combinations(MetaNeighbor_results)
plot_grid(figS5a, figS5b,
          nrow = 2, rel_heights = c(1,0.9), labels = c("A","B"))

ggsave("Figures/suppl_figure_5_validation_MetaNeighbor.pdf", width = 8.3, height = 11.2)




#.................................................................................................................. ####
# Supplementary Figure S6 ####
metadata_full2 <- metadata_full %>% 
  filter(manual_annotation != "unassigned") %>%
  mutate(manual_annotation = factor(manual_annotation, levels = rev(names(ct_colors))))

# ** S6A ####
figS6a <- ggplot(metadata_full2, aes(x = individual, fill = manual_annotation))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = ct_colors)+
  coord_flip()+
  facet_grid(species~., scales = "free", space = "free")+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 11),axis.text.x = element_text(size=10), axis.title.x = element_text(size=12),
        strip.text = element_text(size = 12)
  )+
  labs(y = "Fraction of cells")

# ** S6B ####
figS6b <- ggplot(metadata_full2, aes(x = experiment, fill = manual_annotation))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = ct_colors)+
  coord_flip()+
  facet_grid(orig.ident~., scales = "free", space = "free")+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 11),axis.text.x = element_text(size=10), axis.title.x = element_text(size=12),
        strip.text = element_text(size = 12)
  )+
  labs(y = "Fraction of cells")

# ** S6C ####
figS6c <- ggplot(metadata_full2, aes(x = day, fill = manual_annotation))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = ct_colors)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 11),axis.text.x = element_text(size=10), axis.title.x = element_text(size=12)
  )+
  labs(y = "Fraction of cells")

ct_legend_S6 <- get_legend(figS6a+theme(legend.position = "bottom", legend.title = element_blank())+scale_fill_manual(values = ct_colors, breaks = names(ct_colors)))

# ** Assemble figure ####
figS6 <- plot_grid(
  plot_grid(figS6a,figS6b, labels = c("A","B"), scale = 0.98),
  figS6c,
  ct_legend_S6,
  nrow = 3, rel_heights = c(1,0.4,0.2), labels = c("","C"), scale = c(1,0.98,1)
)

ggsave("Figures/suppl_figure_6_celltype_proportions.pdf",figS6, width = 8.3, height = 6)

