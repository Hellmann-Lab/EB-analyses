library(tidyverse)
library(cowplot)
library(Seurat)

pleio_phast <- readRDS("zenodo/cell_type_specificity/ct_specificity_summary.RDS")
gene_detection_allSpec <- readRDS("zenodo/cell_type_specificity/gene_detection_allSpec.RDS")
seu_list <- readRDS("zenodo/cell_type_assignment/seu_list.RDS")
metadata_full <- readRDS("zenodo/cell_type_assignment/metadata_full.RDS")

#.................................................................................................................. ####
# Figure 3 ####
# ** 3A: UMAPs with expression of example genes ####
metadata <- metadata_full %>% 
  filter(manual_annotation != "unassigned") %>% 
  rownames_to_column("BC")

plot_UMAP_gene <- function(gene){
  expr_values <- lapply(seu_list, function(seu){
    FetchData(object = seu, vars = gene) %>% rownames_to_column("BC")
  }) %>% bind_rows()
  
  md <- left_join(metadata, expr_values)
  
  p <- ggplot(md, aes_string(x = "UMAP_1", y = "UMAP_2", color = gene))+
    geom_point(size = 0.3)+
    #geom_point(data = filter(md, !!sym(gene) > 0), aes_string(color = gene), size = 0.5)+
    scale_color_gradient(low = "lightgrey", high = "darkblue")+
    facet_wrap(~species, nrow = 2, scales = "free")+
    labs(color = "expr")+
    ggtitle(gene)+
    theme_void()+
    theme(strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15), legend.position = "bottom")
  
  return(p)
}

fig3A <- plot_grid(
  plot_UMAP_gene("SOX10"),
  plot_UMAP_gene("ESRG"),
  plot_UMAP_gene("RPL22"),
  nrow = 1
)

fig3A

# ** 3B: Expression fractions of example genes ####
plot_gene_specificity_comb2 <- function(gene_sel){
  df <- gene_detection_allSpec %>% filter(gene == gene_sel)
  p <- ggplot(df, aes(x = celltype_short, y = species, fill = detection))+
    geom_tile(color = "darkgrey")+
    geom_tile(data = filter(df, expressed == 1), color = "black", size = 1.5)+
    scale_fill_gradient(low = "grey100", high = "darkblue", na.value = "grey", limits = c(0,1))+
    scale_x_discrete(position = "top")+
    # ggtitle(gene_sel)+
    theme_bw()+
    theme(axis.text.x.top =  element_text(angle = 45, hjust = 0, vjust = 0, size = 12, colour = c("purple","#0096c7","#0096c7","#BF0D0D","#BF0D0D","#1a7431","#1a7431")),
          axis.text.y = element_text(size = 12),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 15), 
          legend.position = "none",
          strip.text.y = element_text(angle = 0, size = 11))
  return(p)
}

combined_legend_frame <- get_legend(
  ggplot(data.frame(a = 1:3, b = 1:3, binarized = c("expressed", "not expressed", "expressed"), d = c(0,0,1)), aes(x=a,y=b,color=binarized, size = binarized))+
    geom_tile(aes(fill = d), color = "darkgrey")+
    geom_tile(fill = "white")+
    scale_color_manual(values = c("black","darkgrey"))+
    scale_size_manual(values = c(1,0.2))+
    scale_fill_gradient(low = "grey100", high = "darkblue", na.value = "grey", limits = c(0,1))+guides(fill=guide_legend(title="expression fraction"))+
    theme_bw()+theme(legend.position = "right")
)

fig3B <- plot_grid(
  plot_gene_specificity_comb2("SOX10"), 
  plot_gene_specificity_comb2("ESRG"),  
  plot_gene_specificity_comb2("RPL22"),
  combined_legend_frame, 
  nrow = 1, rel_widths = c(1,1,1,0.5)
)

fig3B


# ** 3C: Cell type specificity vs expression conservation ####
fig3C <- pleio_phast %>% 
  filter(species == "human") %>% 
  ggplot(aes(x = as.factor(PD), y = conservation_phylo, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2, notch = T)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  #facet_grid(~species)+
  theme_bw()+
  labs(x = "Cell type specificity (human)", 
       #x = "Pleiotropic degree (human)", 
       y = "Expression conservation")+
  theme(legend.position = "none", axis.text = element_text(size=17), axis.title = element_text(size=19))

fig3C

# ** 3D: cell type specificity vs PhastCons scores (43 primates) ####
fig3D <- pleio_phast %>% 
  filter(species == "human") %>% 
  na.omit() %>% 
  ggplot(aes(x = as.factor(PD), y = fracConsPr, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2, notch = T)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  #facet_grid(~species)+
  theme_bw()+
  labs(x = "Cell type specificity (human)", 
       y = "Fraction of constrained CDS bases per gene")+
  theme(legend.position = "none", axis.text = element_text(size=17), axis.title = element_text(size=19))

fig3D


#.................................................................................................................. ####
# Supplementary Figure S7 ####

# ** S7A ####
figS7A <- ggplot(pleio_phast, aes(x = species, fill = factor(PD)))+
  geom_bar(position = "stack")+
  scale_fill_manual(values = colorspace::lighten(scales::colour_ramp(c("lightblue", "darkblue"))(seq(0, 1, len = 7)), amount = 0.35))+
  labs(fill = "Cell type specificity", y = "Number of genes")+
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        legend.position = "top", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

# ** S7B ####
figS7B <-  pleio_phast %>% 
  filter(species != "human") %>%
  ggplot(aes(x = as.factor(PD), y = conservation_phylo, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2, notch = T)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  facet_grid(~species)+
  theme_bw()+
  labs(x = "Cell type specificity", y = "Expression conservation")+
  theme(legend.position = "none")

# ** S7C ####
figS7C <- pleio_phast %>% 
  filter(species != "human") %>% 
  ggplot(aes(x = as.factor(PD), y = fracConsPr, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2, notch = T)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  facet_grid(~species)+
  theme_bw()+
  labs(x = "Cell type specificity", y = "Fraction of constrained CDS bases")+
  theme(legend.position = "none")

# ** S7D ####
figS7D <- ggplot(pleio_phast, aes(x = as.factor(PD), y = mean_expression, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  scale_y_log10()+
  theme_bw()+
  facet_grid(~species)+
  labs(x = "Pleiotropic degree",y = "Mean expression")+
  theme(axis.title.x = element_blank(), 
        legend.position = "none")

# ** S7E ####
# Selecting genes with similar expression level distribution
per_gene_summary_sample <- pleio_phast %>%
  ungroup() %>% 
  mutate(species = factor(species, levels = c("human","orang","cyno","rhesus"))) %>% 
  mutate(bin = cut(mean_expression, breaks = quantile(mean_expression, probs = seq(0, 1, 0.1)),include.lowest = TRUE)) %>% 
  group_by(species, as.factor(PD), bin) %>%
  mutate(n_genes = n()) %>% 
  group_by(bin) %>% 
  mutate(min_n = min(n_genes)) %>%
  group_by(species, as.factor(PD), bin) %>%
  sample_n(unique(min_n)) 


figS7E <- ggplot(per_gene_summary_sample, aes(x = as.factor(PD), y = mean_expression, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  scale_y_log10()+
  theme_bw()+
  facet_grid(~species)+
  labs(x = "Pleiotropic degree",y = "Mean expression\n(subsampled geneset)")+
  theme(axis.title.x = element_blank(), 
        legend.position = "none")   

# ** S7F ####
figS7F <- ggplot(per_gene_summary_sample, aes(x = as.factor(PD), y = conservation_phylo, fill = PD))+
  geom_boxplot(alpha = 0.7, outlier.size = 0.2)+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  facet_grid(~species)+
  theme_bw()+
  labs(x = "Cell type specificity", y = "Expression conservation\n(subsampled geneset)")+
  theme(legend.position = "none")

# ** Assemble figure ####
plot_grid(
  plot_grid(figS7A, 
            plot_grid(figS7B, figS7C, nrow = 2, labels = c("B","C"), scale = 0.9), 
            #align = "h", axis = "b", 
            nrow = 1, labels = c("A",""), scale = c(1,1), rel_widths = c(1.1,2)), 
  plot_grid(figS7D, figS7E, figS7F, 
            ncol = 1, align = "v", axis = "lr", labels = c("D","E","F"), rel_heights = c(1,1,1.15)),
  
  ncol = 1, rel_heights = c(1.4,2)
)

ggsave("/data/share/htp/EBgrant/analysis_scRNA_allRuns/publication_scripts/Figures/suppl_figure_7_celltype_specificity.pdf",width = 11, height = 13)
