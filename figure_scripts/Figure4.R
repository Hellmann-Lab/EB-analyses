library(tidyverse)
library(cowplot)
library(Seurat)

marker_overlap_full <- readRDS("zenodo/marker_gene_analysis/marker_overlap.rds")
rbo_biotype <- readRDS("zenodo/marker_gene_analysis/rbo_biotype.rds")
F1_macro_bootstrap <- readRDS("zenodo/marker_gene_analysis/F1_macro_bootstrap.rds")

#.................................................................................................................. ####
# Figure 4 ####
# ** 4A: UpsetR plot of marker overlap ####
fig4a <- marker_overlap_full %>% 
  count(celltype,species_list, n_species) %>% 
  
  ggplot(aes(x = species_list, y = n, color = n_species))+
  geom_jitter(width = 0.2, alpha = 0.8)+
  stat_summary(fun = "median", fun.min = "median", fun.max= "median", size= 0.3, geom = "crossbar")+
  scale_x_upset(order_by = "degree", reverse = T, sets=c("human", "orang", "cyno", "rhesus"))+
  scale_color_manual(values = c("#295277","#018571", "#8BA2B2", "#8DD1C6" ))+
  labs(y = "Number of shared markers")+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank())

# ** 4B: Heatmap of conserved and human-specific markers ####

# get expression fractions from filtered sce object
sce_filt <- readRDS("zenodo/marker_gene_analysis/sce_filt.rds")
cnts <- counts(sce_filt)
metadata <- data.frame(colData(sce_filt)) %>% 
  select(species, individual, manual_annotation) %>% 
  rownames_to_column("cell_barcode")

# Calculate gene counts for a given individual and celltype
get_gene_counts <- function(spec, ctype) {
  # Filter the cell barcodes for this individual and celltype
  cells <- metadata %>%
    filter(species == spec, manual_annotation == ctype) %>%
    pull(cell_barcode)
  
  # Subset the count matrix for these cells
  gene_counts <- rowSums(cnts[, cells, drop = FALSE])
  
  # Return a dataframe with gene counts
  tibble(
    gene = rownames(cnts),
    gene_count = gene_counts, 
    ncells_expr = rowSums(cnts[, cells, drop = FALSE] > 0),
    ncells = length(cells)
  )
}

# Apply the function to all combinations of individual and celltype using map()
gene_counts_summary <- metadata %>%
  distinct(species, manual_annotation) %>%
  mutate(data = map2(species, manual_annotation, get_gene_counts)) %>%
  unnest(data)

# Plot heatmaps while highlighting example genes
df_heat <- marker_overlap_full %>% 
  left_join(gene_counts_summary, relationship = "many-to-many") %>% 
  mutate(expr_frac = ncells_expr / ncells, 
         species = factor(species, levels = rev(c("human","orang","cyno","rhesus"))))

# add abbreviations
ct_dictionary = setNames(c("iPSCs","EE","NCI","SMC","CFib","EC","Hepa"),
                         c("iPSCs","early_ectoderm","neural_crest_I","smooth_muscle_cells","cardiac_fibroblasts", "epithelial_cells","hepatocytes"))

df_heat <- mutate(df_heat, 
             celltype_short_expr = ct_dictionary[manual_annotation],
             celltype_short_expr = factor(celltype_short_expr, levels = ct_dictionary))


# shared marker genes in top 100: 
p.cons_markers <- df_heat %>% 
  filter(n_species == 4) %>% 
  ggplot(aes(x =  reorder(gene,as.numeric(as.factor(celltype_short))), y = species, fill = expr_frac))+
  geom_tile()+
  facet_grid(celltype_short_expr ~ ., scales = "free", space = "free")+
  scale_fill_gradient(low = "white", high = "#018571")+
  scale_x_discrete(breaks = c("POU5F1", "HES5", "SOX10", "COL8A1", "DCN", "APOA1", "CDH1"))+
  labs(fill = "Expression\nfraction")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.title = element_blank(), legend.position = "bottom", 
        strip.background = element_blank(), strip.text = element_blank())


# human specific marker genes in top 100: 
p.human_markers <- df_heat %>% 
  filter(n_species == 1 & species_list == "human") %>% 
  ggplot(aes(x = reorder(gene,as.numeric(as.factor(celltype_short))), y = species, fill = expr_frac))+
  geom_tile()+
  facet_grid(celltype_short_expr ~ ., scales = "free", space = "free")+
  scale_fill_gradient(low = "white", high = "#295277")+
  scale_x_discrete(breaks = c("SCGB3A2","SDK2", "BCHE", "NT5E", "ATP7B", "ITLN2", "PAPPA2"))+
  labs(fill = "Expression\nfraction")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(), legend.position = "bottom")

fig4b <- plot_grid(p.cons_markers, p.human_markers, nrow = 1, rel_widths = c(1,1.5), align = "h", axis = "tb")


# ** 4C: RBO by biotype ####
fig4c <- rbo_biotype %>% 
  filter(
    nMarkers == 100, 
    pair %in% c("human_orang","human_cyno", "human_rhesus", "cyno_rhesus")   # only take human comparisons and within clade comparison
  ) %>% 
  
  ggplot(aes(x = reorder(pair, rbo), y = rbo, fill = reorder(biotype, rbo)))+
  geom_boxplot()+
  coord_flip()+
  scale_fill_manual(values = c("#E4DE70", "#468F85", "#E4A070","#76519A" ))+
  labs(y = "RBO per cell type")+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        legend.title = element_blank(), 
        axis.text = element_text(size = 10), 
        legend.position = "top",
        legend.background = element_rect(fill='transparent'))

# ** 4D: F1 of kNN classification ####
spec_colors <-        c("#B0144F", "#92AEC8", "#3AA6A6", "#F2A518")
names(spec_colors) <- c("human","rhesus","cyno", "orang")

fig4d <- ggplot(data = F1_macro_bootstrap, aes(x = as.numeric(NoMarkers))) + 
  geom_line(aes(x = as.numeric(NoMarkers), y = F1_macro, color = TestSpecies, group = Set), linewidth = 1) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high, fill = TestSpecies, group = Set), alpha=0.1)+
  scale_color_manual(values = spec_colors)+
  scale_fill_manual(values = spec_colors)+
  scale_y_continuous(limits = c(0,1))+
  facet_wrap(.~TestSpecies, nrow = 1) +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(#title = "KNN Classification across Species using human markers \ntrained in human clone 29B5 and tested in other clones", 
    x = " Number of Marker Genes",
    y = "macro-averaged F1-Score")

# ** Assemble figure ####
plot_grid(
  plot_grid(fig4a, 
            fig4c,
            nrow = 2, labels = c("A","C"), rel_heights = c(1,0.81), scale = c(0.96, 0.96)),
  plot_grid(fig4b,
            fig4d,
            nrow = 2, labels = c("B","D"), rel_heights = c(1.3,0.8), scale = c(0.96, 0.96)),
  nrow = 1, rel_widths = c(1,1.6)
)




#.................................................................................................................. ####
# Supplementary Figure S8 ####
seu_list <- readRDS("zenodo/cell_type_assignment/seu_list.RDS")
metadata_full <- readRDS("zenodo/cell_type_assignment/metadata_full.RDS")

expr_values <- lapply(seu_list, function(seu){
  FetchData(object = seu, vars = marker_overlap_full$gene) %>% rownames_to_column("BC")
}) %>% bind_rows()

md <- metadata_full %>% 
  rownames_to_column("BC") %>% 
  filter(manual_annotation %in% unique(marker_overlap_full$celltype)) %>% 
  left_join(expr_values)

# custom theme: only frame
theme_frame <- theme_minimal() +
  theme(
    axis.title = element_blank(),          # Remove axis titles
    axis.text = element_blank(),           # Remove axis text (tick labels)
    axis.ticks = element_blank(),          # Remove axis ticks
    panel.grid = element_blank(),          # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Add frame border
    axis.line = element_blank() ,           # Remove the default axis lines
    panel.spacing = unit(0.1, "lines")
  )

plot_UMAP_gene <- function(gene, color_high = "darkblue"){
  ggplot(md, aes_string(x = "UMAP_1", y = "UMAP_2", color = gene))+
      geom_point(size = 0.2)+
      scale_color_gradient(low = "lightgrey", high = color_high)+
      facet_wrap(~species, nrow = 1, scales = "free")+
      theme_frame+
      theme(strip.background = element_blank(),strip.text.x = element_blank(),
        legend.position = "right")
}

# ** S8A ####
# Colored by cell type as reference: 
figS8a <- ggplot(md, aes_string(x = "UMAP_1", y = "UMAP_2", color = "manual_annotation"))+
  geom_point(size = 0.5)+
  scale_color_manual(values = ct_colors)+
  facet_wrap(~species, nrow = 1, scales = "free")+
  labs(color = "cell type")+
  theme_frame+
  theme(strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "right")+
  guides(colour = guide_legend(override.aes = list(size=3)))

# ** S8B ####
# examples of shared markers in all 4 species:
figS8b <- plot_grid(
  plot_UMAP_gene("POU5F1", "#a06cd5"), # iPSC
  plot_UMAP_gene("HES5", "#ddf5ff"), # early ectoderm
  plot_UMAP_gene("SOX10", "#68d8d6"), # neural crest
  plot_UMAP_gene("COL8A1", "#FF4782"), # smooth muscle cells
  plot_UMAP_gene("DCN", "#8e0413"), # cardiac fibroblasts
  plot_UMAP_gene("APOA1", "#1a7431"), # hepatocytes
  plot_UMAP_gene("CDH1", "#81c14b"), # epithelial
  ncol = 1, align = "v", axis = "lr"
)

# ** S8C ####
# Some human specific markers:
figS8c <- plot_grid(
  plot_UMAP_gene("SCGB3A2", "#a06cd5"), # iPSC
  plot_UMAP_gene("SDK2", "#ddf5ff"), # early ectoderm
  plot_UMAP_gene("BCHE", "#68d8d6"), # neural crest
  plot_UMAP_gene("NT5E", "#FF4782"), # smooth muscle cells
  plot_UMAP_gene("ATP7B", "#8e0413"), # cardiac fibroblasts
  plot_UMAP_gene("ITLN2", "#1a7431"), # hepatocytes
  plot_UMAP_gene("PAPPA2", "#81c14b"), # epithelial
  ncol = 1, align = "v", axis = "lr"
)

# Assemble figure
plot_grid(
  figS8a,
  plot_grid(figS8b, figS8c, nrow = 1, labels = c("B","C")), 
  nrow = 2, rel_heights = c(0.3,1), labels = c("A",""))

#.................................................................................................................. ####
# Supplementary Figure S9 ####
crossspecies.perform.ct.data.filt <- readRDS("zenodo/marker_gene_analysis/F1_per_clone.rds")

ggplot(crossspecies.perform.ct.data.filt) + 
  geom_line(aes(x = as.numeric(NoMarkers), y = F1, color = TestSpecies, group = Set), linewidth = 1) +
  facet_wrap(~Celltype) +
  scale_color_manual(values = spec_colors)+
  theme_bw() +
  labs(
    x = "Number of Marker Genes",
    y = "F1-Score")+
  theme(legend.position = c(0.4,0.16),
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.background = element_rect(fill='transparent'))
