library(tidyverse)
library(cowplot)
library(Seurat) 


#.................................................................................................................. ####
# Figure 1 ####

# read in Harmony corrected Seurat object
seu <- readRDS("zenodo/processing_per_species/seu_integrated_allSpecies.RDS")

# Prepare Seurat Object
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
seu <- FindClusters(seu, resolution = 0.1)

# assign germlayer based on classification
seu$germlayer <- case_when(seu$labels %in% c("Early_Ectoderm", "Neural_Crest", "Neurons") ~ "ectoderm",
                           seu$labels %in% c("Endoderm") ~ "endoderm", 
                           seu$labels %in% c("Mesoderm", "Endothelial_Cells") ~ "mesoderm",
                           seu$labels %in% c("Pluripotent_Cells") ~ "pluripotent")

# ** Figure 1D: Marker overlay plot #################################################################################################
genes <- c("SOX2","POU5F1","COL1A1", "ACTA2","EPCAM", "APOA1", "NANOG","SOX10","STMN4")
gene_colors <- setNames(c("deepskyblue2", "cyan4", "#03045e", "#a06cd5", "#c19ee0", "#1a7431","#81c14b","#E00043","darkred"), 
                        c("SOX2","SOX10","STMN4","POU5F1", "NANOG", "APOA1","EPCAM", "COL1A1", "ACTA2"))

# extract normalized count data for candidate genes:
df <- seu@reductions[["umap.harmony"]]@cell.embeddings %>% 
  cbind(FetchData(seu, c(genes, "species"), slot = "data"))
colnames(df)[1:2] <- c("x","y")

# plot gene expression overlay
p <- ggplot(df, aes(x = x, y = y))+
  geom_point(color = "lightgrey", size = 0.1, alpha = 0.2)+
  theme_bw()

#genes_order <- c(3,1,8,7,6,2,4,5)
genes_order <- c(3,1,8,9,7,6,2,4,5)
for(i in 1:9){
  gene <- genes[[i]]
  order <- genes_order[[i]]
  p <- p+
    ggnewscale::new_scale_colour()+
    geom_point(data = df[df[[gene]] > 0,], aes_string(x = "x", y = "y", color = gene), size = 0.5, alpha = 0.4)+
    scale_color_gradient(low = "lightgrey", high = gene_colors[[gene]],
                         guide = guide_colorbar(order = order,  title.position = "right"))
}

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2, "cm")
)

fig1d <- p+theme_void()+
  labs(x = "UMAP1", y = "UMAP2")+
  guides(x = axis, y = axis) +
  theme(axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.1, "inches"))),
        axis.title.x = element_text(hjust = 0.02),
        axis.title.y = element_text(hjust = 0.02, angle = 90),
        legend.direction = "horizontal")

fig1d


# ** Figure 1E: germ layer assignment #######################################################################
germlayer_colors <- setNames(c("#BF0D0D", "#1982C4", "#81c14b" ,"#dabfff"),
                             nm= c("mesoderm","ectoderm", "endoderm", "pluripotent"))

df.rhodes <- seu@reductions[["umap.harmony"]]@cell.embeddings %>% 
  cbind(FetchData(seu, c("species", "labels", "germlayer", "seurat_clusters"))) %>% 
  dplyr::rename(UMAP1 = umapharmony_1, UMAP2 = umapharmony_2) %>% 
  mutate(species = factor(species, levels = c("human","orang", "cynomolgus", "rhesus")))

fig1e <- ggplot(df.rhodes, aes(x = UMAP1, y = UMAP2, color = germlayer))+
  geom_point(size = 0.1, alpha = 0.7)+
  facet_wrap(~species, scales = "free", nrow = 2)+
  scale_color_manual(values = germlayer_colors)+
  theme_void()+
  theme(legend.position = c(0.57,0.52), legend.title = element_blank(), 
        strip.text = element_text(size = 14), legend.text = element_text(size = 12))+
  guides(colour = guide_legend(override.aes = list(size=3)))

fig1e

#.................................................................................................................. ####
# Supplementary Figure S3 ####

# ** S3A ####
# Rhodes cell type assignments
rhodes_colors <- setNames(c("#c49be8", 
                            "#90e0ef", "#407ba7", "#023e8a", 
                            "#ffafcc","#c9184a",
                            "#1a7431"),
                          nm= c("Pluripotent_Cells",  
                                "Early_Ectoderm",  "Neural_Crest", "Neurons",
                                "Mesoderm", "Endothelial_Cells",
                                "Endoderm"))

figS3a <- ggplot(df.rhodes, aes(x = UMAP1, y = UMAP2, color = labels))+
  geom_point(size = 0.1, alpha = 0.7)+
  facet_wrap(~species, scales = "free", nrow = 1)+
  scale_color_manual(values = rhodes_colors)+
  theme_void()+
  theme(legend.position = "right", legend.title = element_blank(), 
        strip.text = element_text(size = 12), legend.text = element_text(size = 11))+
  guides(colour = guide_legend(override.aes = list(size=3)))


# ** S3B #####
# Integrated clusters

figS3b <- ggplot(df.rhodes, aes(x = UMAP1, y = UMAP2, color = seurat_clusters))+
  geom_point(size = 0.1, alpha = 0.7)+
  theme_void()+
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 11))+
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))


# ** S3C ####
# clusters vs Rhodes cell type assignment
figS3c <- ggplot(df.rhodes, aes(x = seurat_clusters, fill = labels))+
  geom_bar(position = "fill")+
  facet_wrap(~species,nrow = 1)+
  scale_fill_manual(values = rhodes_colors)+
  labs(x = "Clusters in integrated data", y = "Fraction of cells")+
  theme_bw()+
  theme(legend.position = "none", strip.text = element_text(size = 12))

# ** Assemble Suppl. Fig. ####
plot_grid(figS3a, 
          plot_grid(
            figS3b, 
            NULL,
            figS3c, 
            nrow = 1, rel_widths = c(1,0.1,3), labels = c("B","C")
          ),
          nrow = 2, labels = c("A",""))

