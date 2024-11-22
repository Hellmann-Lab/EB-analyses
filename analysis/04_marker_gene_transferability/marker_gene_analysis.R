library(tidyverse)
library(cowplot)
library(patchwork)
library(ggupset)
library(scales)
library(rbo)

# ............................................................................................................... ####
# 1) Marker gene overlap #####

# ** Load top marker lists per species ####
marker_list <- list(
  human = readRDS("zenodo/marker_gene_analysis/MarkerSelection/human_all.rds")$Grand,
  orang = readRDS("zenodo/marker_gene_analysis/MarkerSelection/orang_all.rds")$Grand,
  cyno = readRDS("zenodo/marker_gene_analysis/MarkerSelection/cyno_all.rds")$Grand,
  rhesus = readRDS("zenodo/marker_gene_analysis/MarkerSelection/rhesus_all.rds")$Grand
)

all_markers_df <- bind_rows(marker_list, .id = "species")

# ** Top 100 marker overlap comparison #####
# select top100 marker genes by delta expression fraction
marker_list_top100 <- all_markers_df %>% 
  arrange(-dpct) %>% 
  group_by(species,celltype) %>% 
  slice_head(n = 100)

# overlap across all species:
marker_overlap_full <- marker_list_top100 %>% 
  group_by(celltype, gene) %>% 
  summarise(species_list = list(species), 
            n_species = as.factor(length(species)))

# add cell type abbreviations
ct_dictionary = setNames(c("iPSCs","EE","NCI","SMC","CFib","EC","Hepa"),
                         c("iPSCs","early_ectoderm","neural_crest_I","smooth_muscle_cells","cardiac_fibroblasts", "epithelial_cells","hepatocytes"))

marker_overlap_full <- marker_overlap_full %>% 
  mutate(celltype_short = ct_dictionary[celltype],
         celltype_short = factor(celltype_short, levels = ct_dictionary))


# ** Biotype stratification ####

# Transcription factor annotation (Motif in image/Jaspar)
jasp_core <- data.table::fread("grep '>' /data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
                               col.names = c("motif_id","SYMBOL")) %>% 
  dplyr::mutate( motif_id = gsub(">", "", motif_id)) %>% 
  separate_rows(SYMBOL, sep="::") %>% 
  mutate(Evidence = "jaspar_core")
# motifs from Jaspar2022 vertebrate unvalidated
jasp_unval <- data.table::fread("grep '>' /data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/JASPAR2022_UNVALIDATED_non-redundant_pfms_jaspar.txt", 
                                col.names = c("motif_id", "SYMBOL")) %>% 
  dplyr::mutate( motif_id = gsub(">", "", motif_id)) %>% 
  separate_rows(SYMBOL, sep="::") %>% 
  mutate(Evidence = "jaspar_unvalidated")
# motifs from Madsen et al. 2018 (IMAGE)
image <- read_delim("/data/share/htp/perturb-seq/TF_selection_gRNA_design/TF_motifs/IMAGE/utils/Genename_Motif.txt", 
                    delim="\t", col_names = F) 
colnames(image) <- c("SYMBOL", "motif_id", "Evidence")

# human TFs with at least 1 annotated motif in the network
TFs <- bind_rows(jasp_core, jasp_unval, image) %>% 
  group_by(SYMBOL) %>% 
  dplyr::summarise(n_motif = length(motif_id)) %>% 
  dplyr::pull(SYMBOL) %>% 
  unique()

saveRDS(TFs, "zenodo/marker_gene_analysis/TFs.rds")

# protein coding / lncRNA annotation
gtf <- rtracklayer::readGFF("/data/ngs/genomes/refdata-cellranger-6.1.1/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
gene_type <- gtf %>% 
  select(gene_name, gene_type) %>% 
  distinct() %>% 
  filter(Biobase::isUnique(gene_name))
saveRDS(gene_type, "zenodo/marker_gene_analysis/gene_type.rds")

prot_coding_genes <- gene_type$gene_name[gene_type$gene_type == "protein_coding"]
lncRNAs <- gene_type$gene_name[gene_type$gene_type == "lncRNA"] 

# Add gene type annotations
marker_overlap_full <- left_join(marker_overlap_full, gene_type, by = c("gene" = "gene_name")) %>% 
  mutate(isTF = ifelse(gene %in% TFs, T, F))
saveRDS(marker_overlap_full, "zenodo/marker_gene_analysis/marker_overlap.rds")

marker_list <- lapply(marker_list, function(x){
  mutate(x, 
         ProtCoding = ifelse(gene %in% prot_coding_genes, T, F),
         TF = ifelse(gene %in% TFs, T, F),
         lncRNA = ifelse(gene %in% lncRNAs, T, F))
})
marker_df <- bind_rows(marker_list, .id = "species")

# 2) Rank Biased Overlap (RBO) ####
calculate_rbo <- function(SpecA, SpecB, Celltype, nMarkers, biotype){
  comp_list <- lapply(marker_list[c(SpecA, SpecB)], function(df){
    if(biotype != "all"){
      df <- filter(df, !!sym(biotype) == T)
    }
    
    df <- df %>% 
      filter(celltype == Celltype) %>% 
      #arrange(-lfc) %>% 
      arrange(-dpct) %>% 
      pull(gene) %>% 
      .[1:nMarkers]
    
    return(df)
  })
  
  rbo <- rbo::rbo_ext(comp_list[[1]], comp_list[[2]], p = 1)
  
  return(tibble(SpecA = SpecA, 
                SpecB = SpecB,
                Celltype = Celltype, 
                nMarkers = nMarkers,
                biotype = biotype,
                rbo = rbo))
}

# All possible combinations of species, celltype, number of markers and biotype
combinations <- expand.grid(
  SpecA = unique(marker_df$species),
  SpecB =  unique(marker_df$species),
  Celltype =  unique(marker_df$celltype),
  nMarkers = 100,
  biotype = c("all", "ProtCoding","TF","lncRNA")
) %>% 
  filter(SpecA != SpecB) %>% 
  mutate(biotype = as.character(biotype))

rbo_biotype <- pmap_dfr(combinations, calculate_rbo)

# Summarize for pairwise comparisons
rbo_biotype <- rbo_biotype %>% 
  rowwise() %>% 
  mutate(pair = paste0(sort(c(SpecA, SpecB)), collapse = "_"))

saveRDS(rbo_biotype, "zenodo/marker_gene_analysis/rbo_biotype.rds")

# 3) kNN classification #####
sce_filt <- readRDS("zenodo/marker_gene_analysis/sce_filt.rds")
PerfResCTData <- readRDS(file = "zenodo/marker_gene_analysis/F1_celltype.rds")

crossspecies.perform.ct.data <- PerfResCTData %>% 
  dplyr::filter(MarkerSpecies == "human" &
                  TrainSpecies == 'human' &
                  TrainSet == 'clone' &
                  TestSet == "clone" &
                  MarkerLevel == 'all' &
                  MarkerType == 'all' &
                  Biotype %in% c("proteincoding", "tf")) %>%
  dplyr::filter(grepl("human_29B5_vs_*",Set)) %>% 
  dplyr::filter(NoMarkers %in% as.character(seq(1:30)))

# ** macro averaged F1 score with bootstrapping ####

# remove cell types from clones below a cell number threshold
cell_numbers <- colData(sce_filt) %>% data.frame() %>% dplyr::count(individual,manual_annotation, name = "n_cells")

crossspecies.perform.ct.data2 <- crossspecies.perform.ct.data %>% 
  mutate(TrainClone = gsub("_vs.*", "", Set), TestClone = gsub(".*vs_","",Set)) %>% 
  left_join(cell_numbers, by = c("Celltype" = "manual_annotation", "TestClone" = "individual"))


# filter for cell numbers
crossspecies.perform.ct.data.filt <- crossspecies.perform.ct.data2 %>% 
  filter(n_cells > 20)

# save F1 per clone for protein coding genes
crossspecies.perform.ct.data.filt %>% 
  filter(Biotype == "proteincoding") %>% 
  saveRDS("zenodo/marker_gene_analysis/F1_per_clone.rds")


F1_macro_bootstrap_biotype <-  crossspecies.perform.ct.data.filt %>% 
  group_by(TestSpecies,Set, NoMarkers, Biotype) %>% 
  nest() %>%
  mutate(bootstrap_results = map(data, ~ bootstrap_f1(.x))) %>%
  unnest_wider(bootstrap_results) %>% 
  mutate(TestSpecies = factor(TestSpecies, levels = c("human", "orang", "cyno", "rhesus")),
         TestClone = gsub(".*_","",Set))

saveRDS(F1_macro_bootstrap_biotype, "zenodo/marker_gene_analysis/F1_macro_bootstrap.rds")