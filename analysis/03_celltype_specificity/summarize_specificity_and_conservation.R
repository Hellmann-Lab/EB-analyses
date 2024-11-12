library(tidyverse)
library(ape)

# ....................................................................................................................................................... ----
# Presence-absence scoring of gene expression ####

# ** Detection summary tables ####
gene_detection <- list(
  human = readRDS("zenodo/cell_type_specificity/GeneDetection/Summary/human.rds")$Celltype,
  orang = readRDS("zenodo/cell_type_specificity/GeneDetection/Summary/orang.rds")$Celltype,
  cyno = readRDS("zenodo/cell_type_specificity/GeneDetection/Summary/cyno.rds")$Celltype,
  rhesus = readRDS("zenodo/cell_type_specificity/GeneDetection/Summary/rhesus.rds")$Celltype
) %>% bind_rows(.id = "species")

# add cell type abbreviations
ct_dictionary <- setNames(
  c("iPSCs","EE","AP","GC","NCI","NCII","Neu","MDI","SMC","MDII","CFib","MDnR","CPC","EEC","EC","Hepa","unass"),
  nm = c("iPSCs","early_ectoderm","astrocyte_progenitor","glial_cells","neural_crest_I","neural_crest_II","neurons","mesoderm_I","smooth_muscle_cells","mesoderm_II","cardiac_fibroblasts","mesoderm_noRibo","cardiac_progenitor_cells","early_epithelial_cells","epithelial_cells","hepatocytes","unassigned")
)
gene_detection$celltype_short <- ct_dictionary[as.character(gene_detection$Celltype)]
gene_detection$celltype_short <- factor(gene_detection$celltype_short, levels = ct_dictionary)

# ** Detection thresholds per species & cell type ####
detection_thresholds <- list(
  human = readRDS("zenodo/cell_type_specificity/GeneDetection/human_all.rds")$CT,
  orang = readRDS("zenodo/cell_type_specificity/GeneDetection/orang_all.rds")$CT,
  cyno = readRDS("zenodo/cell_type_specificity/GeneDetection/cyno_all.rds")$CT,
  rhesus = readRDS("zenodo/cell_type_specificity/GeneDetection/rhesus_all.rds")$CT
) %>% 
  lapply(function(species){
    lapply(species, function(celltype){
      data.frame(detection_cutoff = celltype$Cutoff)
    }) %>% bind_rows(.id = "celltype")
  }) %>% bind_rows(.id = "species")

# Select the maximum cut-off for each cell type
max_thresh <- detection_thresholds %>% 
  group_by(celltype) %>% 
  summarize(max_cutoff = max(detection_cutoff))

# ** Filter cell types present in all species ####
ct_allSpec <- c("iPSCs", "early_ectoderm","neural_crest_I","smooth_muscle_cells","cardiac_fibroblasts", "epithelial_cells","hepatocytes")
gene_detection_allSpec <- gene_detection %>% filter(Celltype %in% ct_allSpec)

# remove extra columns and add celltype-specific detection thresholds
gene_detection_allSpec <- gene_detection_allSpec %>% 
  rename("celltype" = "Celltype", "gene" = "Gene", "detection" = "Detection") %>% 
  select(gene,species,celltype,celltype_short,detection,Mean,SE) %>% 
  left_join(max_thresh) %>% 
  mutate(expressed = ifelse(detection >= max_cutoff,1,0),
         species = factor(species, levels = rev(c("human","orang","cyno","rhesus"))))

saveRDS(gene_detection_allSpec, "zenodo/cell_type_specificity/gene_detection_allSpec.RDS")

# ....................................................................................................................................................... ----
# Cell type specificity and conservation scores ####

# ** Summarize expressed number of species and celltypes per gene ####
pleio_df <- gene_detection_allSpec %>% 
  filter(expressed == T) %>% 
  group_by(gene,species) %>% 
  mutate(n_celltypes_per_species = length(unique(celltype)),
         expr_celltypes = list(celltype)) %>% 
  group_by(gene,celltype) %>% 
  mutate(n_species_per_celltype = sum(expressed),
         expr_species = list(species))

# ** Phylogeny based correction of conservation scores ####

# *** Load species tree ####
# The tree with phylogenetic distances was obtained from Bininda-Emonds et al. 2007 (10.1038/nature05634)
tree <- read.tree("zenodo/cell_type_specificity/mammaltree.txt") %>%
  drop.tip(.$tip.label[!.$tip.label %in% c("Homo_sapiens","Pongo_abelii","Macaca_mulatta", "Macaca_fascicularis")])
tree$tip.label <- c("rhesus","cyno",  "human","orang")

# *** Summarize branch lengths of all possible subtrees #### 
get_tree_weight <- function(tree, species){
  total_bl <- tree %>% unroot() %>% as_tibble() %>% drop_na(branch.length) %>%  pull(branch.length) %>% sum
  spec_bl <- tree %>% 
    unroot() %>% 
    keep.tip(.$tip.label[.$tip.label %in% species]) %>% 
    as_tibble() %>% 
    drop_na(branch.length) %>% 
    pull(branch.length) %>% 
    sum()
  return(spec_bl / total_bl)
}

species <- sort(c("cyno","rhesus","orang","human"))
tree_sum_species_combinations <- setNames(
  unlist(lapply(1:length(species), function(i) {combn(species, i, FUN = get_tree_weight, tree = tree)})),
  unlist(lapply(1:length(species), function(i) {combn(species, i, FUN = paste, collapse = "_")})))

# ** Summarize expression conservation scores ####

# Per gene & celltype: sum of branch lengths
pleio_df <- pleio_df %>%
  rowwise() %>% 
  mutate(tree_size_celltype = tree_sum_species_combinations[paste(sort(as.character(unlist(expr_species))), collapse = "_")])

# Per gene: conservation score
conservation_phylo_df <- pleio_df %>% 
  dplyr::select(gene, celltype, tree_size_celltype) %>%
  distinct() %>% 
  group_by(gene) %>% 
  summarize(conservation_phylo = mean(tree_size_celltype))

# Pleiotropic degree (per gene & species) + conservation (per gene) 
pleio_df_summarized <- pleio_df %>% 
  group_by(species, gene) %>% 
  summarize(PD = unique(n_celltypes_per_species),
            mean_expression = mean(Mean)) %>% 
  left_join(conservation_phylo_df) %>% 
  mutate(species = factor(species, levels = c("human","orang","cyno","rhesus")))


# ....................................................................................................................................................... ----
# Sequence conservation (Zoonomia phyloP & phastCons) ####

# Summarized to measure constraint per gene in Sullivan et al. 2023 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10259825/) as 
# 1) fraction of CDS bases with mammalian phyloP ≥ 2.27 (fracCdsCons) 
# 2) fraction of CDS bases with primate phastCons ≥ 0.96 (fracConsPr)

fracCdsCons_zoonomia <- readxl::read_excel("NIHMS1897004-supplement-SuppTab_s11-15.xlsx", 
                                           sheet = "Table S14", skip = 1) %>% 
  mutate(gene_name = gsub("/","",gene_name))

# Add scores:
pleio_phast <- pleio_df_summarized %>% 
  left_join(fracCdsCons_zoonomia, by = c("gene" = "gene_name"))

saveRDS(pleio_phast, "zenodo/cell_type_specificity/ct_specificity_summary.RDS")


