
# SETTINGS ----------------------------------------------------------------

# R packages
require(tidyverse)
require(SingleCellExperiment)
require(Seurat)
require(scran)
require(scuttle)
require(quantreg)
require(MASS)
require(fastDummies)
source("scripts/foos.R")

BPPARAM <- BiocParallel::MulticoreParam(workers = 8, tasks = 0L,
                                        stop.on.error = TRUE,
                                        progressbar = FALSE, RNGseed = NULL,
                                        timeout = 30L * 24L * 60L * 60L, exportglobals=TRUE,
                                        log = FALSE, threshold = "INFO", logdir = NA_character_,
                                        resultdir = NA_character_, jobname = "BPJOB",
                                        manager.hostname = NA_character_, manager.port = NA_integer_)

verbose <- TRUE

# INPUT -------------------------------------------------------------------

# Read in annotated Seurat object
seu_list <- readRDS("zenodo/cell_type_assignment/seu_list.RDS")

# Convert to SingleCellExperiment Object
sce_list <- sapply(names(seu_list), function(seu){
  Seurat::as.SingleCellExperiment(seu_list[[seu]])
}, USE.NAMES = T, simplify = F)

sce_all <- cbind(sce_list[[1]], sce_list[[2]], sce_list[[3]], sce_list[[4]])

# columns: species, individual for clone, manual_annotation for cell type
table(sce_all$manual_annotation, sce_all$species)

# abandoned prior filtering based on species and presence of ct, only remove unassigned!
sce_all <- S4Vectors::subset(sce_all, , manual_annotation != "unassigned")

# GENE AND CELL FILTERING (powsimR pipeline)

sce <- sce_all

# kick out empty cells and keep only nonsparse gene
totalC <- ncol(SingleCellExperiment::counts(sce))
totalG <- nrow(SingleCellExperiment::counts(sce))
nonsparse.gene.cut <- dplyr::case_when(totalC <= 1000 ~ 1,
                                       totalC > 1000 & totalC <= 5000 ~ 2,
                                       totalC > 5000 ~ 5)
DetectC <- Matrix::colSums(SingleCellExperiment::counts(sce) > 0) > 10
DetectG <- Matrix::rowSums(SingleCellExperiment::counts(sce) > 0) >= nonsparse.gene.cut
detectC <- sum(DetectC)
detectG <- sum(DetectG)
sce <- sce[DetectG,DetectC]

if(verbose) {
  message(paste0("The provided SingleCellExperiment has ",
                 detectC, " out of ",
                 totalC, " cells and ",
                 detectG, " out of ",
                 totalG, " genes with nonsparse expression."))
}

### CELLS
group <- paste(SingleCellExperiment::colData(sce)[, "manual_annotation"],
               SingleCellExperiment::colData(sce)[, "species"],
               sep = "_")
# MAD calculation based on seq depth, no of genes, spike-ins etc
cell.stats <- scuttle::perCellQCMetrics(sce,
                                        percent_top = NULL,
                                        flatten = TRUE)
cell.filter <- sapply(colnames(cell.stats), function(i) {
  scuttle::isOutlier(cell.stats[,i],
                     nmads = 3, 
                     type = "both",
                     log = TRUE, 
                     batch = group)
}, USE.NAMES = TRUE, simplify = TRUE)
cell.filter <- cbind(cell.filter,
                     final = rowSums(cell.filter, na.rm = T) > 0L)
colnames(cell.filter) <- paste(colnames(cell.filter), "outlier", sep = "_")
cell.annot <- cbind(cell.stats, cell.filter)
### GENES
gene.stats.L <- sapply(unique(group), function(j) {
  batch.tmp <- group == j
  sce_set <- subset(sce, , batch.tmp)
  gene_set <- data.frame(scuttle::perFeatureQCMetrics(sce_set, BPPARAM = BPPARAM))
  colnames(gene_set) <- paste(colnames(gene_set), j, sep = "_")
  gene_set[, grepl(pattern = "detected", colnames(gene_set))]
}, simplify = FALSE, USE.NAMES = TRUE)
gene.stats <- dplyr::bind_cols(gene.stats.L)
colnames(gene.stats) <- paste(colnames(gene.stats), "Detection", sep = "_")
gene.filter <- cbind(gene.stats,
                     Dropout = !apply(gene.stats, 1,
                                      function(x) any(x > 0.01 * 100)),
                     MinExpr = apply(gene.stats, 1,
                                     function(x) all(x > 0.01 * 100)))

### Append to colData and rowData
sce.filt <- sce
SummarizedExperiment::colData(sce.filt) <- cbind(SingleCellExperiment::colData(sce.filt),
                                                 S4Vectors::DataFrame(cell.annot))
SummarizedExperiment::colData(sce.filt)[, "filter"] <- !SummarizedExperiment::colData(sce.filt)[, "final_outlier"]

SummarizedExperiment::rowData(sce.filt) <- cbind(SingleCellExperiment::rowData(sce.filt),
                                                 S4Vectors::DataFrame(gene.filter))
SummarizedExperiment::rowData(sce.filt)[, "filter"] <- !SummarizedExperiment::rowData(sce.filt)[, "Dropout"]

if(verbose) {
  message(paste0(sum(cell.annot[,"final_outlier"]), " out of ",
                 nrow(cell.annot),
                 " cells were determined to be outliers."))
  message(paste0(sum(!gene.filter[,"Dropout"]), " out of ", nrow(gene.filter),
                 " genes were determined to be expressed."))
}

sce.clean <- sce.filt[SummarizedExperiment::rowData(sce.filt)[, "filter"], SummarizedExperiment::colData(sce.filt)[, "filter"]]

sce.clean
table(sce.clean$manual_annotation, sce.clean$species)
saveRDS(sce.clean, file="data/input/sce_clean.rds")


#  SPLITTING OF DATA ------------------------------------------------------

sce.clean <- readRDS("data/input/sce_clean.rds")
table(sce.clean$species, sce.clean$manual_annotation)

pdat <- data.frame(CellType = sce.clean$manual_annotation,
                   Species = sce.clean$species)
ggplot(data = pdat, aes(x = CellType, fill = Species)) +
  geom_bar(stat="count", position=position_dodge(), color = "black") +
  scale_y_log10() +
  coord_flip() +
  theme_minimal(base_size = 16)

#### DOWNSAMPLE AND k ITERATIONS #####

# downsample to even numbers and assign iterations for testing per species
sce_hum_filt <- run_downsample_iter(sce = sce.clean,  maxc = 4700, iter = 4,
                                    SpeciesCol = "human", 
                                    PhenotypeCol = "manual_annotation")
table(sce_hum_filt$manual_annotation, sce_hum_filt$TestIterSet)
saveRDS(sce_hum_filt, file="data/input/sce_hum_iter.rds")
sce_cyno_filt <- run_downsample_iter(sce = sce.clean, maxc = 4700, iter = 4,
                                     SpeciesCol = "cyno",
                                     PhenotypeCol = "manual_annotation")
table(sce_cyno_filt$manual_annotation, sce_cyno_filt$TestIterSet)
saveRDS(sce_cyno_filt, file="data/input/sce_cyno_iter.rds")
sce_rhesus_filt <- run_downsample_iter(sce = sce.clean, maxc = 4700, iter = 4,
                                       SpeciesCol = "rhesus",
                                       PhenotypeCol = "manual_annotation")
table(sce_rhesus_filt$manual_annotation, sce_rhesus_filt$TestIterSet)
saveRDS(sce_rhesus_filt, file="data/input/sce_rhesus_iter.rds")
sce_orang_filt <- run_downsample_iter(sce = sce.clean, maxc = 4700, iter = 4,
                                      SpeciesCol = "orang",
                                      PhenotypeCol = "manual_annotation")
table(sce_orang_filt$manual_annotation, sce_orang_filt$TestIterSet)
saveRDS(sce_orang_filt, file="data/input/sce_orang_iter.rds")

#### CLONE-SPECIFIC DATA SETS #####

# split by clone id per species
# human
sce_hum_29B5 <- subset(sce.clean, , sce.clean$individual == "human_29B5")
table(sce_hum_29B5$manual_annotation)
saveRDS(sce_hum_29B5, file="data/input/sce_hum_29B5.rds")
sce_hum_63Ab2.2 <- subset(sce.clean, , sce.clean$individual == "human_63Ab2.2")
table(sce_hum_63Ab2.2$manual_annotation)
saveRDS(sce_hum_63Ab2.2, file="data/input/sce_hum_63Ab2.2.rds")
# cyno
sce_cyno_56A1 <- subset(sce.clean, , sce.clean$individual == "cyno_56A1")
table(sce_cyno_56A1$manual_annotation)
saveRDS(sce_cyno_56A1, file="data/input/sce_cyno_56A1.rds")
sce_cyno_56B1 <- subset(sce.clean, , sce.clean$individual == "cyno_56B1")
table(sce_cyno_56B1$manual_annotation)
saveRDS(sce_cyno_56B1, file="data/input/sce_cyno_56B1.rds")
sce_cyno_82A3 <- subset(sce.clean, , sce.clean$individual == "cyno_82A3")
table(sce_cyno_82A3$manual_annotation)
saveRDS(sce_cyno_82A3, file="data/input/sce_cyno_82A3.rds")
# rhesus 
sce_rhesus_83Ab1.1 <- subset(sce.clean, , sce.clean$individual == "rhesus_83Ab1.1")
table(sce_rhesus_83Ab1.1$manual_annotation)
saveRDS(sce_rhesus_83Ab1.1, file="data/input/sce_rhesus_83Ab1.1.rds")
sce_rhesus_83D1 <- subset(sce.clean, , sce.clean$individual == "rhesus_83D1")
table(sce_rhesus_83D1$manual_annotation)
saveRDS(sce_rhesus_83D1, file="data/input/sce_rhesus_83D1.rds")
sce_rhesus_87B1 <- subset(sce.clean, , sce.clean$individual == "rhesus_87B1")
table(sce_rhesus_87B1$manual_annotation)
saveRDS(sce_rhesus_87B1, file="data/input/sce_rhesus_87B1.rds")
# orang 
sce_orang_68A20 <- subset(sce.clean, , sce.clean$individual == "orang_68A20")
table(sce_orang_68A20$manual_annotation)
saveRDS(sce_orang_68A20, file="data/input/sce_orang_68A20.rds")
sce_orang_69A1 <- subset(sce.clean, , sce.clean$individual == "orang_69A1")
table(sce_orang_69A1$manual_annotation)
saveRDS(sce_orang_69A1, file="data/input/sce_orang_69A1.rds")

# GLMGAMPOI ----------------------------------------------------


##### SETUP #####

# downsampled iterations
infiles <- c("data/input/sce_hum_iter.rds",
             "data/input/sce_cyno_iter.rds",
             "data/input/sce_rhesus_iter.rds",
             "data/input/sce_orang_iter.rds")
species <- c('human', 'cyno', 'rhesus', 'orang')
iter <- c(1:4, NA)
IterCol <- "TestIterSet"
PhenotypeCol <- "manual_annotation"
outpath <- "data/output/glmGamPoi/"
mem <- 25000
verbose <- TRUE

glm.params.iter <- tidyr::expand_grid(tibble(infiles, 
                                             species,
                                             IterCol,
                                             PhenotypeCol,
                                             outpath,
                                             mem,
                                             verbose),
                                      iter = iter)

write.table(glm.params.iter, 
            file ="data/output/glmGamPoi/param_iter.txt",
            sep = "\t", row.names = T, quote = F, col.names = F, na='none')

# clone-wise data
infiles <- c("data/input/sce_hum_29B5.rds",
             "data/input/sce_hum_63Ab2.2.rds",
             "data/input/sce_cyno_56A1.rds",
             "data/input/sce_cyno_56B1.rds",
             "data/input/sce_cyno_82A3.rds",
             "data/input/sce_rhesus_83Ab1.1.rds",
             "data/input/sce_rhesus_83D1.rds",
             "data/input/sce_rhesus_87B1.rds",
             "data/input/sce_orang_68A20.rds",
             "data/input/sce_orang_69A1.rds")
species <- c('human', 'human', 
             'cyno', 'cyno', 'cyno', 
             'rhesus',  'rhesus',  'rhesus', 
             'orang', 'orang')
clone <- c("29B5", "63Ab2.2", 
           "56A1", "56B1", "82A3",
           "83Ab1.1", "83D1", "87B1",
           "68A20", "69A1")
iter <- paste(species, clone, sep = "_")
IterCol <- "individual"
PhenotypeCol <- "manual_annotation"
outpath <- "data/output/glmGamPoi/"
mem <- 25000
verbose <- TRUE

glm.params.clone <- tibble(infiles, 
                           species,
                           IterCol,
                           PhenotypeCol,
                           outpath,
                           mem,
                           verbose,
                           iter)

write.table(glm.params.clone, 
            file ="data/output/glmGamPoi/param_clone.txt",
            sep = "\t", row.names = T, quote = F, col.names = F, na='none')

# combined
glm.params.all <- rbind(glm.params.iter, glm.params.clone)
write.table(glm.params.all, 
            file ="data/output/glmGamPoi/param_all.txt",
            sep = "\t", row.names = T, quote = F, col.names = F, na='none')


##### SCRIPTS #####

# Rscript
# slurm script
# bash wrapper around slurm script with Rscript


# GENE DETECTION ------------------------------------------------

##### SETUP #####
# downsampled iterations
sce.files <- c("data/input/sce_hum_iter.rds",
               "data/input/sce_cyno_iter.rds",
               "data/input/sce_rhesus_iter.rds",
               "data/input/sce_orang_iter.rds")
species <- c('human', 'cyno', 'rhesus', 'orang')
itercol <- "TestIterSet"
iternumber <- c(1:4, NA)
phenocol <- "manual_annotation"
quantile.se <- 0.9
quantile.mean <- 0.05
pmax <- 0.5
pthreshold <- 0.999
pc <- TRUE
workers <- 4
mem <- 60000
verbose <- TRUE
outpath <- "data/output/GeneDetection/"

detect.params.iter <- tidyr::expand_grid(tibble(sce.files, 
                                                species,
                                                itercol,
                                                phenocol,
                                                quantile.se,
                                                quantile.mean,
                                                pmax,
                                                pthreshold,
                                                outpath,
                                                pc,
                                                workers,
                                                mem,
                                                verbose),
                                         iter = iternumber) %>% 
  dplyr::mutate(itername = ifelse(is.na(iter), "all", iter)) %>% 
  dplyr::mutate(ggp.files = paste0("data/output/glmGamPoi/",
                                   species,
                                   "_",
                                   itername,
                                   ".rds")) %>% 
  dplyr::select(sce.files, ggp.files, phenocol, itercol, iter, quantile.se, quantile.mean, pmax, pthreshold, outpath, species, pc, workers, mem, verbose)

write.table(detect.params.iter, 
            file ="data/output/GeneDetection/param_clone.txt",
            sep = "\t", row.names = T, quote = F, col.names = F, na='none')


# clone-wise data
sce.files <- c("data/input/sce_hum_29B5.rds",
               "data/input/sce_hum_63Ab2.2.rds",
               "data/input/sce_cyno_56A1.rds",
               "data/input/sce_cyno_56B1.rds",
               "data/input/sce_cyno_82A3.rds",
               "data/input/sce_rhesus_83Ab1.1.rds",
               "data/input/sce_rhesus_83D1.rds",
               "data/input/sce_rhesus_87B1.rds",
               "data/input/sce_orang_68A20.rds",
               "data/input/sce_orang_69A1.rds")
species <- c('human', 'human', 
             'cyno', 'cyno', 'cyno', 
             'rhesus',  'rhesus',  'rhesus', 
             'orang', 'orang')
clone <- c("29B5", "63Ab2.2", 
           "56A1", "56B1", "82A3",
           "83Ab1.1", "83D1", "87B1",
           "68A20", "69A1")
iter <- paste(species, clone, sep = "_")
itercol <- "individual"
phenocol <- "manual_annotation"
quantile.se <- 0.9
quantile.mean <- 0.05
pmax <- 0.5
pthreshold <- 0.999
pc <- TRUE
workers <- 4
mem <- 60000
verbose <- TRUE
outpath <- "data/output/GeneDetection/"

detect.params.clone <- tibble(sce.files, 
                              species,
                              itercol,
                              phenocol,
                              quantile.se,
                              quantile.mean,
                              pmax,
                              pthreshold,
                              outpath,
                              pc,
                              workers,
                              mem,
                              verbose,
                              iter) %>% 
  dplyr::mutate(ggp.files = paste0("data/output/glmGamPoi/",
                                   species,
                                   "_",
                                   iter,
                                   ".rds")) %>% 
  dplyr::select(sce.files, ggp.files, phenocol, itercol, iter, quantile.se, quantile.mean, pmax, pthreshold, outpath, species, pc, workers, mem, verbose)

write.table(detect.params.clone, 
            file ="data/output/GeneDetection/param_clone.txt",
            sep = "\t", row.names = T, quote = F, col.names = F, na='none')

# combined
detect.params.all <- rbind(detect.params.iter, detect.params.clone)

write.table(detect.params.all, 
            file ="data/output/GeneDetection/param_all.txt",
            sep = "\t", row.names = T, quote = F, col.names = F, na='none')


##### SCRIPTS #####


# GENE DETECTION SUMMARY --------------------------------------------------

# cyno
FitFiles.Iter <- c("data/output/GeneDetection/cyno_1.rds",
                   "data/output/GeneDetection/cyno_2.rds",
                   "data/output/GeneDetection/cyno_3.rds",
                   "data/output/GeneDetection/cyno_4.rds")
names(FitFiles.Iter) <- c("cyno_1", "cyno_2", "cyno_3", "cyno_4")
FitFiles.All <- c("data/output/GeneDetection/cyno_all.rds")
FitFiles.Clone <- c("data/output/GeneDetection/cyno_cyno_56A1.rds",
                    "data/output/GeneDetection/cyno_cyno_82A3.rds")
names(FitFiles.Clone) <- c("cyno_cyno_56A1", "cyno_cyno_82A3")
DetectResList <- list("Iteration" = sapply(names(FitFiles.Iter), function(i){
  readRDS(file = FitFiles.Iter[[i]])
}, USE.NAMES = T, simplify = F),
"Clone" = sapply(names(FitFiles.Clone), function(i){
  readRDS(file = FitFiles.Clone[[i]])
}, USE.NAMES = T, simplify = F),
"All" = readRDS(file = FitFiles.All)
) 

TmpC <- process_gene_detection(DetectResList = DetectResList, 
                               Type = "Clone",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))
TmpK <- process_gene_detection(DetectResList = DetectResList, 
                               Type = "Subsampling",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))

TmpGrand <- dplyr::left_join(TmpK$Grand, 
                             TmpC$Grand %>% dplyr::select(Gene, Expressed_c, Reproducible_Beta_c),
                             by = "Gene")

TmpCT <- dplyr::left_join(TmpK$Celltype, 
                          TmpC$Celltype %>% dplyr::select(Celltype, Gene, Expressed_c, Reproducible_Beta_c),
                          by = c("Celltype", "Gene"))

Tmp <- list("Grand" = TmpGrand,
            "Celltype" = TmpCT,
            "Setings" = list("Clone" = TmpC$Settings, 
                             "Subsampling" = TmpK$Settings))

saveRDS(Tmp, file = "data/output/GeneDetection/Summary/cyno.rds")

#rhesus
FitFiles.Iter <- c("data/output/GeneDetection/rhesus_1.rds",
                   "data/output/GeneDetection/rhesus_2.rds",
                   "data/output/GeneDetection/rhesus_3.rds",
                   "data/output/GeneDetection/rhesus_4.rds")
names(FitFiles.Iter) <- c("rhesus_1", "rhesus_2", "rhesus_3", "rhesus_4")
FitFiles.All <- c("data/output/GeneDetection/rhesus_all.rds")
FitFiles.Clone <- c("data/output/GeneDetection/rhesus_rhesus_83Ab1.1.rds",
                    "data/output/GeneDetection/rhesus_rhesus_83D1.rds",
                    "data/output/GeneDetection/rhesus_rhesus_87B1.rds")
names(FitFiles.Clone) <- c("rhesus_rhesus_83Ab1.1", "rhesus_rhesus_83D1", "rhesus_rhesus_87B1")
DetectResList <- list("Iteration" = sapply(names(FitFiles.Iter), function(i){
  readRDS(file = FitFiles.Iter[[i]])
}, USE.NAMES = T, simplify = F),
"Clone" = sapply(names(FitFiles.Clone), function(i){
  readRDS(file = FitFiles.Clone[[i]])
}, USE.NAMES = T, simplify = F),
"All" = readRDS(file = FitFiles.All))

TmpC <- process_gene_detection(DetectResList = DetectResList, 
                               Type = "Clone",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))
TmpK <- process_gene_detection(DetectResList = DetectResList,
                               Type = "Subsampling",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))

TmpGrand <- dplyr::left_join(TmpK$Grand, 
                             TmpC$Grand %>% dplyr::select(Gene, Expressed_c, Reproducible_Beta_c),
                             by = "Gene")

TmpCT <- dplyr::left_join(TmpK$Celltype, 
                          TmpC$Celltype %>% dplyr::select(Celltype, Gene, Expressed_c, Reproducible_Beta_c),
                          by = c("Celltype", "Gene"))

Tmp <- list("Grand" = TmpGrand, 
            "Celltype" = TmpCT,
            "Setings" = list("Clone" = TmpC$Settings, 
                             "Subsampling" = TmpK$Settings))

saveRDS(Tmp, file = "data/output/GeneDetection/Summary/rhesus.rds")

#orang
FitFiles.Iter <- c("data/output/GeneDetection/orang_1.rds",
                   "data/output/GeneDetection/orang_2.rds",
                   "data/output/GeneDetection/orang_3.rds",
                   "data/output/GeneDetection/orang_4.rds")
names(FitFiles.Iter) <- c("orang_1", "orang_2", "orang_3", "orang_4")
FitFiles.All <- c("data/output/GeneDetection/orang_all.rds")
FitFiles.Clone <- c("data/output/GeneDetection/orang_orang_68A20.rds",
                    "data/output/GeneDetection/orang_orang_69A1.rds")
names(FitFiles.Clone) <- c("orang_orang_68A20", "orang_orang_69A1")
DetectResList <- list("Iteration" = sapply(names(FitFiles.Iter), function(i){
  readRDS(file = FitFiles.Iter[[i]])
}, USE.NAMES = T, simplify = F),
"Clone" = sapply(names(FitFiles.Clone), function(i){
  readRDS(file = FitFiles.Clone[[i]])
}, USE.NAMES = T, simplify = F),
"All" = readRDS(file = FitFiles.All))

TmpC <- process_gene_detection(DetectResList = DetectResList,
                               Type = "Clone",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))
TmpK <- process_gene_detection(DetectResList = DetectResList,
                               Type = "Subsampling",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))

TmpGrand <- dplyr::left_join(TmpK$Grand, 
                             TmpC$Grand %>% dplyr::select(Gene, Expressed_c, Reproducible_Beta_c),
                             by = "Gene")

TmpCT <- dplyr::left_join(TmpK$Celltype, 
                          TmpC$Celltype %>% dplyr::select(Celltype, Gene, Expressed_c, Reproducible_Beta_c),
                          by = c("Celltype", "Gene"))

Tmp <- list("Grand" = TmpGrand, 
            "Celltype" = TmpCT,
            "Setings" = list("Clone" = TmpC$Settings,
                             "Subsampling" = TmpK$Settings))

saveRDS(Tmp, file = "data/output/GeneDetection/Summary/orang.rds")

#human
FitFiles.Iter <- c("data/output/GeneDetection/human_1.rds",
                   "data/output/GeneDetection/human_2.rds",
                   "data/output/GeneDetection/human_3.rds",
                   "data/output/GeneDetection/human_4.rds")
names(FitFiles.Iter) <- c("human_1", "human_2", "human_3", "human_4")
FitFiles.All <- c("data/output/GeneDetection/human_all.rds")
FitFiles.Clone <- c("data/output/GeneDetection/human_human_29B5.rds",
                    "data/output/GeneDetection/human_human_63Ab2.2.rds")
names(FitFiles.Clone) <- c("human_human_29B5", "human_human_63Ab2.2")
DetectResList <- list("Iteration" = sapply(names(FitFiles.Iter), function(i){
  readRDS(file = FitFiles.Iter[[i]])
}, USE.NAMES = T, simplify = F),
"Clone" = sapply(names(FitFiles.Clone), function(i){
  readRDS(file = FitFiles.Clone[[i]])
}, USE.NAMES = T, simplify = F),
"All" = readRDS(file = FitFiles.All))

TmpC <- process_gene_detection(DetectResList = DetectResList, 
                               Type = "Clone",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))
TmpK <- process_gene_detection(DetectResList = DetectResList, 
                               Type = "Subsampling",
                               Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                               Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                    "CT" = c("Subsampling" = 'within', "Clone" = 'any')))

TmpGrand <- dplyr::left_join(TmpK$Grand, 
                             TmpC$Grand %>% dplyr::select(Gene, Expressed_c, Reproducible_Beta_c),
                             by = "Gene")

TmpCT <- dplyr::left_join(TmpK$Celltype, 
                          TmpC$Celltype %>% dplyr::select(Celltype, Gene, Expressed_c, Reproducible_Beta_c),
                          by = c("Celltype", "Gene"))

Tmp <- list("Grand" = TmpGrand,
            "Celltype" = TmpCT,
            "Setings" = list("Clone" = TmpC$Settings, 
                             "Subsampling" = TmpK$Settings))

saveRDS(Tmp, file = "data/output/GeneDetection/Summary/human.rds")

# GENE CUTOFF COMPARISON --------------------------------------------------

GeneDetectResFiles <- list.files(path = "data/output/GeneDetection/", pattern = "*.rds$", full.names = TRUE)
GeneDetectResNames <- list.files(path = "data/output/GeneDetection/", pattern = "*.rds$")
GeneDetectResNames <- gsub(GeneDetectResNames, pattern = ".rds", replacement = "")
names(GeneDetectResFiles) <- GeneDetectResNames

CutoffData.L <- sapply(GeneDetectResNames, function(i){
  GeneDetectRes <- readRDS(file = GeneDetectResFiles[i])
  process_cutoff(GeneDetectRes = GeneDetectRes)
}, simplify = F, USE.NAMES = T)

GrandCutoff <- sapply(names(CutoffData.L), function(i){
  CutoffData.L[[i]]$Grand
}, USE.NAMES = TRUE, simplify = T)
GrandCutoffData <- tibble::tibble(Detection = GrandCutoff,
                                  Set = names(GrandCutoff)) %>% 
  tidyr::separate(col = Set, into = c("Species", "Sample", "Clone"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Name = dplyr::case_when(is.na(Clone) ~ Sample,
                                        !is.na(Clone) ~ Clone),
                Type = dplyr::case_when(is.na(Clone) ~ "Subsampling",
                                        !is.na(Clone) ~ "Clone"))

ggplot(data = GrandCutoffData %>% filter(Type == "Subsampling"),
       aes(x = Sample, y = Detection, fill = Species )) +
  geom_bar(stat="identity", position=position_dodge(), color = "black") +
  facet_wrap(~ Species, scales = "free_x" ) + 
  ylim(c(0,0.1)) + 
  theme_minimal(base_size = 16)

ggplot(data = GrandCutoffData %>% filter(Type == "Clone"),
       aes(x = Clone, y = Detection, fill = Species )) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~ Species, scales = "free_x" ) + 
  ylim(c(0,0.1)) + 
  theme_minimal()

CelltypeCutoff <- sapply(names(CutoffData.L), function(i){
  CutoffData.L[[i]]$Celltype
}, USE.NAMES = TRUE, simplify = F)
CelltypeCutoffData <- data.table::rbindlist(CelltypeCutoff, idcol = "Set") %>% 
  tibble::as_tibble() %>% 
  tidyr::separate(col = Set, into = c("Species", "Sample", "Clone"), sep = "_", remove = FALSE) %>% 
  dplyr::mutate(Name = dplyr::case_when(is.na(Clone) ~ Sample,
                                        !is.na(Clone) ~ Clone),
                Type = dplyr::case_when(is.na(Clone) ~ "Subsampling",
                                        !is.na(Clone) ~ "Clone")) %>% 
  dplyr::rename(Detection = min.spct, Celltype = celltype)

p <- ggplot(data = CelltypeCutoffData %>% filter(Type == "Subsampling" & Name == "all" & Celltype %in% c("cardiac_fibroblasts",
                                                                                                         "early_ectoderm",
                                                                                                         "epithelial_cells",
                                                                                                         "hepatocytes",
                                                                                                         "iPSCs",
                                                                                                         "neural_crest_I",
                                                                                                         "smooth_muscle_cells")),
            aes(x = Celltype, y = Detection, fill = Celltype)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black") +
  scale_fill_manual(values= c("cardiac_fibroblasts" = "tomato3",
                              "early_ectoderm" = "lightblue",
                              "epithelial_cells" = "palegreen3",
                              "hepatocytes" = "darkgreen",
                              "iPSCs" = "grey40",
                              "neural_crest_I" = "dodgerblue4",
                              "smooth_muscle_cells" = "tomato1")) +
  facet_wrap(~Species) +
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "", y="Detection Threshold")
p

ggsave(p, file = "report/detect_threshold_ct_species.png", dpi = 600, width = 250, height = 150, units = "mm")


ctset <- c("cardiac_fibroblasts",
           "early_ectoderm",
           "epithelial_cells",
           "hepatocytes",
           "iPSCs",
           "neural_crest_I",
           "smooth_muscle_cells")

glset <- c("mesoderm",
           "ectoderm",
           "endoderm",
           "endoderm",
           'pluripotent',
           "ectoderm",
           "mesoderm")

ggplot(data = CelltypeCutoffData %>% filter(Type == "Clone"),
       aes(x = Clone, y = Detection, fill = Celltype )) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~ Species, scales = "free_x" ) + 
  # coord_flip() +
  theme_minimal(base_size = 16)


unique(CelltypeCutoffData$Celltype)


i = 1


GeneDetectRes <- readRDS(file = GeneDetectResFiles[1])
Dat <- GeneDetectRes[["CT"]][["smooth_muscle_cells"]][["Data"]] %>% 
  dplyr::filter(Detection < 0.5)

ggplot(data = Dat, aes(x = Detection, y = Expressed)) + geom_point() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) +
  scale_y_continuous(limits=c(0,1), breaks = c(0,1)) +
  theme_bw(base_size = 14) +
  labs(y = "Putative \nexpression status") 

ggsave(file = "report/lr_fit.png", dpi = 600, width = 100, height = 80, units = "mm")