# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse));

option_list = list(
  make_option(c("--trainsce"), type="character", help="input rds train sce object with full path", metavar="character"),
  make_option(c("--trainsplitcol"), type="character", help="input train set column name", metavar="character"),
  make_option(c("--testsce"), type="character", help="input rds test sce object with full path", metavar="character"),
  make_option(c("--testsplitcol"), type="character", help="input test set column name", metavar="character"),
  make_option(c("--markers"), type="character", help="input rds markers object with full path", metavar="character"),
  make_option(c("--level"), type="character", help="consider celltype or germlayer", metavar="character"),
  make_option(c("--type"), type="character", help="unique, multiple or all per level", metavar="character"),
  make_option(c("--biotype"), type="character", help="input rds gene biotype object (gene names character vector) with full path", metavar="character"),
  make_option(c("--phenocol"), type="character", help="phenotype column name", metavar="character"),
  make_option(c("--nogenes"), type="character", help="input rds with vector number of top marker genes to use", metavar="character"),
  make_option(c("--k"), type="numeric", help="size of k for knn classification", default = 3, metavar="character"),
  make_option(c("--pc"), type="logical", help="parallel computation", default = T, metavar="character"),
  make_option(c("--workers"), type="numeric", help="number of cores", default = 2, metavar="character"),
  make_option(c("--outpath"), type="character", help="full output path", metavar="character"),
  make_option(c("--outname"), type="character", help="name for output .rds file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# SETTINGS ----------------------------------------------------------------

# source powsimR functions
source("/data/share/htp/EBgrant/Beate/scripts/foos.R")

# required R libraries
.libPaths()
suppressPackageStartupMessages(require(methods))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(FNN))
suppressPackageStartupMessages(require(caret))


# INPUT -------------------------------------------------------------------

# read in option arguments

# sce object
if(file.exists(opt$trainsce)){
  trainsce <- readRDS(file = opt$trainsce)
} else {
  message("Input train sce object could not be read. Abort.")
  invisible(gc())
  q()
}
if(file.exists(opt$testsce)){
  testsce <- readRDS(file = opt$testsce)
} else {
  message("Input test sce object could not be read. Abort.")
  invisible(gc())
  q()
}

print(trainsce)
print(testsce)

sce <- list('train' = trainsce, 'test' = testsce)

# markers object
if(file.exists(opt$markers)){
  markers <- readRDS(file = opt$markers)
} else {
  message("Input marker object could not be read. Abort.")
  invisible(gc())
  q()
}

# biotype object
if(file.exists(opt$biotype)){
  biotype <- readRDS(file = opt$biotype)
} else {
  message("Input biotype object could not be read. Abort.")
  invisible(gc())
  q()
}

# phenotype
phenocol <- opt$phenocol

# set to use for training/testing (clone / subsampling / all)
trainsplitcol <- opt$trainsplitcol
testsplitcol <- opt$testsplitcol

splitcol <- list('train' = trainsplitcol, 'test' = testsplitcol)

# level
ct_level <- opt$level

# type
type <- opt$type

# number of top marker genes to predict with
# nogenes <- opt$nogenes
if(file.exists(opt$nogenes)){
  nogenes <- readRDS(file = opt$nogenes)
} else {
  message("Input nogenes object could not be read. Abort.")
  invisible(gc())
  q()
}

# k for knn classification 
k_knn <- opt$k

# parallel computation
if(isTRUE(opt$pc)){
  BPPARAM <- BiocParallel::MulticoreParam(workers = opt$workers, tasks = 0L,
                                          stop.on.error = TRUE,
                                          progressbar = FALSE, RNGseed = NULL,
                                          timeout = 30L * 24L * 60L * 60L, exportglobals=TRUE,
                                          log = FALSE, threshold = "INFO", logdir = NA_character_,
                                          resultdir = NA_character_, jobname = "BPJOB",
                                          manager.hostname = NA_character_, manager.port = NA_integer_)
   
  } else {
    BPPARAM <- BiocParallel::bpparam('SerialParam')
}

# output object
if(dir.exists(opt$outpath)){
  Outpath <- opt$outpath
} else {
  message("Directory of output path cannot be found. Creating...")
  dir.create(file.path(opt$outpath))
}
Outname <- opt$outname

message(paste0("Output File of run_predict_marker(): ", 
               Outpath, "/", Outname, 
               ".rds"))

# CALCULATION--------------------------------------------------------------

# nogenes <- 1:4

# PerformanceList <- lapply(1:length(nogenes), function(i){
#   print(i)
#   predict_celltype_markers(scelist = sce,
#                            splitlist = splitcol,
#                            phenocol = phenocol,
#                            markers = markers,
#                            ct_level = ct_level,
#                            type = type,
#                            biotype = biotype,
#                            nogenes = nogenes[i],
#                            k_knn = k_knn)
# })

PerformanceList <- BiocParallel::bplapply(1:length(nogenes), function(i){
  predict_celltype_markers(scelist = sce,
                           splitlist = splitcol,
                           phenocol = phenocol,
                           markers = markers,
                           ct_level = ct_level,
                           type = type,
                           biotype = biotype,
                           nogenes = nogenes[i],
                           k_knn = k_knn)
}, BPPARAM = BPPARAM)
names(PerformanceList) <- nogenes

# str(PerformanceList)

Overall <- lapply(names(PerformanceList), function(i){
  PerformanceList[[i]]$Overall
})
names(Overall) <- names(PerformanceList)
Overall <- data.table::rbindlist(Overall, idcol = "nogenes")
CellType <- lapply(names(PerformanceList), function(i){
  PerformanceList[[i]]$CellType
})
names(CellType) <- names(PerformanceList)
CellType <- data.table::rbindlist(CellType, idcol = "nogenes")
MarkerData <- lapply(names(PerformanceList), function(i){
  PerformanceList[[i]]$Marker$Data
})
names(MarkerData) <- names(PerformanceList)
MarkerData <- data.table::rbindlist(MarkerData, idcol = "nogenes")
MarkerNumber <- sapply(names(PerformanceList), function(i){
  PerformanceList[[i]]$Marker$Number
}, simplify = F, USE.NAMES = T)
names(MarkerNumber) <- names(PerformanceList)
MarkerNumber <- tibble(nogenes = names(MarkerNumber),
                       nomarkers = unlist(unname(MarkerNumber)))

Out <- list("Overall" = Overall,
            "CellType" = CellType,
            "Marker" = list("Data" = MarkerData,
                            "Number" = MarkerNumber))

# PerformanceData <- predict_celltype_markers(sce = sce,
#                                             splitcol = splitcol,
#                                             phenocol = phenocol,
#                                             markers = markers,
#                                             ct_level = ct_level,
#                                             type = type,
#                                             biotype = biotype,
#                                             nogenes = nogenes,
#                                             k_knn = k_knn)



# OUTPUT ------------------------------------------------------------------


saveRDS(Out, file = paste0(Outpath, "/", Outname,  ".rds"))

# FIN ---------------------------------------------------------------------

sessionInfo()

gc()

q()
