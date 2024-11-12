# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse));

option_list = list(
  make_option(c("--sce"), type="character", help="input rds sce object with full path", metavar="character"),
  make_option(c("--ggp"), type="character", help="input rds glmgampoi result object with full path", metavar="character"),
  make_option(c("--phenotype"), type="character", help="phenotype column name", metavar="character"),
  make_option(c("--iteration"), type="character", help="iteration set column name", metavar="character"),
  make_option(c("--iter"), type="numeric", help="iteration number, can also be NA then takes all", metavar="character"),
  make_option(c("--qmean"), type="numeric", help="quantile for mean cutoff (lowly expressed genes)", default = 0.1, metavar="character"),
  make_option(c("--qse"), type="numeric", help="quantile for standard error cutoff (unreliably estimated genes)", default = 0.75, metavar="character"),
  make_option(c("--pmax"), type="numeric", help="maximum detection to consider in logisitic regression", default = 0.5, metavar="character"),
  make_option(c("--pthreshold"), type="numeric", help="probability for logisitic regression", default = 0.999, metavar="character"),
  make_option(c("--pc"), type="logical", help="parallelisation or not", default = TRUE, metavar="character"),
  make_option(c("--workers"), type="numeric", help="number of workers for parallel computation", default = 2, metavar="character"),
  make_option(c("--outpath"), type="character", help="full output path", metavar="character"),
  make_option(c("--outname"), type="character", help="name for output .rds file", metavar="character"),
  make_option(c("--verbose"), type="logical", help="verbosity", default = TRUE, metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

# SETTINGS ----------------------------------------------------------------

# source powsimR functions
source("analysis/03_celltype_specificity/foos.R")

# required R libraries
.libPaths()
suppressPackageStartupMessages(require(methods))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(scuttle))
suppressPackageStartupMessages(require(BiocParallel))
suppressPackageStartupMessages(require(glmGamPoi))
suppressPackageStartupMessages(require(logistf))


# INPUT -------------------------------------------------------------------

# read in option arguments

# sce object
if(file.exists(opt$ggp)){
  sce <- readRDS(file = opt$sce)
} else {
  message("Input sce object could not be read. Abort.")
  invisible(gc())
  q()
}

print(sce)

# ggp object
if(file.exists(opt$ggp)){
  ggp <- readRDS(file = opt$ggp)
} else {
  message("Input ggp object could not be read. Abort.")
  invisible(gc())
  q()
}

str(ggp)

# phenotype
Phenotype <- opt$phenotype

# iteration
Iteration <- opt$iteration
Iter <- opt$iter
if(Iter == "none"){
  Iter <- NA
  outext <- 'all'
} else {
  outext <- Iter
}

# parallel
pc <- opt$pc
workers <- opt$workers

# quantile cutoffs
Quantiles <- list("SE" = opt$qse, "Mean" = opt$qmean)

# pmax detection
PMax <- opt$pmax

# pthreshold
Pthreshold <- opt$pthreshold

# output object
if(dir.exists(opt$outpath)){
  Outpath <- opt$outpath
} else {
  message("Directory of output path cannot be found. Creating...")
  dir.create(file.path(opt$outpath))
}
Outname <- opt$outname

verbose <- opt$verbose

message(paste0("Output File of run_gene_detection(): ", Outpath, "/", Outname, "_", outext,  ".rds"))

# CALCULATION--------------------------------------------------------------

coldata.sce <- SingleCellExperiment::colData(sce)
if(Iteration == "TestIterSet" && "Sampling" %in% colnames(SingleCellExperiment::colData(sce))){
  sce <- sce[, SingleCellExperiment::colData(sce)[, "Sampling"]]
}

print(sce)

DetectRes <- run_gene_detection(sce = sce,
                                ggpfit = ggp,
                                PhenotypeCol = Phenotype,
                                IterCol = Iteration,
                                Iter = Iter,
                                Quantile = Quantiles,
                                PThreshold = Pthreshold, 
                                PMax = PMax,
                                pc = pc,
                                workers = workers,
                                verbose = verbose)

# OUTPUT ------------------------------------------------------------------


saveRDS(DetectRes, file = paste0(Outpath, "/", Outname, "_", outext, ".rds"))

# FIN ---------------------------------------------------------------------

sessionInfo()

gc()

q()
