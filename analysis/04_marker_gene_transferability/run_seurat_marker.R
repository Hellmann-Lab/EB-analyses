# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse));

option_list = list(
  make_option(c("--input"), type="character", help="input rds sce object with full path", metavar="character"),
  make_option(c("--phenotype"), type="character", help="phenotype column name", metavar="character"),
  make_option(c("--iteration"), type="character", help="iteration set column name", metavar="character"),
  make_option(c("--iter"), type="numeric", help="iteration number, can also be NA then takes all", metavar="character"),
  make_option(c("--mincell"), type="numeric", help="minimum number of cells per group", default = 10, metavar="character"),
  make_option(c("--maxcell"), type="numeric", help="maximum number of cells per group", default = 250, metavar="character"),
  make_option(c("--pc"), type="logical", help="parallelisation or not", default = TRUE, metavar="character"),
  make_option(c("--workers"), type="numeric", help="number of workers for parallel computation", default = 2, metavar="character"),
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

suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(scuttle))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(future))


# INPUT -------------------------------------------------------------------

# read in option arguments

# sce object
if(file.exists(opt$input)){
  sce <- readRDS(file = opt$input)
} else {
  message("Input object could not be read. Abort.")
  invisible(gc())
  q()
}

print(sce)

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

# min cells
MinCells <- opt$mincell

# max cells 
MaxCells <- opt$maxcell

# output object
if(dir.exists(opt$outpath)){
  Outpath <- opt$outpath
} else {
  message("Directory of output path cannot be found. Creating...")
  dir.create(file.path(opt$outpath))
}
Outname <- opt$outname

message(paste0("Output File of run_Seurat(): ", Outpath, "/", Outname, "_", outext,  ".rds"))

# CALCULATION--------------------------------------------------------------

coldata.sce <- SingleCellExperiment::colData(sce)
if(Iteration == "TestIterSet" && "Sampling" %in% colnames(SingleCellExperiment::colData(sce))){
  sce <- sce[, SingleCellExperiment::colData(sce)[, "Sampling"]]
}

print(sce)

SeuratRes <- run_Seurat(sce, 
                        IterCol = Iteration, 
                        Iter = Iter,
                        PhenotypeCol = Phenotype, 
                        MinCells = MinCells,
                        MaxCells = MaxCells,
                        pc = pc,
                        workers = workers)

# OUTPUT ------------------------------------------------------------------


saveRDS(SeuratRes, file = paste0(Outpath, "/", Outname, "_", outext, ".rds"))

# FIN ---------------------------------------------------------------------

sessionInfo()

gc()

q()
