# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse));

option_list = list(
  make_option(c("--input"), type="character", help="input rds sce object with full path", metavar="character"),
  make_option(c("--phenotype"), type="character", help="phenotype column name", metavar="character"),
  make_option(c("--iteration"), type="character", help="iteration set column name", metavar="character"),
  make_option(c("--iter"), type="numeric", help="iteration number, can also be NA then takes all", metavar="character"),
  make_option(c("--tau"), type="character", help="taus sequence for umi or nonumi data", default = "umi", metavar="character"),
  make_option(c("--mtc"), type="character", help="mutliple testing correction method for p.adjust()", default = "BH", metavar="character"),
  make_option(c("--clean"), type="logical", help="clean output or not", default = TRUE, metavar="character"),
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
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(scuttle))
suppressPackageStartupMessages(require(quantreg))
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(fastDummies))


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

# other
clean <- opt$clean
tau <- opt$tau
MTC <- opt$mtc

# output object
if(dir.exists(opt$outpath)){
  Outpath <- opt$outpath
} else {
  message("Directory of output path cannot be found. Creating...")
  dir.create(file.path(opt$outpath))
}
Outname <- opt$outname

message(paste0("Output File of run_ZIQRank_Combination(): ", Outpath, "/", Outname, "_", outext,  ".rds"))

# CALCULATION--------------------------------------------------------------

coldata.sce <- SingleCellExperiment::colData(sce)
if(Iteration == "TestIterSet" && "Sampling" %in% colnames(SingleCellExperiment::colData(sce))){
  sce <- sce[, SingleCellExperiment::colData(sce)[, "Sampling"]]
}

# sce <- sce[sample(rownames(sce), 1000, replace = F), ]

print(sce)

ZIQRankRes <- run_ZIQRank_Combination(sce, 
                                      IterCol = Iteration, 
                                      Iter = Iter,
                                      PhenotypeCol = Phenotype, 
                                      tau=tau, 
                                      MTC = MTC,
                                      clean = clean,
                                      pc = pc,
                                      workers = workers)

# OUTPUT ------------------------------------------------------------------


saveRDS(ZIQRankRes, file = paste0(Outpath, "/", Outname, "_", outext, ".rds"))

# FIN ---------------------------------------------------------------------

sessionInfo()

gc()

q()
