# OPTPARSE ----------------------------------------------------------------

#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse));

option_list = list(
  make_option(c("--prediction"), type="character", help="input rds prediction object with full path", metavar="character"),
  make_option(c("--metric"), type="character", help="Metric column name", metavar="character"),
  make_option(c("--exprdata"), type="character", help="species set of expression data", metavar="character"),
  make_option(c("--markerspecies"), type="character", help="species set of marker genes", metavar="character"),
  make_option(c("--classification"), type="character", help="expression data type (iteration or clone)", metavar="character"),
  make_option(c("--celltype"), type="character", help="over all cells or celltype", metavar="character"),
  make_option(c("--level"), type="character", help="consider celltype or germlayer", metavar="character"),
  make_option(c("--type"), type="character", help="unique, multiple or all per level", metavar="character"),
  make_option(c("--biotype"), type="character", help="gene biotype of markers", metavar="character"),
  make_option(c("--lowern"), type="numeric", help="minimal number of marker genes", default = 1, metavar="character"),
  make_option(c("--uppern"), type="numeric", help="maximal number of marker genes", default = 50, metavar="character"),
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
suppressPackageStartupMessages(require(stats))
suppressPackageStartupMessages(require(investr))


# INPUT -------------------------------------------------------------------

# read in option arguments

# prediction object
if(file.exists(opt$prediction)){
  prediction <- readRDS(file = opt$prediction)
} else {
  message("Input prediction object could not be read. Abort.")
  invisible(gc())
  q()
}

head(prediction)

if(opt$celltype == "overall"){
  input <- prediction %>% 
    dplyr::filter(exprdata == opt$exprdata &
                    markerspecies == opt$markerspecies &
                    classification == opt$classification &
                    biotype == opt$biotype &
                    ct.level == opt$level &
                    type == opt$type &
                    nogenes %in% c(opt$lowern:opt$uppern))
} else {
  input <- prediction %>% 
    dplyr::filter(Celltype == opt$celltype &
                    exprdata == opt$exprdata &
                    markerspecies == opt$markerspecies &
                    classification == opt$classification &
                    biotype == opt$biotype &
                    ct.level == opt$level &
                    type == opt$type &
                    nogenes %in% c(opt$lowern:opt$uppern))
}

head(input)

# output object
if(dir.exists(opt$outpath)){
  Outpath <- opt$outpath
} else {
  message("Directory of output path cannot be found. Creating...")
  dir.create(file.path(opt$outpath))
}
Outname <- opt$outname

message(paste0("Output File of summarise_predict(): ", 
               Outpath, "/", Outname, 
               ".rds"))



# STEP BY STEP ------------------------------------------------------------

# yval = opt$metric
# xval = "nogenes"
# 
# # SORT INPUT
# dat.srt <- stats::sortedXyData(input[, xval, drop = T], 
#                                input[, yval, drop = T], 
#                                input)
# 
# print(head(dat.srt))
# 
# dat.val <- tibble(x = input[, xval, drop = T], 
#                   y = input[, yval, drop = T]) %>% 
#   dplyr::group_by(x) %>% 
#   dplyr::summarise(ymean = mean(y, na.rm = T),
#                    yse = sd(y , na.rm = T)/sqrt(n()))
# 
# print(head(dat.val))
# 
# # RUN Assymptotic FIT
# print(nls(y ~ SSasymp(x, Asym, R0, lrc), data = dat.srt))
# unconstrained.ssasymp.fit <- nls(y ~ SSasymp(x, Asym, R0, lrc), data = dat.srt)
# 
# print(summary(unconstrained.ssasymp.fit))
# 
# 
# # ESTIMATE FIT VALUES
# cimat1 <- investr::predFit(unconstrained.ssasymp.fit, 
#                            interval = "confidence", 
#                            level = 0.95) %>% 
#   data.frame() %>% 
#   dplyr::rename(yhat = fit, lowerCI = lwr, upperCI = upr)
# 
# pimat1 <- investr::predFit(unconstrained.ssasymp.fit, 
#                            interval = "prediction",
#                            level = 0.95) %>% 
#   data.frame() %>% 
#   dplyr::rename(yhat = fit, lowerPI = lwr, upperPI = upr) %>% 
#   dplyr::select(lowerPI, upperPI)
# 
# asymptote.val <- stats::coef(unconstrained.ssasymp.fit)[1]
# r0.val <- stats::coef(unconstrained.ssasymp.fit)[2]
# start.val <- stats::fitted(unconstrained.ssasymp.fit)[1]
# half.val <- start.val + ((asymptote.val - start.val)/2)
# 
# nogenes.max <- ceiling(NLSstClosestX(dat.srt, asymptote.val))
# nogenes.half <- ceiling(NLSstClosestX(dat.srt, half.val))
# 
# unconstrained.ssasymp.fit$Fit <- dplyr::bind_cols(dat.val, 
#                                                   cimat1,
#                                                   pimat1)
# 
# print(head(unconstrained.ssasymp.fit$Fit))
# 
# # OUTPUT
# res <- data.frame(point = c("asymptote", "r0", 'r1', 'r0.5'),
#                   xval = c(nogenes.max, 0, 1, nogenes.half),
#                   yval = c(asymptote.val, r0.val, start.val, half.val),
#                   percyval = c(100, 
#                                (r0.val/asymptote.val) * 100,
#                                (start.val/asymptote.val) * 100,
#                                (half.val/asymptote.val) * 100)
# ) 
# 
# out <- list("Data" = res,
#             "Fit" = unconstrained.ssasymp.fit)
# 
# str(out)

# CALCULATION--------------------------------------------------------------

SummaryData <- summarise_predict(input = input,
                                 yval = opt$metric,
                                 xval = "nogenes")

str(SummaryData)

# OUTPUT ------------------------------------------------------------------


saveRDS(SummaryData, file = paste0(Outpath, "/", Outname,  ".rds"))

# FIN ---------------------------------------------------------------------

sessionInfo()

gc()

q()
