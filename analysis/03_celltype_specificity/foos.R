
# powsimR prep ------------------------------------------------------------

create_estimate_inputs <- function(sce,
                                   Phenotype,
                                   phenotype.args,
                                   Aggregation,
                                   aggregation.args,
                                   Marker,
                                   marker.args,
                                   Batch,
                                   batch.args,
                                   Nuisance,
                                   nuisance.args,
                                   Noise,
                                   noise.args,
                                   Estimation,
                                   estimation.args,
                                   Selection,
                                   selection.args,
                                   parallel,
                                   BPPARAM,
                                   verbose,
                                   clean,
                                   infilepath,
                                   inname,
                                   outfilepath,
                                   outname,
                                   ncores,
                                   mem) {
  
  if(inname == outname){
    stop(message("Please choose unique input and output names"))
  }
  
  dir.create(file.path(infilepath, inname), recursive = T)
  
  # sce
  saveRDS(sce, file = paste0(infilepath, "/", inname, "/", inname, ".rds"))
  
  # phenotype
  if(missing(Phenotype)){
    Phenotype <- FALSE
  } else {
    if(!is.logical(Phenotype)){
      Phenotype <- FALSE
    }
  }
  if(!missing(phenotype.args)){
    saveRDS(phenotype.args, 
            file = paste0(infilepath, "/", inname, "/", "phenotype", ".rds"))
  } else {
    phenotype.args <- NULL
  }
  
  
  # aggregation
  if(missing(Aggregation)){
    Aggregation <- FALSE
  } else {
    if(!is.logical(Aggregation)){
      Aggregation <- FALSE
    }
  }
  if(!missing(aggregation.args)){
    saveRDS(aggregation.args, 
            file = paste0(infilepath, "/", inname, "/", "aggregation", ".rds"))
  } else {
    aggregation.args <- NULL
  }
  
  # marker
  if(missing(Marker)){
    Marker <- FALSE
  } else {
    if(!is.logical(Marker)){
      Marker <- FALSE
    }
  }
  if(!missing(marker.args)){
    saveRDS(marker.args, 
            file = paste0(infilepath, "/", inname, "/", "marker", ".rds"))
  } else {
    marker.args <- NULL
  }
  
  # batch
  if(missing(Batch)){
    Batch <- FALSE
  } else {
    if(!is.logical(Batch)){
      Batch <- FALSE
    }
  }
  if(!missing(batch.args)){
    saveRDS(batch.args, 
            file = paste0(infilepath, "/", inname, "/", "batch", ".rds"))
  } else {
    batch.args <- NULL
  }
  
  # nuisance
  if(missing(Nuisance)){
    Nuisance <- FALSE
  } else {
    if(!is.logical(Nuisance)){
      Nuisance <- FALSE
    }
  }
  if(!missing(nuisance.args)){
    saveRDS(nuisance.args, 
            file = paste0(infilepath, "/", inname, "/", "nuisance", ".rds"))
  } else {
    nuisance.args <- NULL
  }
  
  # noise
  if(missing(Noise)){
    Noise <- FALSE
  } else {
    if(!is.logical(Noise)){
      Noise <- FALSE
    }
  }
  if(!missing(noise.args)){
    saveRDS(noise.args, 
            file = paste0(infilepath, "/", inname, "/", "noise", ".rds"))
  } else {
    noise.args <- NULL
  }
  
  # estimation
  if(missing(Estimation)){
    Estimation <- TRUE
  } else {
    if(!is.logical(Estimation)){
      Estimation <- TRUE
    }
  }
  if(!missing(estimation.args)){
    saveRDS(estimation.args, 
            file = paste0(infilepath, "/", inname, "/", "estimation", ".rds"))
  } else {
    estimation.args <- NULL
  }
  
  # selection
  if(missing(Selection)){
    Selection <- FALSE
  } else {
    if(!is.logical(Selection)){
      Selection <- FALSE
    }
  }
  if(!missing(selection.args)){
    saveRDS(selection.args, 
            file = paste0(infilepath, "/", inname, "/", "selection", ".rds"))
  } else {
    selection.args <- NULL
  }
  
  
  # parallel
  if(missing(parallel)){
    parallel <- FALSE
  } else {
    if(!is.logical(parallel)){
      parallel <- FALSE
    }
  }
  if(!missing(BPPARAM)){
    saveRDS(BPPARAM, 
            file = paste0(infilepath, "/", inname, "/", "BPPARAM", ".rds"))
  } else {
    BPPARAM <- NULL
  }
  
  if(missing(verbose)){
    verbose <- TRUE
  } else {
    if(!is.logical(verbose)){
      verbose <- TRUE
    }
  }
  
  if(missing(clean)){
    clean <- TRUE
  } else {
    if(!is.logical(clean)){
      clean <- TRUE
    }
  }
  
  dir.create(file.path(outfilepath, outname), recursive = TRUE)
  
  # parameters text file
  param.dat <- data.frame(sce =  paste0(infilepath, "/", inname, "/", inname, ".rds"),
                          Phenotype = Phenotype,
                          phenotype.args = ifelse(is.null(phenotype.args), 
                                                  NA,
                                                  paste0(infilepath, "/", inname, "/", "phenotype", ".rds")),
                          Aggregation = Aggregation,
                          aggregation.args = ifelse(is.null(aggregation.args), 
                                                    NA,
                                                    paste0(infilepath, "/", inname, "/", "aggregation", ".rds")),
                          Marker = Marker,
                          marker.args = ifelse(is.null(marker.args), 
                                               NA,
                                               paste0(infilepath, "/", inname, "/", "marker", ".rds")),
                          Batch = Batch,
                          batch.args = ifelse(is.null(batch.args), 
                                              NA,
                                              paste0(infilepath, "/", inname, "/", "batch", ".rds")),
                          Nuisance = Nuisance,
                          nuisance.args = ifelse(is.null(nuisance.args), 
                                                 NA,
                                                 paste0(infilepath, "/", inname, "/", "nuisance", ".rds")),
                          Noise = Noise,
                          noise.args = ifelse(is.null(noise.args), 
                                              NA,
                                              paste0(infilepath, "/", inname, "/", "noise", ".rds")), 
                          Estimation = Estimation,
                          estimation.args = ifelse(is.null(estimation.args), 
                                                   NA,
                                                   paste0(infilepath, "/", inname, "/", "estimation", ".rds")), 
                          Selection = Selection,
                          selection.args = ifelse(is.null(selection.args), 
                                                  NA,
                                                  paste0(infilepath, "/", inname, "/", "selection", ".rds")), 
                          parallel = parallel,
                          BPPARAM = ifelse(is.null(BPPARAM), 
                                           NA,
                                           paste0(infilepath, "/", inname, "/", "BPPARAM", ".rds")),  
                          verbose = verbose,
                          clean = clean,
                          outfilepath = paste0(outfilepath, "/", outname),
                          outname = outname,
                          ncores = ncores,
                          mem = mem)
  
  write.table(param.dat, 
              file = paste0(infilepath, "/", inname, "/", "param", ".txt"),
              sep = "\t", row.names = T, quote = F, col.names = F, na='none')
  
}

# sampling and iteration split --------------------------------------------

run_downsample_iter <- function(sce,
                                SpeciesCol = NULL,
                                PhenotypeCol,
                                maxc = 5000,
                                iter = 3,
                                verbose = TRUE) {
  
  if(!is.null(SpeciesCol)){
    sce <- subset(sce, , species == SpeciesCol)
  }
  
  # annotation information
  celltype <- SingleCellExperiment::colData(sce)[, PhenotypeCol]
  
  totalcells <- ncol(sce)
  totalpop <- maxc * iter
  if(totalpop > totalcells){
    maxc <- floor(totalcells / iter)
  }
  
  BC <- colnames(sce)
  
  samp_no <- tibble(phenotype = celltype,
                    stringsAsFactors = F) %>%
    dplyr::group_by(phenotype) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::mutate(new_n = round(proportion * maxc)) %>%
    dplyr::mutate(new_n = ifelse(proportion < 0.01, n, new_n)) %>% 
    dplyr::mutate(total_n = new_n * iter) %>% 
    dplyr::mutate(total_n = ifelse(n < total_n, new_n, total_n))
  
  
  BC_samp <- data.frame(BC = BC,
                        phenotype = celltype,
                        stringsAsFactors = F) %>%
    dplyr::group_by(phenotype) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::left_join(samp_no, by = c('phenotype' = "phenotype")) %>%
    dplyr::mutate(samp = purrr::map2(data, total_n, dplyr::sample_n)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(samp) %>%
    dplyr::pull(BC)
  
  Sampling <- BC %in% BC_samp
  SummarizedExperiment::colData(sce)[, "Sampling"] <- Sampling
  
  sce_samp <- sce[, colnames(sce) %in% BC_samp]
  samp.dat <- data.frame(BC_id = colnames(sce_samp), 
                         phenotype = SingleCellExperiment::colData(sce_samp)[, PhenotypeCol],
                         stringsAsFactors = F) 
  samp.L <- split(x = samp.dat, f = samp.dat$phenotype) 
  samp.L <- lapply(1:length(samp.L), function(i){
    samp.K <- rep_len(x = 1:iter, length.out = nrow(samp.L[[i]]))
    samp.L[[i]][, 'K'] <- sample(samp.K, size = length(samp.K), replace = F)
    samp.L[[i]][, c("BC_id", "K")]
  })
  samp.dat <- do.call('rbind', samp.L)
  all.dat <- data.frame(BC_id = BC) %>% 
    dplyr::left_join(samp.dat, by = "BC_id")
  
  SummarizedExperiment::colData(sce)[, "TestIterSet"] <- all.dat$K
  
  invisible(gc())
  return(sce)
}

# glmgampoi ---------------------------------------------------------------

run_scran_clust <- function(sce,
                            PhenotypeCol,
                            BPPARAM,
                            verbose) {
  
  if(verbose){
    message(paste0("Normalisation: Using scran with pooling over clusters!"))
  }
  
  # identify clusters and define pool sizes accordingly
  ncelltypes <- min(table(SingleCellExperiment::colData(sce)[, PhenotypeCol]))
  
  min.size <- floor(ncelltypes*0.75)
  
  clusters <- scran::quickCluster(sce,
                                  min.size = min.size,
                                  BPPARAM = BPPARAM)
  
  if(ncol(sce)<=14) {
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, by = 1))))
  }
  if(ncol(sce)>14 & ncol(sce)<=50) {
    sizes <- unique(c(round(seq(from=2, to=min(table(clusters))-1, length.out=10))))
  }
  if(ncol(sce)>50 & ncol(sce)<=5000) {
    sizes <- unique(c(round(seq(from=10, to=min(table(clusters))-1, length.out=20))))
  }
  if(trunc(ncol(sce))>5000) {
    sizes <- unique(c(round(seq(from=20, to=min(table(clusters))-1, length.out=20))))
  }
  
  
  # calculate size factors
  sf <- suppressWarnings(scuttle::pooledSizeFactors(x = sce, sizes = sizes,
                                                    clusters = clusters, ref.clust = NULL, max.cluster.size = 3000,
                                                    positive = TRUE, scaling = NULL, min.mean = NULL,
                                                    subset.row = NULL, BPPARAM = BPPARAM))
  names(sf) <- colnames(sce)
  
  invisible(gc())
  
  return(sf)
}


run_glm_gampoi <- function(sce, 
                          IterCol = "TestIterSet",
                          Iter = 1,
                          PhenotypeCol,
                          pc = FALSE,
                          workers = NULL,
                          verbose = TRUE){
  
  if(!is.na(Iter) && !is.null(IterCol)){
    TestIterSet <- SingleCellExperiment::colData(sce)[, IterCol]
    sce <- SummarizedExperiment::subset(sce, , TestIterSet == Iter)
  } 
  
  celltypes <- SingleCellExperiment::colData(sce)[, PhenotypeCol]
  modelmatrix <- stats::model.matrix(~ 0 + celltypes)
  
  # setting parameter on_disk of glm_gp() depending on memory consumption
  mat.mem <- format(object.size(SummarizedExperiment::assay(sce, 'counts')), units = "GB", digits = 2)
  mat.mem <- as.numeric(sapply(strsplit(mat.mem, " "), "[[", 1))
  pc.mem <- memuse::Sys.meminfo()$freeram@size
  on_disk <- ifelse(mat.mem*100 > pc.mem, TRUE, FALSE) # since format() with as.matrix() call would give a warning message and might actually break on smaller ram systems, i removed as.matrix() and multiplied by 100 instead of 1.25 so that there should be enough memory if glmgampoi internally does some conversion to dense matrix.
  
  # glm gam poi fit
  if(verbose){
    message(paste0('Fit a Gamma-Poisson Generalized Linear Model ...'))
  }
  if(verbose){
    message(paste0('... full model'))
  }
  fit.mod <- glmGamPoi::glm_gp(data = sce,
                           design = modelmatrix,
                           col_data = NULL,
                           reference_level = NULL,
                           offset = 0,
                           size_factors = SingleCellExperiment::sizeFactors(sce),
                           overdispersion = TRUE,
                           overdispersion_shrinkage = TRUE,
                           do_cox_reid_adjustment = FALSE, 
                           subsample = FALSE,
                           on_disk = on_disk,
                           verbose = verbose)
  if(verbose){
    message(paste0('... intercept model'))
  }
  fit.null <- glmGamPoi::glm_gp(data = sce,
                               design = ~1,
                               col_data = NULL,
                               reference_level = NULL,
                               offset = 0,
                               size_factors = SingleCellExperiment::sizeFactors(sce),
                               overdispersion = TRUE,
                               overdispersion_shrinkage = TRUE,
                               do_cox_reid_adjustment = FALSE, 
                               subsample = FALSE,
                               on_disk = on_disk,
                               verbose = verbose)
  
  # predictions
  if(verbose){
    message(paste0('Predict Expression using Gamma-Poisson Generalized Linear Model ...'))
  }
  if(verbose){
    message(paste0('... full model'))
  }
  pred.mod.response <- glmGamPoi:::predict.glmGamPoi(object = fit.mod,
                                            type = "response",
                                            se.fit = TRUE,
                                            on_disk = on_disk,
                                            verbose = verbose)
  pred.mod.link <- glmGamPoi:::predict.glmGamPoi(object = fit.mod,
                                                 type = "link",
                                                 se.fit = TRUE,
                                                 on_disk = on_disk,
                                                 verbose = verbose)
  if(verbose){
    message(paste0('... intercept model'))
  }
  pred.null.response <- glmGamPoi:::predict.glmGamPoi(object = fit.null,
                                                     type = "response",
                                                     se.fit = TRUE,
                                                     on_disk = on_disk,
                                                     verbose = verbose)
  pred.null.link <- glmGamPoi:::predict.glmGamPoi(object = fit.null,
                                                 type = "link",
                                                 se.fit = TRUE,
                                                 on_disk = on_disk,
                                                 verbose = verbose)
  
  
  # residuals
  if(verbose){
    message(paste0('Compute randomized quantile residuals of Gamma-Poisson Generalized Linear Model ...'))
  }
  if(verbose){
    message(paste0('... full model'))
  }
  ResidualData.mod <- transformGamPoi::residual_transform(data = fit.mod,
                                                      residual_type = 'randomized_quantile',
                                                      overdispersion = TRUE,
                                                      overdispersion_shrinkage = TRUE,
                                                      on_disk = on_disk,
                                                      verbose = verbose)
  ResidualData.mod <- as(ResidualData.mod, "TsparseMatrix")
  
  
  if(verbose){
    message(paste0('... intercept model'))
  }
  ResidualData.null <- transformGamPoi::residual_transform(data = fit.null,
                                                      residual_type = 'randomized_quantile',
                                                      overdispersion = TRUE,
                                                      overdispersion_shrinkage = TRUE,
                                                      on_disk = on_disk,
                                                      verbose = verbose)
  ResidualData.null <- as(ResidualData.null, "TsparseMatrix")
  
  out <- list("Fit" = list("Full" = fit.mod, 
                           "Null" = fit.null),
              "Residuals" = list("Full" = ResidualData.mod, 
                                 "Null" = ResidualData.null),
              "Prediction" = list("Full" = list("Response" = pred.mod.response,
                                                "Link" = pred.mod.link), 
                                  "Null" = list("Response" = pred.null.response,
                                                "Link" = pred.null.link))
                )
  
  invisible(gc())
  
  return(out)
  
  
}


# gene detection ----------------------------------------------------------


run_gene_detection <- function(sce,
                               ggpfit,
                               PhenotypeCol = "manual_annotation",
                               IterCol = "TestIterSet",
                               Iter = 1,
                               Quantile = list("SE" = 0.5, "Mean" = 0.1),
                               PMax = 0.5,
                               PThreshold = 0.9999, 
                               pc = TRUE,
                               workers = 2,
                               verbose = TRUE){
  
  if(isTRUE(pc)){
    BPPARAM = BiocParallel::MulticoreParam(workers = workers, tasks = 0L,
                                           stop.on.error = TRUE,
                                           progressbar = FALSE, RNGseed = NULL,
                                           timeout = 30L * 24L * 60L * 60L, exportglobals=TRUE,
                                           log = FALSE, threshold = "INFO", logdir = NA_character_,
                                           resultdir = NA_character_, jobname = "BPJOB",
                                           manager.hostname = NA_character_, manager.port = NA_integer_)
  } else {
    BPPARAM = BiocParallel::SerialParam()
  }
  
  # GOAL
  # find genes that are considered expressed above background
  # use glm gampoi fitting results and estimated gene dropout in a logistic regression model
  # find cut point for dropout where genes do not include zero in confidence interval anymore
  # do this once for specific glm with model matrix (celltype specific) and once without celltype annotation ('grand')
  # also record beta coefficient and s.e. of beta coefficient in output
  
  # subset sce to iteration specified
  if(!is.na(Iter) && !is.null(IterCol)){
    TestIterSet <- SingleCellExperiment::colData(sce)[, IterCol]
    sce <- SummarizedExperiment::subset(sce, , TestIterSet == Iter)
  } 
  
  phenotype <- SingleCellExperiment::colData(sce)[, PhenotypeCol]
  phenotypes <- unique(phenotype)
  
  # grand
  # mean expression and standard error
  l.ggp.mu <- ggpfit$Prediction$Null$Link$fit
  r.ggp.mu <- ggpfit$Prediction$Null$Response$fit
  l.ggp.se <- ggpfit$Prediction$Null$Link$se.fit
  r.ggp.se <- ggpfit$Prediction$Null$Response$se.fit
  l.gmedian.grand <- matrixStats::rowMedians(l.ggp.mu, na.rm = TRUE)
  r.gmedian.grand <- matrixStats::rowMedians(r.ggp.mu, na.rm = TRUE)
  l.gse.grand <- matrixStats::rowMedians(l.ggp.se, na.rm = TRUE)
  r.gse.grand <- matrixStats::rowMedians(r.ggp.se, na.rm = TRUE)
  # dropout
  mat_0 <- SingleCellExperiment::counts(sce) == 0
  gdrop.grand <- (ncol(mat_0) - Matrix::rowSums(!mat_0))/ncol(mat_0)
  # data for fitting
  FitData <- tibble(Gene = names(l.gmedian.grand),
                    LMean = l.gmedian.grand,
                    LSE = l.gse.grand,
                    Mean = r.gmedian.grand,
                    SE = r.gse.grand,
                    Detection = 1-gdrop.grand) %>% 
    dplyr::mutate(SESize = SE/Mean)
  QSESize <- 2^quantile(log2(FitData$SESize[is.finite(FitData$SESize)]), Quantile$SE, na.rm = T)
  QLMean <- quantile(FitData$LMean[is.finite(FitData$LMean)], Quantile$Mean, na.rm = T)
  QMean <- 2^quantile(FitData$LMean[is.finite(FitData$LMean)], Quantile$Mean, na.rm = T)
  FitData.Grand <- FitData %>% 
    dplyr::mutate(Estimable = ifelse(Mean<QMean, 0, 1)) %>% 
    dplyr::mutate(Reliable = ifelse(SESize>QSESize, 0, 1)) %>% 
    dplyr::mutate(Reliable = ifelse(is.na(Reliable), 0, Reliable)) %>% 
    dplyr::mutate(Expressed = case_when(Estimable ==1 & Reliable == 1 ~ 1,
                                        TRUE ~ 0))
  
  # per cell type
  FitData.CT <- BiocParallel::bplapply(1:length(phenotypes), function(p){
    pheno <- phenotypes[p]
    sce.i <- SummarizedExperiment::subset(sce, , phenotype == pheno)
    bc <- colnames(sce.i)
    if(length(bc)>=10){
      # mean expression and standard error
      l.ggp.mu <- ggpfit$Prediction$Full$Link$fit[,bc]
      r.ggp.mu <- ggpfit$Prediction$Full$Response$fit[,bc]
      l.ggp.se <- ggpfit$Prediction$Full$Link$se.fit[,bc]
      r.ggp.se <- ggpfit$Prediction$Full$Response$se.fit[,bc]
      l.gmedian <- matrixStats::rowMedians(l.ggp.mu, na.rm = TRUE)
      r.gmedian <- matrixStats::rowMedians(r.ggp.mu, na.rm = TRUE)
      l.gse <- matrixStats::rowMedians(l.ggp.se, na.rm = TRUE)
      r.gse <- matrixStats::rowMedians(r.ggp.se, na.rm = TRUE)
      # dropout
      mat_0 <- SingleCellExperiment::counts(sce.i) == 0
      gdrop <- (ncol(mat_0) - Matrix::rowSums(!mat_0))/ncol(mat_0)
      # data for fitting
      FitData <- tibble(Gene = names(l.gmedian),
                        LMean = l.gmedian,
                        LSE = l.gse,
                        Mean = r.gmedian,
                        SE = r.gse,
                        Detection = 1-gdrop) %>% 
        dplyr::mutate(SESize = SE/Mean)
      QSESize <- 2^quantile(log2(FitData$SESize[is.finite(FitData$SESize)]), Quantile$SE, na.rm = T)
      QMean <- 2^quantile(FitData$LMean[is.finite(FitData$LMean)], Quantile$Mean, na.rm = T)
      FitData <- FitData %>% 
        dplyr::mutate(Estimable = ifelse(Mean<QMean, 0, 1)) %>% 
        dplyr::mutate(Reliable = ifelse(SESize>QSESize, 0, 1)) %>% 
        dplyr::mutate(Reliable = ifelse(is.na(Reliable), 0, Reliable)) %>% 
        dplyr::mutate(Expressed = case_when(Estimable ==1 & Reliable == 1 ~ 1,
                                            TRUE ~ 0))
    } else {
      FitData <- NULL
    }
    
    return(FitData)
  }, BPPARAM = BPPARAM)
  names(FitData.CT) <- phenotypes
  
  FitData.CT <- FitData.CT[!sapply(FitData.CT, is.null)]
  
  # logistic regression
  # grand
  FitResults.Grand <- fit_logreg(FitData = FitData.Grand, 
                                 PX = PThreshold,
                                 PMax = PMax)
  
  # per cell type
  FitResults.CT <- sapply(names(FitData.CT), function(i){
    fit_logreg(FitData = FitData.CT[[i]],
               PX = PThreshold,
               PMax = PMax)
  }, USE.NAMES = T, simplify = F)
  
  Out <- list("Grand" = FitResults.Grand,
              "CT" = FitResults.CT)
  
  invisible(gc())
  
  return(Out)

}

# fit logistic regression with Firth penalization
fit_logreg <- function(FitData, PX = 0.9999, PMax = 0.5){
  
  FitData.Red <- FitData %>% dplyr::filter(Detection < PMax)
  
  # Firth’s modified score procedure in logistic regression analysis
  logreg.model <- try(logistf::logistf(Expressed ~ Detection, data=FitData.Red), 
                      silent = TRUE)
  
  if(inherits(logreg.model, "try-error")){
    message("Logistic regression model fit failed. Setting default cutoff 0.25 for detection.")
    Cutoff <- 0.25
    FitData <- FitData %>% 
      dplyr::mutate(Fitted = ifelse(Detection >= Cutoff, 1, 0))
    res <- list("Data" = FitData,
                "Cutoff" = Cutoff)
  } else{
    Cutoff <- probX(p = PX, model = logreg.model)
    FitData <- FitData %>% 
      dplyr::mutate(Fitted = ifelse(Detection >= Cutoff, 1, 0))
    res <- list("Data" = FitData,
                "Cutoff" = Cutoff)
  }
  return(res)
}

# Function to get Detection value for a given Expression probability assuming normal
probX = function(p, model) {
  as.numeric((qnorm(p) - coef(model)[1])/coef(model)[2])
}


process_gene_detection <- function(DetectResList,
                                   Type = c("Subsampling", "Clone"),
                                   Detection = list("All" = c("Subsampling" = 'any', "Clone" = 'all'),
                                                    "CT" = c("Subsampling" = 'any', "Clone" = 'all')),
                                   Beta_Interval = list("All" = c("Subsampling" = 'within', "Clone" = 'any'),
                                                       "CT" = c("Subsampling" = 'within', "Clone" = 'any'))){
  
  if(Type == "Subsampling"){

    # over all cells and subsamplings
    Data.All <- DetectResList$All$Grand$Data %>% 
      dplyr::rename(Expressed_g = Fitted) %>% 
      dplyr::select(Gene, LMean, LSE, Mean, SE, Detection, Expressed_g) %>% 
      dplyr::mutate(LMean_lower = LMean - 1.96 * LSE,
                    LMean_upper = LMean + 1.96 * LSE,
                    .after = SE)
    # detection
    # over all cells per subsampling
    Data.Iter.All.L <- lapply(names(DetectResList$Iteration), function(i){
      DetectResList$Iteration[[i]]$Grand$Data
    })
    # number of iterations
    n_iter <- length(Data.Iter.All.L)
    # detection
    detect.type <- Detection$All[Type]
    Detect.Data.Iter.All <- data.table::rbindlist(Data.Iter.All.L, idcol = "Iteration") %>% 
      tibble::as_tibble() %>% 
      tidyr::pivot_wider(id_cols = Gene,
                         names_from = Iteration, 
                         values_from = Fitted, 
                         names_prefix="K", 
                         values_fill  = 0) %>% 
      dplyr::mutate(Expressed_k = rowSums(across(where(is.numeric)))/n_iter) %>% 
      dplyr::select(Gene, Expressed_k)
    if(detect.type == "all"){
      Detect.Data.Iter.All <- Detect.Data.Iter.All %>% 
        dplyr::mutate(Expressed_k = ifelse(Expressed_k == 1, 1, 0))
    }
    
    # effect size
    iv.type <- Beta_Interval$All[Type]
    ES.Data.Iter.All <- data.table::rbindlist(Data.Iter.All.L, idcol = "Iteration") %>% 
      tibble::as_tibble() %>% 
      dplyr::select(Iteration, Gene, LMean, LSE, Mean, SE) %>% 
      dplyr::mutate(LMean_lower = LMean - 1.96 * LSE,
                    LMean_upper = LMean + 1.96 * LSE,
                    .after = SE) %>% 
      dplyr::left_join(Data.All %>% 
                         dplyr::select(Gene, LMean_lower, LMean_upper), 
                       by = "Gene",
                       suffix = c(".i", ".g")) %>% 
      mutate(Fitted = ivs::iv_overlaps(ivs::iv(LMean_lower.i, LMean_upper.i), 
                                       ivs::iv(LMean_lower.g, LMean_upper.g),
                                       missing = FALSE,
                                       type = iv.type)) %>% 
      tidyr::pivot_wider(id_cols = Gene,
                         names_from = Iteration, 
                         values_from = Fitted, 
                         names_prefix="K",
                         values_fill  = FALSE) %>% 
      dplyr::mutate(Reproducible_Beta_k = rowSums(across(where(is.logical)))/n_iter) %>% 
      dplyr::select(Gene, Reproducible_Beta_k)
    
    Data.Grand <- Data.All %>% 
      dplyr::left_join(Detect.Data.Iter.All, by = "Gene") %>% 
      dplyr::left_join(ES.Data.Iter.All, by = "Gene")
    
    # cell-type specific
    CT <- names(DetectResList$All$CT)
    Data.CT.L <- sapply(CT, function(ct){
      # over all cells of the celltype
      Data.All <- DetectResList$All$CT[[ct]]$Data %>% 
        dplyr::rename(Expressed_g = Fitted) %>% 
        dplyr::select(Gene, LMean, LSE, Mean, SE, Detection, Expressed_g) %>% 
        dplyr::mutate(LMean_lower = LMean - 1.96 * LSE,
                      LMean_upper = LMean + 1.96 * LSE,
                      .after = SE)
      # detection
      detect.type <- Detection$CT[Type]
      # iterations
      Data.Iter.CT.L <- lapply(names(DetectResList$Iteration), function(i){
        DetectResList$Iteration[[i]]$CT[[ct]]$Data
      })
      # number of iterations
      n_iter <- length(Data.Iter.CT.L)
      # detection
      Detect.Data.Iter.CT <- data.table::rbindlist(Data.Iter.CT.L, idcol = "Iteration") %>% 
        tibble::as_tibble() %>% 
        tidyr::pivot_wider(id_cols = Gene, 
                           names_from = Iteration, 
                           values_from = Fitted, 
                           names_prefix="K", 
                           values_fill  = 0) %>% 
        dplyr::mutate(Expressed_k = rowSums(across(where(is.numeric)))/n_iter) %>% 
        dplyr::select(Gene, Expressed_k)
      if(detect.type == "all"){
        Detect.Data.Iter.CT <- Detect.Data.Iter.CT %>% 
          dplyr::mutate(Expressed_k = ifelse(Expressed_k == 1, 1, 0))
      }
      
      # effect size
      iv.type <- Beta_Interval$CT[Type]
      ES.Data.Iter.CT <- data.table::rbindlist(Data.Iter.CT.L, idcol = "Iteration") %>% 
        tibble::as_tibble() %>% 
        dplyr::select(Iteration, Gene, LMean, LSE, Mean, SE) %>% 
        dplyr::mutate(LMean_lower = LMean - 1.96 * LSE,
                      LMean_upper = LMean + 1.96 * LSE,
                      .after = SE) %>% 
        dplyr::left_join(Data.All %>%  
                           dplyr::select(Gene, LMean_lower, LMean_upper), 
                         by = "Gene", 
                         suffix = c(".i", ".g")) %>% 
        mutate(Fitted = ivs::iv_overlaps(ivs::iv(LMean_lower.i, LMean_upper.i), 
                                         ivs::iv(LMean_lower.g, LMean_upper.g),
                                         missing = FALSE,
                                         type = iv.type)) %>% 
        tidyr::pivot_wider(id_cols = Gene,
                           names_from = Iteration, 
                           values_from = Fitted, 
                           names_prefix="K",
                           values_fill  = FALSE) %>% 
        dplyr::mutate(Reproducible_Beta_k = rowSums(across(where(is.logical)))/n_iter) %>% 
        dplyr::select(Gene, Reproducible_Beta_k)
      
      Data.CT <- Data.All %>% 
        dplyr::left_join(Detect.Data.Iter.CT, by = "Gene") %>% 
        dplyr::left_join(ES.Data.Iter.CT, by = "Gene")
      
      return(Data.CT)
      
    }, USE.NAMES = T, simplify = F)
    
    Data.CT <- data.table::rbindlist(Data.CT.L, idcol = "Celltype") %>% 
      tibble::as_tibble()

  }

  if(Type == "Clone"){
    clones <- names(DetectResList$Clone)
    # number of clones
    nclones <- length(clones)
    
    # over all cells and clones
    Data.All <- DetectResList$All$Grand$Data %>% 
      dplyr::rename(Expressed_g = Fitted) %>% 
      dplyr::select(Gene, LMean, LSE, Mean, SE, Detection, Expressed_g) %>% 
      dplyr::mutate(LMean_lower = LMean - LSE,
                    LMean_upper = LMean + LSE,
                    .after = SE)
    # over all cells per clone
    Data.Clone.All.L <- sapply(clones, function(clone){
      DetectResList$Clone[[clone]]$Grand$Data
    }, USE.NAMES = T, simplify = F)

    # detection
    detect.type <- Detection$All[Type]
    Detect.Data.Clone.All <- data.table::rbindlist(Data.Clone.All.L, idcol = "Clone") %>% 
      tibble::as_tibble() %>% 
      tidyr::pivot_wider(id_cols = Gene, 
                         names_from = Clone, 
                         values_from = Fitted, 
                         names_prefix="", 
                         values_fill  = 0) %>% 
      dplyr::mutate(Expressed_c = rowSums(across(where(is.numeric)))/nclones) %>% 
      dplyr::select(Gene, Expressed_c)
    
    if(detect.type == "all"){
      Detect.Data.Clone.All <- Detect.Data.Clone.All %>% 
        dplyr::mutate(Expressed_c = ifelse(Expressed_c == 1, 1, 0))
    }
    
    # effect size
    iv.type <- Beta_Interval$All[Type]
    ES.Data.Clone.All <- data.table::rbindlist(Data.Clone.All.L, idcol = "Clone") %>% 
      tibble::as_tibble() %>% 
      dplyr::select(Clone, Gene, LMean, LSE, Mean, SE) %>% 
      dplyr::mutate(LMean_lower = LMean - LSE,
                    LMean_upper = LMean + LSE,
                    .after = SE) %>% 
      dplyr::left_join(Data.All %>%  
                         dplyr::select(Gene, LMean_lower, LMean_upper), 
                       by = "Gene", suffix = c(".i", ".g")) %>% 
      mutate(Fitted = ivs::iv_overlaps(ivs::iv(LMean_lower.i, LMean_upper.i), 
                                       ivs::iv(LMean_lower.g, LMean_upper.g),
                                       missing = FALSE,
                                       type = iv.type)) %>% 
      tidyr::pivot_wider(id_cols = Gene,
                         names_from = Clone, 
                         values_from = Fitted, 
                         names_prefix="",
                         values_fill  = FALSE) %>% 
      dplyr::mutate(Reproducible_Beta_c = rowSums(across(where(is.logical)))/nclones) %>% 
      dplyr::select(Gene, Reproducible_Beta_c)
      
    
    Data.Grand <- Data.All %>% 
      dplyr::left_join(Detect.Data.Clone.All, by = "Gene") %>% 
      dplyr::left_join(ES.Data.Clone.All, by = "Gene")
    
    # celltype specific
    CT <- names(DetectResList$All$CT)
    Data.CT.L <- sapply(CT, function(ct){
      # all cells of celltype
      Data.All <- DetectResList$All$CT[[ct]]$Data %>% 
        dplyr::rename(Expressed_g = Fitted) %>% 
        dplyr::select(Gene, LMean, LSE, Mean, SE, Detection, Expressed_g) %>% 
        dplyr::mutate(LMean_lower = LMean - 1.96 * LSE,
                      LMean_upper = LMean + 1.96 * LSE,
                      .after = SE)

      # celltype specific per clone
      Data.CT.Clone.L <- sapply(clones, function(clone){
        DetectResList$Clone[[clone]]$CT[[ct]]$Data
      }, USE.NAMES = T, simplify = F)
      
      # detection
      detect.type <- Detection$CT[Type]
      Detect.Data.CT.Clone <- data.table::rbindlist(Data.CT.Clone.L, 
                                                    idcol = "Clone") %>% 
        tibble::as_tibble() %>% 
        tidyr::pivot_wider(id_cols = Gene, 
                           names_from = Clone, 
                           values_from = Fitted, 
                           names_prefix="", 
                           values_fill  = 0) %>% 
        dplyr::mutate(Expressed_c = rowSums(across(where(is.numeric)))/nclones) %>% 
        dplyr::select(Gene, Expressed_c)
      if(detect.type == "all"){
        Detect.Data.CT.Clone <- Detect.Data.CT.Clone %>% 
          dplyr::mutate(Expressed_c = ifelse(Expressed_c == 1, 1, 0))
      }
      
      # effect size
      iv.type <- Beta_Interval$CT[Type]
      ES.Data.CT.Clone <- data.table::rbindlist(Data.CT.Clone.L, idcol = "Clone") %>% 
        tibble::as_tibble() %>% 
        dplyr::select(Clone, Gene, LMean, LSE, Mean, SE) %>% 
        dplyr::mutate(LMean_lower = LMean - 1.96 * LSE,
                      LMean_upper = LMean + 1.96 * LSE,
                      .after = SE) %>% 
        dplyr::left_join(Data.All %>%  
                           dplyr::select(Gene, LMean_lower, LMean_upper), 
                         by = "Gene", 
                         suffix = c(".i", ".g")) %>% 
        mutate(Fitted = ivs::iv_overlaps(ivs::iv(LMean_lower.i, LMean_upper.i), 
                                         ivs::iv(LMean_lower.g, LMean_upper.g),
                                         missing = FALSE,
                                         type = iv.type)) %>% 
        tidyr::pivot_wider(id_cols = Gene,
                           names_from = Clone, 
                           values_from = Fitted, 
                           names_prefix="",
                           values_fill  = FALSE) %>% 
        dplyr::mutate(Reproducible_Beta_c = rowSums(across(where(is.logical)))/nclones) %>% 
        dplyr::select(Gene, Reproducible_Beta_c)
      
      Data.All <- Data.All %>% 
        dplyr::left_join(Detect.Data.CT.Clone, by = "Gene") %>% 
        dplyr::left_join(ES.Data.CT.Clone, by = "Gene")
      
      return(Data.All)
      
    }, USE.NAMES = T, simplify = F)
    
    Data.CT <- data.table::rbindlist(Data.CT.L, idcol = "Celltype") %>% 
      tibble::as_tibble()
    
  }
  
  # output
  Out <- list("Grand" = Data.Grand,
              "Celltype" = Data.CT,
              "Settings" = list("Type" = Type,
                              "Args" = list("Detection" = sapply(names(Detection), function(i){ Detection[[i]][Type] }, USE.NAMES = T, simplify = F),
                                            "Beta_Interval" = sapply(names(Beta_Interval), function(i){ Beta_Interval[[i]][Type] }, USE.NAMES = T, simplify = F))))
  
  return(Out)  
}


process_cutoff <- function(GeneDetectRes){
  # process cutoff res per cell type
  CutoffRes.CT.L <- sapply(names(GeneDetectRes$CT), function(ct){
    GeneDetectRes$CT[[ct]]$Cutoff
  }, USE.NAMES = T, simplify = F)
  CutoffRes.CT <- data.frame(do.call('rbind', CutoffRes.CT.L))
  colnames(CutoffRes.CT) <- "min.spct"
  CutoffRes.CT[,"celltype"] <- rownames(CutoffRes.CT)
  CutoffRes.CT <- tibble::as_tibble(CutoffRes.CT)
  #process cutoff res over all cells
  CutoffRes.Grand <- GeneDetectRes$Grand$Cutoff
  
  Out <- list("Grand" = CutoffRes.Grand, "Celltype" = CutoffRes.CT)
}

# Seurat ------------------------------------------------------------------

run_Seurat <- function(sce, 
                       IterCol = "TestIterSet",
                       Iter = 1,
                       PhenotypeCol,
                       MinCells,
                       MaxCells,
                       pc = FALSE,
                       workers = NULL){
  

  if(!is.na(Iter) && !is.null(IterCol)){
    TestIterSet <- SingleCellExperiment::colData(sce)[, IterCol]
    sce <- SummarizedExperiment::subset(sce, , TestIterSet == Iter)
  } 

  SingleCellExperiment::colLabels(sce) <- SingleCellExperiment::colData(sce)[, PhenotypeCol]
  
  seurat.dat <- as_seurat_data(sce)
  
  if(isTRUE(pc) && !is.null(workers)){
    future::plan("multicore", workers = workers)
    options(future.globals.maxSize = 1000 * 1024^2)
  }
  rawresult <- Seurat::FindAllMarkers(
    seurat.dat,
    test.use = "LR",
    min.cells.group = MinCells,
    max.cells.per.ident = MaxCells,
    return.thresh = 1,
    logfc.threshold = 0
  )
  
  rawresult <- tibble::as_tibble(rawresult)
  result <- dplyr::select(
    rawresult,
    gene = gene,
    celltype = cluster,
    spct = pct.1,
    opct = pct.2,
    lfc = avg_log2FC,
    pval = p_val,
    pvaladj = p_val_adj
  )
  
  result <- tibble::as_tibble(result)

  out <- list(result = result, 
              rawresult = rawresult,
              pars = list(Iter = Iter, 
                          Phenotype = PhenotypeCol, 
                          MinCells = MinCells, 
                          pc = pc, 
                          workers = workers))
  
  return(out)
  
}

as_seurat_data <- function(sce) {
  
  stopifnot(is(sce, "SingleCellExperiment"))
  
  sce2 <- sce
  
  # Seurat requires the object to have column names.
  if (is.null(colnames(sce2))) {
    colnames(sce2) <- paste0("cell_", seq_len(ncol(sce2)))
  }
  
  # DelayedArrays are not supported by Seurat.
  if (is(counts(sce2), "DelayedArray")) {
    counts(sce2) <- as(counts(sce2), "dgCMatrix")
  }
  
  # logcounts are converted to a normal matrix as they are not sparse.
  if (is(logcounts(sce2), "DelayedArray")) {
    logcounts(sce2) <- as(logcounts(sce2), "Matrix")
  }
  
  out <- Seurat::as.Seurat(sce2)
  
  # Fix up simulated + Zeisel data for COSG.
  old_assay_name <- names(out@assays)[[1]]
  out <- switch(
    old_assay_name,
    "RNA" = out,
    "endogenous" = SeuratObject::RenameAssays(out, endogenous = "RNA"),
    "originalexp" = SeuratObject::RenameAssays(out, originalexp = "RNA"),
    stop("Unknown assay name")
  )
  
  # Add the cluster ids to the Seurat object.
  Seurat::Idents(out) <- colLabels(sce2)
  
  return(out)
}

# ZIQRANK -----------------------------------------------------------------

run_ZIQRank_Combination <- function(sce, 
                                    IterCol = NULL, 
                                    Iter = NA,
                                    PhenotypeCol, 
                                    tau='umi', 
                                    MTC = "BH",
                                    clean = TRUE,
                                    pc = TRUE,
                                    workers = 2){
  
  if(!is.na(Iter) && !is.null(IterCol)){
    TestIterSet <- SingleCellExperiment::colData(sce)[, IterCol]
    sce <- SummarizedExperiment::subset(sce, , TestIterSet == Iter)
  } 
  
  if(tau == 'umi'){
    taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  }
  if(tau == "nonumi"){
    taus <- seq(0.05, 0.95, by=0.05)
  }
  
  cannot <- data.frame("Phenotype" = SingleCellExperiment::colData(sce)[, PhenotypeCol])
  
  gene_id <- rownames(sce)
  
  if(isTRUE(pc)){
    BPPARAM = BiocParallel::MulticoreParam(workers = workers, tasks = 0L,
                                           stop.on.error = TRUE,
                                           progressbar = FALSE, RNGseed = NULL,
                                           timeout = 30L * 24L * 60L * 60L, exportglobals=TRUE,
                                           log = FALSE, threshold = "INFO", logdir = NA_character_,
                                           resultdir = NA_character_, jobname = "BPJOB",
                                           manager.hostname = NA_character_, manager.port = NA_integer_)
  } else {
    BPPARAM = BiocParallel::SerialParam()
  }
  
  fits <- BiocParallel::bplapply(1:length(gene_id), function(gene_i){
    # create input for ZIQRank
    Dat0 <- cbind(cannot, 
                  expression = as.vector(SingleCellExperiment::logcounts(sce)[gene_id[gene_i],]))
    
    # run ZIQRank
    ziq.result <- try(ZIQRank(formula.logistic = expression ~ Phenotype,
                          formula.quantile = expression ~ Phenotype, 
                          C = "Phenotype", 
                          y_CorD = "C", 
                          data = Dat0,
                          taus = taus, 
                          seed = 2020), 
                      silent = TRUE)
    
    # Combine P-values
    if(inherits(ziq.result, "try-error")){
      ziq.result <- NULL
      combine.result <- NULL
      res <- NULL
    } else{
      combine.result <- try(Combination(input = ziq.result,
                                        method="MinP", 
                                        taus=taus, 
                                        M=10000), 
                            silent = TRUE)
      if(inherits(ziq.result, "try-error")){
        combine.result <- NULL
        res <- ziq.result
      } else{
        res <- c(ziq.result,
                 "pvalue.combine" = combine.result)
      }
    }
    return(res)
  }, BPPARAM = BPPARAM)
  names(fits) <- gene_id
  
  # MTC
  if(!is.null(MTC)){
    pvals <- t(sapply(names(fits), function(i){
      fit <- fits[[i]]
      if(is.null(fit)){
        LRP <- NA
        minP <- NA
      } else {
        log.pval <- fit$pvalue.logistic
        if(is.null(log.pval)){
          LRP <- NA
          minP <- NA
        } else {
          LRP <- log.pval
          combine.pval <- fit$pvalue.combine
          if(is.null(combine.pval)){
            minP <- NA
          } else {
            minP <- combine.pval
          }
        }
      }
      dat <- c(LRP, minP)
      return(dat)
    }, simplify = TRUE))
    qvals <- tibble::as_tibble(data.frame("gene" = names(fits),
                                          "LRP" = pvals[, 1],
                             "minP" = pvals[, 2],
                    "LRP.adj" = p.adjust(pvals[, 1], method = MTC),
                    "minP.adj" = p.adjust(pvals[, 2], method = MTC)))
    
  } else {
    qvals <- NULL
  }
  
  # output
  if(isTRUE(clean)){
    out <- list("Test" = qvals)
  } else {
    out <- list("Fit" = fits,
                "Test" = qvals)
  }
  
  invisible(gc())
  
  return(out)
  
}

 
ZIQRank <- function(formula.logistic, 
                    formula.quantile, 
                    C, 
                    y_CorD="C", 
                    data, 
                    taus=c(0.1, 0.25, 0.5, 0.75, 0.9), 
                    seed=2020){
  
  ## formulas
  
  # logistic model
  
  # arrange logistic model
  mf.logistic = model.frame(formula.logistic, data=data)
  y = model.response(mf.logistic, "numeric")
  b = 1*(y > 0)
  formula.logistic = update(formula.logistic, b ~ .)
  data.logistic = cbind(data, b)
  mf.logistic = model.frame(formula.logistic, data=data.logistic)
  
  # locate C in logistic model
  namesx = all.vars(formula.logistic)[-1]
  condition.loc = which(namesx %in% C) 
  
  # arrange logistic null model if C is not a single covariate
  if (length(condition.loc) > 1){
    mul.logistic = T
    namesx.null = setdiff(namesx, C)
    if (length(namesx.null) == 0){
      formula.logistic.null = as.formula( "b ~ 1" )
    } else formula.logistic.null = as.formula( paste( "b ~", paste(namesx.null, collapse = "+") ) )
    mf.logistic.null = model.frame(formula.logistic.null, data=data.logistic)
  } else mul.logistic = F
  
  
  # quantile model
  
  # determine elements in quantile model, create the "positive subset"
  namey = all.vars(formula.quantile)[1]
  namesx = all.vars(formula.quantile)[-1]
  namesx.score = setdiff(namesx, C)
  data.quantile = data[b==1, ]
  if (y_CorD == "D"){ # perturbation if response is count
    set.seed(seed)
    data.quantile[, namey] = dither(data.quantile[, namey], type = "right", value = 1)
  } 
  
  # extract C
  formula.quantile = as.formula( paste( namey, "~", paste(C, collapse = "+") ) )
  mf.quantile = model.frame(formula.quantile, data=data.quantile)
  c = model.matrix(attr(mf.quantile, "terms"), data=mf.quantile)[, -1]
  if (is.null(dim(c))){ # determine whether C is a single covariate, also continuous or binary
    single_CorB=T
  } else single_CorB=F
  
  # arrange quantile null model, and extract Z
  if (length(namesx.score) == 0){
    formula.quantile = as.formula( paste( namey, "~ 1" ) )
  } else formula.quantile = as.formula( paste( namey, "~", paste(namesx.score, collapse = "+") ) )
  mf.quantile = model.frame(formula.quantile, data=data.quantile)
  z = model.matrix(attr(mf.quantile, "terms"), data=mf.quantile)
  
  
  
  ## set up parameters
  m = length(y) # total sample size
  width = length(taus) # size of tau
  zerorate = length(which(b == 0)) / m # rate of 0's
  
  
  
  ## compute p-values from the marginal tests
  
  if (single_CorB == T){ # when C is a single covariate, either continuous or binary
    
    # logistic, wald test
    mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
    pvalue.logistic = summary(mod.logistic)$coef[condition.loc+1, 4]
    
    # estimate quantiles of y|y>0 | H0
    rq0 = rq(mf.quantile, tau=taus)
    qpred0 = predict(rq0)
    
    # project C on the space of intercept and Z
    C.star = c - z %*% solve( (t(z) %*% z) ) %*% t(z) %*% c
    
    # compute the rank-score test stats, and its covariance matrix
    RS = unlist( lapply(1:width, function(kk){ sum( (taus[kk] - (data.quantile[, namey] < as.matrix(qpred0, ncol=width)[, kk]))*C.star ) / sqrt(m) }) )
    
    if (width == 1){
      cov.RS = taus*(1 - taus)
    } else {
      cov.RS = matrix(0, ncol=width, nrow=width)
      for (kk in 1:(width-1)){
        for (ll in (kk+1):width){
          cov.RS[kk, ll] = min(taus[kk], taus[ll]) - taus[kk]*taus[ll]
        }
      }
      cov.RS = cov.RS + t(cov.RS) + diag(taus*(1 - taus))
    }
    
    Sigma.hat = cov.RS * sum( C.star^2 ) / m
    if (width == 1){
      sigma.hat = sqrt( Sigma.hat )
    } else {
      sigma.hat = sqrt( diag(Sigma.hat) )
    }
    
    # marginal p-value in quantile regression 
    pvalue.quantile = 2*( 1 - pnorm( abs( RS / sigma.hat ) ) ) 
    
  } else { # when C is a set of covariates, or a single covariate with multiple categories
    
    # logistic, score test
    if (mul.logistic != T){
      mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
      pvalue.logistic = anova(mod.logistic, test="Rao")$`Pr(>Chi)`[condition.loc+1] 
    } else {
      mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
      mod.logistic.null = glm(mf.logistic.null, family=binomial(link = 'logit'))
      pvalue.logistic = anova(mod.logistic.null, mod.logistic, test="Rao")$`Pr(>Chi)`[2]
    }
    
    # estimate quantiles of y|y>0 | H0
    rq0 = rq(mf.quantile, tau=taus)
    qpred0 = predict(rq0)
    
    # project C on the space of intercept and Z
    C.star = c - z %*% solve( (t(z) %*% z) ) %*% t(z) %*% c
    
    # compute the rank-score test stats, and its covariance matrix
    RS = lapply(1:width, function(kk){ apply( (taus[kk] - (data.quantile[, namey] < as.matrix(qpred0, ncol=width)[, kk]))*C.star, 2, sum ) / sqrt(m) })
    
    df = ncol(C.star)
    tmp = t(C.star) %*% C.star / m
    
    var.RS = NULL
    for (kk in 1:width){
      var.RS[[kk]] = taus[kk] * (1 - taus[kk]) * tmp
    }
    
    Sigma.hat  = NULL
    for (kk in 1:width){
      temp = NULL
      for (ll in 1:width){
        temp = cbind( temp, ( min(taus[kk], taus[ll]) - taus[kk]*taus[ll] )*tmp )
      }
      Sigma.hat = rbind(Sigma.hat, temp)
    }
    
    # marginal p-value in quantile regression 
    pvalue.quantile = NULL
    for (kk in 1:width){
      stat = t( RS[[kk]] ) %*% solve(var.RS[[kk]]) %*% RS[[kk]]
      pvalue.quantile[kk] = 1 - pchisq(stat, df=df)
    }
    
  }
  
  return(list(pvalue.logistic=pvalue.logistic, 
              pvalue.quantile=pvalue.quantile, 
              Sigma.hat=Sigma.hat,
              zerorate=zerorate,
              taus=taus))
  
}

Combination <- function(input,
                        method="MinP", 
                        taus=c(0.1, 0.25, 0.5, 0.75, 0.9), 
                        M=10000){
  
  ## check
  
  # choose either MinP or Cauchy
  if (method != "MinP" & method != "Cauchy"){
    stop("Please choose 'MinP' or 'Cauchy', no other options.")
  }
  
  # taus and ind should match
  if (!all(taus %in% input$taus)){
    stop("taus should be a subset of that taus used to produce input.")
  }
  
  
  ## whether from single_CorB=T or not
  if (!is.null( ncol(input$Sigma.hat) )){
    if (length(input$pvalue.quantile) != ncol(input$Sigma.hat)){
      single_CorB = F
      df = ncol(input$Sigma.hat) / length(input$pvalue.quantile)
    } else single_CorB = T
  } else single_CorB = T
  
  
  ## compute the aggregate p-value
  width = length(taus)
  
  ind = match(taus, input$taus)
  pvalue.quantile = input$pvalue.quantile[ind]
  
  
  if (method == "MinP"){
    
    # t.obs in the minp test
    t.obs = min(input$pvalue.logistic, pvalue.quantile)
    
    if (single_CorB != F){
      
      if (is.null( ncol(input$Sigma.hat) )){
        Sigma.hat = input$Sigma.hat[ind]
      } else Sigma.hat = input$Sigma.hat[ind, ind]
      
      
      # the (1 - t.obs/2)th percentile of the statistics in quantile regression, normal
      if (is.null( ncol(input$Sigma.hat) )){
        sigma.hat = sqrt( Sigma.hat )
      } else sigma.hat = sqrt( diag(Sigma.hat) )
      qmin.quantile = qnorm((1-t.obs/2), mean=0, sd=sigma.hat)
      
      # MC, to estimate the probability that the absolute of joint statistics in quantile regression < each threshold
      beta.sim = mvrnorm(n=M, mu=rep(0, width), Sigma=Sigma.hat)
      prob.quantile = mean( apply(beta.sim, 1, function(z){ all( abs(z) < qmin.quantile ) }) )
      
    } else {
      
      index = unlist( lapply(ind, function(kk){ ((kk-1)*df+1):(kk*df) }) ) 
      Sigma.hat = input$Sigma.hat[index, index]
      
      # the (1 - t.obs)th percentile of the statistics in quantile regression, chisq
      qmin.quantile = rep(qchisq(1-t.obs, df=df), width)
      
      # MC, to estimate the probability that the joint statistics in quantile regression < each threshold
      beta.sim = mvrnorm(n=M, mu=rep(0, width*df), Sigma=Sigma.hat)
      
      prob.quantile = mean( apply(beta.sim, 1, function(z){ 
        obs.sim = NULL
        for (kk in 1:width){
          obs.sim[kk] = t( z[((kk-1)*df+1):(kk*df)] ) %*% solve(Sigma.hat[((kk-1)*df+1):(kk*df), ((kk-1)*df+1):(kk*df)]) %*% z[((kk-1)*df+1):(kk*df)]
        }
        all( obs.sim < qmin.quantile )
      }) )
      
    }
    
    # final p-value
    pvalue = 1 - (1 - t.obs) * prob.quantile
    
  } else {
    
    # weights for quantile tests
    w = taus*(taus <= 0.5) + (1-taus)*(taus > 0.5)
    w = w / sum(w) * (1 - input$zerorate)
    
    stats.cauchy = input$zerorate * tan( (0.5-input$pvalue.logistic)*pi ) + sum ( w * tan( (0.5-pvalue.quantile)*pi ) )
    
    # final p-value
    pvalue = 1 - pcauchy(stats.cauchy)
    
  }
  
  return(pvalue)
  
}




# MARKER GENE SELECTION ---------------------------------------------------

process_marker_testing <- function(DetectRes,
                                   SeuratRes,
                                   ZIQRankRes,
                                   CLData,
                                   Selection = list("LRP.adj" = 0.01, 
                                                    "minP.adj" = 0.01,
                                                    'pvaladj' = 0.05,
                                                    'lfc' = 10^-2,
                                                    'dpct' = 10^-2,
                                                    'slmean' = log(10^-2))){
  # process cutoff res
  CutoffRes.CT.L <- sapply(names(DetectRes$CT), function(ct){
    DetectRes$CT[[ct]]$Cutoff
  }, USE.NAMES = T, simplify = F)
  CutoffRes.CT <- data.frame(do.call('rbind', CutoffRes.CT.L))
  colnames(CutoffRes.CT) <- "min.spct"
  CutoffRes.CT[,"celltype"] <- rownames(CutoffRes.CT)
  CutoffRes.CT <- tibble::as_tibble(CutoffRes.CT)
  
  # process detect res
  DetectRes.CT.L <- sapply(names(DetectRes$CT), function(ct){
    DetectRes$CT[[ct]]$Data[, c("Gene", "LMean", "LSE")]
  }, USE.NAMES = T, simplify = F)
  DetectRes.CT <- data.table::rbindlist(DetectRes.CT.L, idcol = 'celltype') %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(gene = Gene,
                  slmean = LMean,
                  slse = LSE) %>% 
    dplyr::left_join(CLData, by = 'celltype')   # add germlayer information
  
  # combine ziqrankres and seuratres and gene estimates
  DETestRes <- dplyr::left_join(SeuratRes$result, ZIQRankRes$Test, by = "gene") %>% 
    dplyr::mutate(dpct = spct - opct, .after = opct) 

  # annotate ZIQRank and Seurat Test Result
  DETestRes <- DETestRes %>% 
    dplyr::mutate(sigZIQRank = dplyr::case_when(LRP.adj <= Selection[["LRP.adj"]] &
                                                minP.adj <= Selection[["minP.adj"]] ~ TRUE,
                                                .default = FALSE)) %>% 
    dplyr::mutate(sigSeurat = dplyr::case_when(pvaladj < Selection[["pvaladj"]] ~ TRUE,
                                               .default = FALSE)) %>% 
    dplyr::select(-LRP, -LRP.adj, -minP, -minP.adj, -pval, -pvaladj)

  # detected genes with spct above cutoff
  DETestRes <- DETestRes %>% 
  dplyr::left_join(CutoffRes.CT, by = "celltype") %>% 
    dplyr::mutate(detected = dplyr::case_when(spct >= min.spct ~ TRUE,
                                              .default = FALSE)) %>% 
    dplyr::select(-min.spct)
  
  # expressed genes with slmean above threshold
  DETestRes <- DETestRes %>% 
    dplyr::left_join(DetectRes.CT, 
                     by = c('gene' = 'gene', 'celltype' = 'celltype')) %>% 
    dplyr::mutate(expressed = dplyr::case_when(slmean > Selection[["slmean"]] ~ TRUE,
                                               .default = FALSE))
  
  # indicator when gene is both expressed and detected
  DETestRes <- DETestRes %>% 
    dplyr::mutate(NMarked = dplyr::case_when(detected == TRUE & expressed == TRUE ~ TRUE,
                                             .default = FALSE))
    
  # annotate marker gene expression direction (up, down) with minimal lfc and difference in detection
  DETestRes <- DETestRes %>% 
    dplyr::mutate(direction = dplyr::case_when(dpct > Selection[["dpct"]] & lfc > Selection[["lfc"]] ~ "up", 
                                               dpct < -Selection[["dpct"]] & lfc < -Selection[["lfc"]] ~ "down",
                                               .default = "none"))
  
  # number of cell types marker per direction
  DETestRes <- DETestRes %>% 
    dplyr::group_by(gene, direction) %>%
    dplyr::mutate(DNMarked = sum(NMarked, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(DNMarked = dplyr::case_when(direction == "up" ~ DNMarked,
                                              direction == "down" ~ -DNMarked,
                                              direction == 'none' ~ 0))
  
  # define marker type (unique, multiple) based on detected/expressed criterium
  # per cell type and germlayer
  tmp <- DETestRes %>% 
    dplyr::group_by(gene) %>%
    dplyr::mutate(CMarked = sum(NMarked, na.rm = T)) %>% 
    dplyr::mutate(ctype = dplyr::case_when(CMarked == 0 ~ "none",
                                           CMarked == 1 ~ "unique",
                                           CMarked >= 2 ~ "multiple")) %>% 
    dplyr::select(-CMarked)
  
  tmp2 <- tmp %>% 
    dplyr::group_by(gene, germlayer) %>%
    dplyr::summarise(GGMarked = sum(NMarked, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(GMarked = dplyr::case_when(GGMarked != 0 ~ TRUE,
                                              .default = FALSE)) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(GNMarked = sum(GMarked, na.rm = T)) %>% 
    dplyr::mutate(gtype = dplyr::case_when(GNMarked == 0 ~ "none",
                                           GNMarked == 1 ~ "unique",
                                           GNMarked >= 2 ~ "multiple")) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(gene, gtype)
  
  MarkerRes <- dplyr::left_join(tmp, tmp2, by = 'gene') %>% 
    dplyr::select(gene, NMarked,
                  celltype, ctype,
                  germlayer, gtype,
                  direction, DNMarked,  
                  spct, detected, opct, dpct,
                  slmean, expressed, slse,
                  sigZIQRank, sigSeurat, lfc)
  
  # output object
  Out <- list("Markers" = MarkerRes,
              "Settings" = Selection)
  
  return(Out)
  
}


# PREDICTION PERFORMANCE --------------------------------------------------

predict_celltype_markers <- function(scelist,
                                     splitlist,
                                     phenocol,
                                     markers,
                                     ct_level,
                                     type,
                                     biotype,
                                     nogenes,
                                     k_knn){
  
  # marker genes
  MarkerList <- sapply(names(markers$Celltype), function(ct){
    markerdat <- markers$Celltype[[ct]]
    if(ct_level %in% c('celltype', 'germlayer')){
      if(type == "all"){
        type <- c("unique", "multiple")
      }
      if(ct_level == "celltype"){
        markerdat <- markerdat %>% 
          dplyr::filter(ctype %in% type)
      }
      if(ct_level == "germlayer"){
        markerdat <- markerdat %>% 
          dplyr::filter(gtype %in% type)
      }
    }
    if(!is.null(biotype)){
      markerdat <- markerdat %>% 
        dplyr::filter(gene %in% biotype)
    }
    if(nogenes <= nrow(markerdat)){
      markergenes <- as.vector(markerdat[1:nogenes, "gene", drop = T])
    } else {
      markergenes <- as.vector(markerdat[, "gene", drop = T])
    }
    return(markergenes)
  }, simplify = F, USE.NAMES = T)
  
  MarkerData <- tibble::tibble(celltype = rep(names(MarkerList), sapply(MarkerList, length)),
                               gene = unlist(unname(MarkerList)))
  
  MarkerGenes <- unique(unlist(unname(MarkerList)))
  
  TotalMarkerGenes <- length(MarkerGenes)

  sce.marker <- list('train' =  scelist$train[rownames(scelist$train) %in% MarkerGenes,], 
                     'test' =  scelist$test[rownames(scelist$test) %in% MarkerGenes,])
  
  splitgroup <- list('train' = SingleCellExperiment::colData(scelist$train)[, splitlist$train],
                     'test' = SingleCellExperiment::colData(scelist$test)[, splitlist$test])
   
  splitgroups.train <- unique(splitgroup$train[!is.na(splitgroup$train)])
  if(is.logical(splitgroup$train)){
    splitgroups.train <- TRUE
  }
  splitgroups.test <- unique(splitgroup$test[!is.na(splitgroup$test)])
  if(is.logical(splitgroups.test)){
    splitgroups.test <- TRUE
  }

  splitgroups <- tidyr::expand_grid(train = as.character(splitgroups.train), 
                                    test = as.character(splitgroups.test)) %>% 
    dplyr::mutate(trainname = dplyr::case_when(train == "TRUE" ~ 'all', .default = as.character(train)),
                  testname = dplyr::case_when(test == "TRUE" ~ 'all', .default = as.character(test))) %>% 
    tidyr::unite("set", c(trainname, testname), sep = "_vs_")
  class(splitgroups[, 'train', drop = T]) <- class(splitgroups.train)
  class(splitgroups[, 'test', drop = T]) <- class(splitgroups.test)
  
  # KNN classification
  predict.res <- lapply(1:nrow(splitgroups), function(i){
    sce.train <- subset(sce.marker[['train']], , splitgroup[['train']] == splitgroups[i, 'train', drop = T])
    sce.test <- subset(sce.marker[['test']], , splitgroup[['test']] == splitgroups[i, 'test', drop = T] )

    # Collect all levels of the variable from both datasets and combine into one vector
    complete_factor_levels <- c(SingleCellExperiment::colData(sce.train)[, phenocol, drop = T] %>% factor %>% levels, 
                                SingleCellExperiment::colData(sce.test)[, phenocol, drop = T] %>% factor %>% levels) %>% unique
    
    # Assign the complete factor levels to the variable in both datasets
    train.annot <- factor(SingleCellExperiment::colData(sce.train)[, phenocol], levels = complete_factor_levels)
    test.annot <- factor(SingleCellExperiment::colData(sce.test)[, phenocol], levels = complete_factor_levels)

    dat.train <- t(as.matrix(SingleCellExperiment::logcounts(sce.train)))
    dat.test <- t(as.matrix(SingleCellExperiment::logcounts(sce.test)))
    
    res.knn <- try(FNN::knn(dat.train, 
                    dat.test,
                    train.annot, 
                    k = k_knn, 
                    prob=TRUE),
               silent = T)
    if(inherits(res.knn, "try-error")){
      res <- NULL
    } else{
      res <- list(true.annot = test.annot,
                  predict.annot = res.knn,
                  test.annot = train.annot)
    }
    
    return(res)
  })
  names(predict.res) <- splitgroups[, 'set', drop = T]
  # kick out failed knns
  predict.res <- predict.res[!sapply(predict.res, is.null)]
  
  # performance metrics
  metric.res <- lapply(seq_along(predict.res), function(i){
    cm <- suppressWarnings(caret::confusionMatrix(data = predict.res[[i]]$predict.annot,
                                                  reference = predict.res[[i]]$true.annot,
                                                  dnn = c("Prediction", "Truth"),
                                                  mode = "prec_recall"))
    return(cm)
  })
  names(metric.res) <- names(predict.res)
  
  # do not look at accuracy, use only precision, recall, and f1 score, alternative to accuracy is kappa!
  # https://machinelearningmastery.com/tactics-to-combat-imbalanced-classes-in-your-machine-learning-dataset/
  
  # wrangle metrics
  # overall
  overall.res <- sapply(names(predict.res), function(i){
    .run.external.calc(Prediction = predict.res[[i]]$predict.annot,
                       Truth = predict.res[[i]]$true.annot)
  }, USE.NAMES = T, simplify = F)
  overall.res <- data.table::rbindlist(overall.res, idcol = "Set")
  # kappa (accuracy for unbalanced data)
  kappa.res <- sapply(names(metric.res), function(i){
    data.frame("Kappa" = metric.res[[i]]$overall["Kappa", drop='F'], row.names = NULL)
  }, simplify = F, USE.NAMES = T) 
  kappa.res <- data.table::rbindlist(kappa.res, idcol = "Set")

  overall.res <- dplyr::full_join(overall.res, kappa.res, by = "Set")
  
  # per cell type
  # recall, precision, f1 score, etc
  ct.res <- sapply(names(metric.res), function(i){
    Metrics <- cbind.data.frame("Celltype" = gsub(pattern = "Class: ", 
                                                  replacement = "", 
                                                  rownames(metric.res[[i]]$byClass)),
                                metric.res[[i]]$byClass)
    rownames(Metrics) <- NULL
    return(Metrics)
  }, simplify = F, USE.NAMES = T) 
  ct.res <- data.table::rbindlist(ct.res, idcol = "Set")
  
  Out <- list("Overall" = overall.res,
              "CellType" = ct.res,
              "Marker" = list("Data" = MarkerData,
                              "Number" = TotalMarkerGenes))
  
  return(Out)
  
}

.run.external.calc <- function(Prediction,
                               Truth){
  
  # tabulate prediction cell type and true cell type
  tbl = table(Prediction, Truth)
  # tbl = table(Truth, Prediction)
  conv_df = as.data.frame.matrix(tbl)
  
  # compute marginals of confusion matrix
  tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))
  tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))
  tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))
  fp = tp_plus_fp - tp
  fn = tp_plus_fn - tp
  tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn
  
  # adjusted rand
  prod_comb = (tp_plus_fp * tp_plus_fn) / gmp::asNumeric(gmp::chooseZ(length(Truth), 2))
  mean_comb = (tp_plus_fp + tp_plus_fn)/2
  
  # purity
  tmp_pur = apply(conv_df, 1, max)
  res_purity = sum(tmp_pur)/length(Truth)
  
  # entropy
  tmp_entropy = sum(apply(conv_df, 2, function(x) .entropy.formula(x)))
  res_entropy = -(1/(sum(tbl) * log2(length(unique(Truth))))) * tmp_entropy
  
  # nmi, nvi
  mutual_information = 0
  joint_entropy = 0
  for (i in 1:nrow(conv_df)) {
    for (j in 1:ncol(conv_df)) {
      if (conv_df[i, j] > 0) {
        joint_entropy = joint_entropy + (-((conv_df[i, j]/sum(tbl)) * log2(conv_df[i, j]/sum(tbl))))
        mutual_information = mutual_information + ((conv_df[i, j]/sum(tbl)) * log2(as.numeric(gmp::as.bigz(as.numeric(sum(tbl)) * as.numeric(conv_df[i, j]))/gmp::as.bigz(as.numeric(sum(conv_df[i, ])) * as.numeric(sum(conv_df[, j]))))))
      }
    }
  }
  entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x)/sum(tbl)) * log2(sum(x)/sum(tbl))))
  entr_class = sum(apply(conv_df, 2, function(x) -(sum(x)/sum(tbl)) * log2(sum(x)/sum(tbl))))
  unq_true = unique(Truth)
  unq_clust = unique(Prediction)
  if (length(unq_true) == 1 && length(unq_clust) == 1) {
    NMI = 1
    NVI = 1
  } else {
    NMI = (mutual_information/((entr_cluster + entr_class)/2))
    NVI = 1 - (mutual_information/joint_entropy)
  }
  VAR_INFO = (entr_cluster + entr_class) - 2 * mutual_information
  
  # precision
  prec = tp/(tp + fp)
  
  # recall
  rec = tp/(tp + fn)
  
  # specificity
  specif <- tn/(tn + fp)
  
  # sensitivity
  sensitive <- tp/(tp + fn)
  
  # F measure
  fmeasure <- 2 * ((prec * rec)/(prec + rec))
  
  # rand index
  RI <- (tp + tn)/(tp + fp + fn + tn)
  # adjusted rand index
  ARI <- (tp - prod_comb)/(mean_comb - prod_comb)
  
  # jaccard index
  JI <- tp/(tp + fp + fn)
  
  # fowlkes mallows index
  FMI <- sqrt((tp/((tp + fp))) * (tp/(tp + fn)))
  
  # mirkin metric
  MM <- 2 * (fp + fn)
  
  # return object
  res <- list(specif, sensitive, prec, rec, fmeasure, RI, ARI, JI, FMI, MM, res_purity, res_entropy, NMI, VAR_INFO, NVI)
  names(res) <- c('Specificity', 'Sensitivity', 'Precision', 'Recall', 'F1-Score',
                  'RI', 'ARI', 'JI', 'FMI', 'MM', 'Purity', 'Entropy', 'NMI', 'VI', 'NVI')
  
  invisible(gc())
  
  return(res)
  
}

.entropy.formula <- function (x_vec)
{
  vec = rep(NA, length(x_vec))
  for (i in 1:length(x_vec)) {
    if (x_vec[i] == 0) {
      vec[i] = 0
    }
    else {
      vec[i] = ((x_vec[i]) * log2(x_vec[i]/sum(x_vec)))
    }
  }
  return(vec)
}

summarise_predict <- function(input,
                              yval = "F1",
                              xval = "nogenes"){
  
  # SORT INPUT
  dat.srt <- stats::sortedXyData(input[, xval, drop = T], 
                          input[, yval, drop = T], 
                          input)
  
  dat.val <- tibble(x = input[, xval, drop = T], 
                 y = input[, yval, drop = T]) %>% 
    dplyr::group_by(x) %>% 
    dplyr::summarise(ymean = mean(y, na.rm = T),
                     yse = sd(y , na.rm = T)/sqrt(n()))
  
  # RUN Asymptotic FIT
  unconstrained.ssasymp.fit <- try(nls(y ~ SSasymp(x, Asym, R0, lrc),
                                       data = dat.srt),
                      silent = TRUE)
  
  if(inherits(unconstrained.ssasymp.fit, "try-error")){
    out <- NULL
  } else {
    fit.vals <- stats::fitted(unconstrained.ssasymp.fit)
    # ESTIMATE FIT VALUES
    cimat1 <- investr::predFit(unconstrained.ssasymp.fit, 
                               interval = "confidence", 
                               level = 0.95) %>% 
      data.frame() %>% 
      dplyr::rename(yhat = fit, lowerCI = lwr, upperCI = upr)
    
    pimat1 <- investr::predFit(unconstrained.ssasymp.fit, 
                               interval = "prediction",
                               level = 0.95) %>% 
      data.frame() %>% 
      dplyr::rename(yhat = fit, lowerPI = lwr, upperPI = upr) %>% 
      dplyr::select(lowerPI, upperPI)
    
    asymptote.val <- stats::coef(unconstrained.ssasymp.fit)[1]
    r0.val <- stats::coef(unconstrained.ssasymp.fit)[2]
    start.val <- stats::fitted(unconstrained.ssasymp.fit)[1]
    # knee point
    knee.val <- KneeArrower::findCutoff(x =1:length(fit.vals),
                                          y = fit.vals,
                                          method = "first",
                                          frac.of.steepest.slope = 0.1)$y
    half.val <- KneeArrower::findCutoff(x =1:length(fit.vals),
                                        y = fit.vals,
                                        method = "first",
                                        frac.of.steepest.slope = 0.5)$y
    nogenes.max <- floor(NLSstClosestX(dat.srt, asymptote.val))
    nogenes.half <- ceiling(NLSstClosestX(dat.srt, half.val))
    nogenes.knee <- ceiling(NLSstClosestX(dat.srt, knee.val))
    
    unconstrained.ssasymp.fit$Fit <- dplyr::bind_cols(dat.val, 
                                                      cimat1,
                                                      pimat1)
    
    # OUTPUT
    res <- data.frame(point = c("start", 'half', 'knee', "asymptote"),
                      xval = c(1, nogenes.half, nogenes.knee, nogenes.max),
                      yval = c(start.val, half.val, knee.val, asymptote.val),
                      percyval = c((start.val/asymptote.val) * 100,
                                   (half.val/asymptote.val) * 100,
                                   (knee.val/asymptote.val) * 100,
                                   100)
    ) 
    
    out <- list("Data" = res,
                "Fit" = unconstrained.ssasymp.fit)
  }
  
  return(out)
}
