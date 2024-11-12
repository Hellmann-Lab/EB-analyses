# Primate embryoid body analysis
This repository contains the necessary code to reproduce the analysis for the manuscript "Identification and comparison of orthologous cell types from primate embryoid bodies shows limits of marker gene transferability".

## Analysis
Raw sequencing files and count matrices are available at [link to GEO]. Processed files and intermediate analysis results can be downloaded from [link to Zenodo].

### scRNA-seq processing
Code used for the initial processing of the scRNA-seq data per species is collected in `analysis/01_scRNAseq_processing`:

- `cellranger`: scripts for alignment of sequencing data to different reference genomes
- `CellBender.sh`: removal of background RNA
- `species_demultiplexing`
- `processing_per_species.Rmd`: Workflow to prepare the count data of each species including the following steps:
  * QC and filtering
  * Doublet detection
  * Normalization
  * Preliminary cell type classification with a reference dataset of human EBs ([Rhodes et al. 2022](https://doi.org/10.7554/eLife.71361))
  * Data integration per species: combining data from different experiments
  * Data integration across species: combining cells from all species for exploratory analysis

### Cell type annotation
To assign comparable cell types across species we used a two step approach: 1) We identify orthologous cell clusters across species based on reciprocal classification, and 2) we curate and annotate the clusters to cell types using an interactive Shiny app. 

Code to perform the cluster matching across species (1) can be found in `analysis/02_celltype_annotation`. Briefly, this includes the following steps:

- High resolution clustering per species
- Score similarity of clusters across and within species based on reciprocal classification with SingleR
- Hierarchical clustering to group similar clusters across species into orthologous clusters
- Merge orthologous clusters with very similar expression profiles

Furthermore, we provide an interactive [Shiny app](https://shiny.bio.lmu.de/Cross_Species_CellType/), that can be used to explore different parameter choices for the orthologous cell type assignment pipeline, identify and visualize marker gene expression to aid cell type annotation and explore our cross-species dataset of embryoid body differentiation.

### Cell type specificity and expression conservation
Scripts to perform the cell type specificity and expression conservation analysis per gene are collected in `analysis/03_celltype_specificity`. 

`run_glmgampoi.R` is used to estimate distributional characteristics per gene.

`analysis_gene_detection.R` and `run_gene_detection.R` are used to determine species- and cell type-specific thresholds, when to call a gene expressed based on the observed expression fraction. 

`summarize_specificity_and_conservation.R` is used to assess the cell type specificity and expression conservation of each gene based on the binarized expression patterns.

### Marker gene analysis
Scripts for marker gene identification and comparison across species can be found in `analysis/04_marker_gene_transferability`. 

`run_ziqrank.R` and `run_seurat_marker.R` are used to detect marker genes. 

kNN classification and evaluation is performed with `run_predict_marker.R` and `run_predict_summary.R`. 

`marker_gene_analysis.R` contains the comparison of marker gene overlap across species and summarizes the kNN classification results. 

## Figures
Scripts to recreate paper figures in `figure_scripts`:

- `Figure1.R`: Figure 1DE and related supplementary figure S3
- `Figure2.R`: Figure 2 and related supplementary figures S4-6
- `Figure3.R`: Figure 3 and related supplementary figure S7
- `Figure4.R`: Figure 4 and related supplementary figures S8-9

