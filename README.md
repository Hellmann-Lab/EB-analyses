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

[Shiny app](https://shiny.bio.lmu.de/Cross_Species_CellType/)


### Cell type specificity analysis

### Marker gene analysis


## Figures

