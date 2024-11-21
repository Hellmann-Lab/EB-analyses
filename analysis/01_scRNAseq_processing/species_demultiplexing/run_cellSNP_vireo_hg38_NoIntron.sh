#!/bin/bash
#SBATCH --error=cs.%J.err
#SBATCH --output=cs.%J.out
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hellmann@bio.lmu.de

sample=$1
echo "processing $sample"

wd=/data/home/EBgrant/scRNA_run3/cellSNP/$sample
mkdir $wd 

#cellsnp
cellranger_outs=/data/home/EBgrant/scRNA_run3/cellranger/$sample/outs
bam=$cellranger_outs/possorted_genome_bam.bam
cellBCs=$cellranger_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
cellSNP_dir=$wd/cellsnp_NoIntron

REGION_VCF=/data/home/EBgrant/scRNA_run1/cellSNP/reference_vcf/species_variants.recode.vcf

cd $wd

source /opt/anaconda3/bin/activate
conda activate vireo

echo "cellsnp-lite -s $bam -b $cellBCs -O $cellSNP_dir -R $REGION_VCF --minMAF 0 --minCOUNT 10 --cellTAG CB --UMItag UB --gzip
"
cellsnp-lite -s $bam -b $cellBCs -O $cellSNP_dir -R $REGION_VCF --minMAF 0 --minCOUNT 10 --cellTAG CB --UMItag UB --gzip


echo "vireo -c $cellSNP_dir -d $REGION_VCF -o $wd/vireo_NoIntron"
vireo -c $cellSNP_dir -d $REGION_VCF -o $wd/vireo_NoIntron
