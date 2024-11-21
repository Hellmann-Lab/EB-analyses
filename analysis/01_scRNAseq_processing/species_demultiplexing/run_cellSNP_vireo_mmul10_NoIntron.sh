#!/bin/bash
#SBATCH --error=cs.%J.err
#SBATCH --output=cs.%J.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G

sample=$1
echo "processing $sample"

wd=/data/home/EBgrant/scRNA_run3/cellSNP/mmul10_NoIntron_$sample
mkdir $wd 

#cellsnp
cellranger_outs=/data/home/EBgrant/scRNA_run3/cellranger/mmul10_$sample/outs
bam=$cellranger_outs/possorted_genome_bam.bam
cellBCs=/data/home/EBgrant/scRNA_run3/cellSNP/hg38_NoIntron_$sample/cyno_cells.txt
REGION_VCF=/data/home/EBgrant/scRNA_run1/cellSNP/create_reference/cyno_rhesus_variants_mmul10_early.cellSNP/cellSNP.base.vcf.gz

#vireo
n_donor=2
cd $wd

source /opt/anaconda3/bin/activate
conda activate vireo

mkdir $wd/cellsnp
echo "cellsnp-lite -s $bam -b $cellBCs -O $wd/cellsnp -R $REGION_VCF --minMAF 0 --minCOUNT 10 --cellTAG CB --UMItag UB --gzip"
cellsnp-lite -s $bam -b $cellBCs -O $wd/cellsnp -R $REGION_VCF --minMAF 0.05 --minCOUNT 10 --cellTAG CB --UMItag UB --gzip

echo "vireo -c $wd/cellsnp -N $n_donor -o $wd/vireo"
vireo -c $wd/cellsnp -N $n_donor -o $wd/vireo


