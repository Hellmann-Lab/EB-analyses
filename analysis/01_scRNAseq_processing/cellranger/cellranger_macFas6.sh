cellranger count --id=macFas6_run1_1 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run1/EB_10x_scRNA_seq/fastq_trimmed/  \
                      --include-introns=false \
                      --sample=EB_early_pool \
                      --localcores=8 \
                      --localmem=33

cellranger count --id=macFas6_run1_2 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run1/EB_10x_scRNA_seq/fastq_trimmed/  \
                      --include-introns=false \
                      --sample=EB_late_pool \
                      --localcores=8 \
                      --localmem=33

cellranger count --id=macFas6_run2_1 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run2/fastq/deML/  \
                      --include-introns=false \
                      --sample=EB_run2_lane1 \
                      --localcores=8 \
                      --localmem=33

cellranger count --id=macFas6_run2_2 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run2/fastq/deML/  \
                      --include-introns=false \
                      --sample=EB_run2_lane2 \
                      --localcores=8 \
                      --localmem=33

cellranger count --id=macFas6_d8_1_1 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run3/fastq/  \
                      --include-introns=false \
                      --sample=22049a001_01 \
                      --localcores=8 \
                      --localmem=33

cellranger count --id=macFas6_d8_1_2 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run3/fastq/  \
                      --include-introns=false \
                      --sample=22049a002_01 \
                      --localcores=8 \
                      --localmem=33
                    
cellranger count --id=macFas6_d8_2_1 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run3/fastq/  \
                      --include-introns=false \
                      --sample=22049a003_01 \
                      --localcores=8 \
                      --localmem=33
                    
cellranger count --id=macFas6_d8_2_2 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run3/fastq/  \
                      --include-introns=false \
                      --sample=22049a004_01 \
                      --localcores=8 \
                      --localmem=33
                    
cellranger count --id=macFas6_d16_1 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run3/fastq/  \
                      --include-introns=false \
                      --sample=22049a005_01 \
                      --localcores=8 \
                      --localmem=33

cellranger count --id=macFas6_d16_2 \
                      --transcriptome=/data/home/EBgrant/genome_data_new/macFas6 \
                      --fastqs=/data/home/EBgrant/scRNA_run3/fastq/  \
                      --include-introns=false \
                      --sample=22049a006_01 \
                      --localcores=8 \
                      --localmem=33