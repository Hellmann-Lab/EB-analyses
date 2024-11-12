# settings
slurmbash="slurm/bash/";
slurmerr="slurm/err/";
slurmout="slurm/out/";

cd $slurmbash;

cat data/output/GeneDetection/param_all.txt \
| while read ID sce ggp \
phenotypecol itercol iter \
qse qmean pmax pthreshold \
outpath outname \
pc workers mem verbose;

do
# MAKING THE HEADER
echo '#!/bin/bash' > $ID.detect.sh
echo '#SBATCH -n 1' >> $ID.detect.sh
echo '#SBATCH --chdir='$slurmbash'' >> $ID.detect.sh
echo '#SBATCH --error='$slurmerr/$ID'.detect.%J.err' >> $ID.detect.sh
echo '#SBATCH --output='$slurmout/$ID'.detect.%J.out' >> $ID.detect.sh
echo '#SBATCH --cpus-per-task='$workers'' >> $ID.detect.sh
echo '#SBATCH --mem='$mem'' >> $ID.detect.sh
echo '#SBATCH --exclude=gorilla1' >> $ID.detect.sh
echo '#SBATCH -p low' >> $ID.detect.sh

# THE ACTUAL COMMANDS
echo "echo $ID" >> $ID.detect.sh
echo "srun Rscript analysis/03_celltype_specificity/run_gene_detection.R  --sce=$sce --ggp=$ggp --phenotype=$phenotypecol --iteration=$itercol --iter=$iter --qmean=$qmean --qse=$qse --pmax=$pmax --pthreshold=$pthreshold --pc=$pc --workers=$workers --verbose=$verbose --outpath=$outpath --outname=$outname --verbose=$verbose" >> $ID.detect.sh

sbatch $ID.detect.sh

done
