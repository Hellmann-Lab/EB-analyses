# settings
slurmbash="/data/share/htp/EBgrant/Beate/slurm/bash/";
slurmerr="/data/share/htp/EBgrant/Beate/slurm/err/";
slurmout="/data/share/htp/EBgrant/Beate/slurm/out/";

cd $slurmbash;

cat /data/share/htp/EBgrant/Beate/data/output/Seurat/param_all.txt | while read ID input \
outname itercol \
phenotypecol \
mincell maxcell \
outfilepath \
pc ncores mem \
iter;

do
# MAKING THE HEADER
echo '#!/bin/bash' > $ID.seurat.sh
echo '#SBATCH -n 1' >> $ID.seurat.sh
echo '#SBATCH --chdir='$slurmbash'' >> $ID.seurat.sh
echo '#SBATCH --error='$slurmerr/$ID'.seurat.%J.err' >> $ID.seurat.sh
echo '#SBATCH --output='$slurmout/$ID'.seurat.%J.out' >> $ID.seurat.sh
echo '#SBATCH --cpus-per-task='$ncores'' >> $ID.seurat.sh
echo '#SBATCH --mem='$mem'' >> $ID.seurat.sh
echo '#SBATCH -p low' >> $ID.seurat.sh
echo '#SBATCH --exclude=gorilla1,gorilla4' >> $ID.seurat.sh

# THE ACTUAL COMMANDS
echo "echo $ID" >> $ID.seurat.sh
echo "srun Rscript /data/share/htp/EBgrant/Beate/scripts/run_seurat_marker.R  --input=$input --phenotype=$phenotypecol --iteration=$itercol --iter=$iter --mincell=$mincell --maxcell=$maxcell --pc=$pc --workers=$ncores --outpath=$outfilepath --outname=$outname" >> $ID.seurat.sh

sbatch $ID.seurat.sh

done
