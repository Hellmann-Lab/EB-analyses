# settings
slurmbash="/data/share/htp/EBgrant/Beate/slurm/bash/";
slurmerr="/data/share/htp/EBgrant/Beate/slurm/err/";
slurmout="/data/share/htp/EBgrant/Beate/slurm/out/";

cd $slurmbash;

cat /data/share/htp/EBgrant/Beate/data/output/glmGamPoi/param_clone.txt | while read ID input \
outname itercol \
phenotypecol \
outfilepath \
mem verbose \
iter;

do
# MAKING THE HEADER
echo '#!/bin/bash' > $ID.glm.sh
echo '#SBATCH -n 1' >> $ID.glm.sh
echo '#SBATCH --chdir='$slurmbash'' >> $ID.glm.sh
echo '#SBATCH --error='$slurmerr/$ID'.glm.%J.err' >> $ID.glm.sh
echo '#SBATCH --output='$slurmout/$ID'.glm.%J.out' >> $ID.glm.sh
echo '#SBATCH --cpus-per-task=2' >> $ID.glm.sh
echo '#SBATCH --mem='$mem'' >> $ID.glm.sh
echo '#SBATCH --exclude=gorilla1' >> $ID.glm.sh
echo '#SBATCH -p normal' >> $ID.glm.sh

# THE ACTUAL COMMANDS
echo "echo $ID" >> $ID.glm.sh
echo "srun Rscript /data/share/htp/EBgrant/Beate/scripts/run_glmgampoi.R  --input=$input --phenotype=$phenotypecol --iteration=$itercol --iter=$iter --verbose=$verbose --outpath=$outfilepath --outname=$outname" >> $ID.glm.sh

sbatch $ID.glm.sh

done
