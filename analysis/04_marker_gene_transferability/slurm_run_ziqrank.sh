# settings
slurmbash="/data/share/htp/EBgrant/Beate/slurm/bash/";
slurmerr="/data/share/htp/EBgrant/Beate/slurm/err/";
slurmout="/data/share/htp/EBgrant/Beate/slurm/out/";

cd $slurmbash;

cat /data/share/htp/EBgrant/Beate/data/output/ZIQRank/param_clone.txt | while read ID input \
outname itercol \
phenotypecol tau mtc \
clean outfilepath \
pc ncores mem \
iter;

do
# MAKING THE HEADER
echo '#!/bin/bash' > $ID.ziq.sh
echo '#SBATCH -n 1' >> $ID.ziq.sh
echo '#SBATCH --chdir='$slurmbash'' >> $ID.ziq.sh
echo '#SBATCH --error='$slurmerr/$ID'.ziq.%J.err' >> $ID.ziq.sh
echo '#SBATCH --output='$slurmout/$ID'.ziq.%J.out' >> $ID.ziq.sh
echo '#SBATCH --cpus-per-task='$ncores'' >> $ID.ziq.sh
echo '#SBATCH --mem='$mem'' >> $ID.ziq.sh
echo '#SBATCH --nodelist=gorilla3' >> $ID.ziq.sh
echo '#SBATCH -p normal' >> $ID.ziq.sh

# THE ACTUAL COMMANDS
echo "echo $ID" >> $ID.ziq.sh
echo "srun Rscript /data/share/htp/EBgrant/Beate/scripts/run_ziqrank.R  --input=$input --phenotype=$phenotypecol --iteration=$itercol --iter=$iter --tau=$tau --mtc=$mtc --pc=$pc --workers=$ncores --clean=$clean --outpath=$outfilepath --outname=$outname" >> $ID.ziq.sh

sbatch $ID.ziq.sh

done
