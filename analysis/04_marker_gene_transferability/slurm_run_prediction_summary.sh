# settings
slurmbash="/data/share/htp/EBgrant/Beate/slurm/bash/";
slurmerr="/data/share/htp/EBgrant/Beate/slurm/err/";
slurmout="/data/share/htp/EBgrant/Beate/slurm/out/";

cd $slurmbash;

mem=20000 ;
ncores=1 ;

cat /data/share/htp/EBgrant/Beate/data/output/MarkerSelection/PredictivePerformance/param_summary_f1.txt | while read ID \
predfile exprdata markerspecies \
classification celltype \
biotype ctlevel mtype \
ln un metric \
outpath outname ;

do
# MAKING THE HEADER
echo '#!/bin/bash' > $ID.predict.sh
echo '#SBATCH -n 1' >> $ID.predict.sh
echo '#SBATCH --chdir='$slurmbash'' >> $ID.predict.sh
echo '#SBATCH --error='$slurmerr/$ID'.predict.%J.err' >> $ID.predict.sh
echo '#SBATCH --output='$slurmout/$ID'.predict.%J.out' >> $ID.predict.sh
echo '#SBATCH --cpus-per-task='$ncores'' >> $ID.predict.sh
echo '#SBATCH --mem='$mem'' >> $ID.predict.sh
echo '#SBATCH -p low' >> $ID.predict.sh
echo '#SBATCH --exclude=gorilla1,gorilla4' >> $ID.predict.sh
echo '#SBATCH --nodelist=gorilla2' >> $ID.predict.sh

# THE ACTUAL COMMANDS
echo -e "echo $ID" >> $ID.predict.sh
echo -e "srun Rscript /data/share/htp/EBgrant/Beate/scripts/run_predict_summary.R  --prediction=$predfile --metric=$metric --exprdata=$exprdata --markerspecies=$markerspecies --celltype=$celltype --classification=$classification --level=$ctlevel --type=$mtype --biotype=$biotype --lowern=$ln --uppern=$un --outpath=$outpath --outname=$outname" >> $ID.predict.sh

sbatch $ID.predict.sh

done
