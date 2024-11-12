# settings
slurmbash="/data/share/htp/EBgrant/Beate/slurm/bash/";
slurmerr="/data/share/htp/EBgrant/Beate/slurm/err/";
slurmout="/data/share/htp/EBgrant/Beate/slurm/out/";

cd $slurmbash;

mem=20000 ;
ncores=1 ;

cat /data/share/htp/EBgrant/Beate/data/output/MarkerSelection/PredictivePerformance/param_prediction.txt | while read ID \
trainsce trainsplitcol \
testsce testsplitcol \
phenocol markers \
biotype \
ctlevel mtype \
nogenes k \
outpath outname \
pc workers ;

do
# MAKING THE HEADER
echo '#!/bin/bash' > $ID.predict.sh
echo '#SBATCH -n 1' >> $ID.predict.sh
echo '#SBATCH --chdir='$slurmbash'' >> $ID.predict.sh
echo '#SBATCH --error='$slurmerr/$ID'.predict.%J.err' >> $ID.predict.sh
echo '#SBATCH --output='$slurmout/$ID'.predict.%J.out' >> $ID.predict.sh
echo '#SBATCH --cpus-per-task='$workers'' >> $ID.predict.sh
echo '#SBATCH --mem='$mem'' >> $ID.predict.sh
echo '#SBATCH -p high' >> $ID.predict.sh
echo '#SBATCH --exclude=gorilla1,gorilla4' >> $ID.predict.sh

# THE ACTUAL COMMANDS
echo "echo $ID" >> $ID.predict.sh
echo "srun Rscript /data/share/htp/EBgrant/Beate/scripts/run_predict_marker.R  --trainsce=$trainsce --trainsplitcol=$trainsplitcol --testsce=$testsce --testsplitcol=$testsplitcol --markers=$markers --level=$ctlevel --type=$mtype --biotype=$biotype --phenocol=$phenocol --nogenes=$nogenes --k=$k --pc=$pc --workers=$workers --outpath=$outpath --outname=$outname" >> $ID.predict.sh

sbatch $ID.predict.sh

done
