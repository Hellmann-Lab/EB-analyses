#!/bin/bash

CELLRANGER=$1
CB_OUT=$2

OBS_CELLS=$(zcat $CELLRANGER/filtered_feature_bc_matrix/barcodes.tsv.gz | wc -l | awk '{print $1}')
EXP_CELLS=$(($OBS_CELLS - 1000))
TOTAL_DROP=$(($OBS_CELLS + 15000))
echo "Expected cells: $EXP_CELLS"
echo "Total droplets included: $TOTAL_DROP"

echo "Running CellBender"
cellbender remove-background \
	     --input $CELLRANGER/raw_feature_bc_matrix.h5 \
	     --output $CB_OUT \
	     --expected-cells $EXP_CELLS \
	     --total-droplets-included $TOTAL_DROP \
	     --fpr 0.01 \
	     --epochs 150
