#!/bin/bash
OUTPUT_DIR=$1
NETWORK=$2
REGULATORS=$3
GENES=$4
PROMOTER=$5
FLAG=$6

## Prepare score bins
printf "Binning promoters based on network scores ... "
python CODE/parse_network_scores.py -a $NETWORK -r $REGULATORS -t $GENES -o ${OUTPUT_DIR}/motif_inference/network_scores/
python CODE/parse_quantized_bins.py -n 20 -i ${OUTPUT_DIR}/motif_inference/network_scores/ -o ${OUTPUT_DIR}/motif_inference/network_bins/
printf "DONE\n"

## Infer FIRE motifs using SLURM array scheme
rm -f ${OUTPUT_DIR}/motif_inference/motif_inference.log
touch ${OUTPUT_DIR}/motif_inference/motif_inference.log

printf "Inferring DNA binding motifs using FIRE ... "
num_regulators=$( wc -l ${REGULATORS} )
sbatch --array=1-${num_regulators}%48 CODE/infer_motifs.sh $REGULATORS $PROMOTER ${OUTPUT_DIR}/motif_inference/network_bins/ ${OUTPUT_DIR}/motif_inference/motif_inference.log

## Check if all motifs are ready
bash CODE/check_inference_status.sh ${OUTPUT_DIR}/motif_inference/motif_inference.log $REGULATORS $FLAG
