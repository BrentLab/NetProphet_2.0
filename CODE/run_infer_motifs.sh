#!/bin/bash
OUTPUT_DIR=$1
NETWORK=$2
REGULATORS=$3
GENES=$4
TF_LIST=$5
PROMOTER=$6
FLAG=$7

##Prepare score bins
python CODE/parse_network_scores.py -a $NETWORK -r $REGULATORS -t $GENES -o ${OUTPUT_DIR}/motif_inference/network_scores/
python CODE/parse_quantized_bins.py -n 20 -i ${OUTPUT_DIR}/motif_inference/network_scores/ -o ${OUTPUT_DIR}/motif_inference/network_bins/

##Infer FIRE motifs
rm -f ${OUTPUT_DIR}/motif_inference/motif_inference.log
touch ${OUTPUT_DIR}/motif_inference/motif_inference.log

for f in ${TF_LIST}/*; do
	TF_SUBLIST=$f
	sbatch CODE/infer_motifs.sh $TF_SUBLIST $PROMOTER ${OUTPUT_DIR}/motif_inference/network_bins/ ${OUTPUT_DIR}/motif_inference/motif_inference.log
done

##Check if all motifs are ready
bash CODE/check_inference_status.sh ${OUTPUT_DIR}/motif_inference/motif_inference.log $REGULATORS $FLAG
