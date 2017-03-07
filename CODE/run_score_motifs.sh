#!/bin/bash 
OUTPUT_DIR=$1
MOTIFS_DIR=$2
REGULATORS=$3
TF_LIST=$4
PROMOTERS=$5
MOTIFS_LIST=$6
FLAG=$7

##Parse FIRE results
python CODE/parse_motif_summary.py -a True -i $MOTIFS_DIR -o $MOTIFS_LIST
python CODE/convert_fire2meme.py -i $MOTIFS_LIST -o $OUTPUT_DIR/motif_inference/motifs_pfm/

##Score scer FIRE motifs
rm -f ${OUTPUT_DIR}/../LOG/motif_scoring.log
touch ${OUTPUT_DIR}/../LOG/motif_scoring.log

for f in ${TF_LIST}/*; do
	TF_SUBLIST=$f;
	sbatch CODE/score_motifs.sh $TF_SUBLIST $OUTPUT_DIR/motif_inference/motifs_pfm/ $PROMOTERS $OUTPUT_DIR/motif_inference/motifs_score ${OUTPUT_DIR}/../LOG/motif_scoring.log;
done

##Check if all motifs are ready
bash CODE/check_inference_status.sh ${OUTPUT_DIR}/../LOG/motif_scoring.log $REGULATORS $FLAG