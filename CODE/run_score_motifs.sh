#!/bin/bash 
OUTPUT_DIR=$1
MOTIFS_DIR=$2
REGULATORS=$3
PROMOTERS=$4
MOTIFS_LIST=$5
FLAG=$6

## Parse FIRE results
printf "Parsing motif inference results ... "
python CODE/parse_motif_summary.py -a True -i $MOTIFS_DIR -o $MOTIFS_LIST
python CODE/convert_fire2meme.py -i $MOTIFS_LIST -o $OUTPUT_DIR/motif_inference/motifs_pfm/
printf "DONE\n"

## Score FIRE motifs
rm -f ${OUTPUT_DIR}/motif_inference/motif_scoring.log
touch ${OUTPUT_DIR}/motif_inference/motif_scoring.log

printf "Score promoters using MEME-FIMO ... "
num_regulators=$( wc -l ${REGULATORS} )
sbatch --array=1-${num_regulators}%48 CODE/score_motifs.sh $REGULATORS $OUTPUT_DIR/motif_inference/motifs_pfm/ $PROMOTERS $OUTPUT_DIR/motif_inference/motifs_score ${OUTPUT_DIR}/motif_inference/motif_scoring.log

## Check if all motifs are ready
bash CODE/check_inference_status.sh ${OUTPUT_DIR}/motif_inference/motif_scoring.log $REGULATORS $FLAG