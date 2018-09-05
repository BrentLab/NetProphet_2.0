#!/bin/bash 
function score_motifs {
	FN_REGULATORS=$1	# list of tf names
	FN_TF_PWM=$2		# directory of tf pwm
	FN_PROMOTERS=$3		# promoter sequence file
	OUT_FIMO=$4			# directory of fimo alignment output 
	LOG_FILE=$5
	while read regulator; do
		if [[ ! -z ${regulator} ]]; then
			if [ -f $FN_TF_PWM/$regulator ]; then
				fimo -o $OUT_FIMO/$regulator --thresh 5e-3 $FN_TF_PWM/$regulator $FN_PROMOTERS
				sed ' 1d ' $OUT_FIMO/$regulator/fimo.txt | cut -f 1,2,7 > $OUT_FIMO/$regulator/temp.txt
				ruby ./CODE/estimate_affinity.rb -i $OUT_FIMO/$regulator/temp.txt > $OUT_FIMO/${regulator}.summary
			else
				printf "No motif inferred for %s\n" $regulator
			fi
			echo $regulator >> $LOG_FILE
		fi
	done < $FN_REGULATORS
}


OUTPUT_DIR=$1
MOTIFS_DIR=$2
REGULATORS=$3
PROMOTERS=$4
MOTIFS_LIST=$5
FLAG=$6
USE_SERIAL=$7

## Parse FIRE results
printf "Parsing motif inference results ... "
python CODE/parse_motif_summary.py -a True -i $MOTIFS_DIR -o $MOTIFS_LIST
python CODE/convert_fire2meme.py -i $MOTIFS_LIST -o $OUTPUT_DIR/motif_inference/motifs_pfm/
printf "DONE\n"

## Score FIRE motifs
rm -f ${OUTPUT_DIR}/motif_inference/motif_scoring.log
touch ${OUTPUT_DIR}/motif_inference/motif_scoring.log

printf "Score promoters using MEME-FIMO ... "
if $USE_SERIAL; then
	score_motifs $REGULATORS $OUTPUT_DIR/motif_inference/motifs_pfm/ $PROMOTERS $OUTPUT_DIR/motif_inference/motifs_score ${OUTPUT_DIR}/motif_inference/motif_scoring.log
else
	num_regulators=$( wc -l ${REGULATORS} | cut -d" " -f1 )
	sbatch --array=1-${num_regulators}%48 CODE/score_motifs.sh $REGULATORS $OUTPUT_DIR/motif_inference/motifs_pfm/ $PROMOTERS $OUTPUT_DIR/motif_inference/motifs_score ${OUTPUT_DIR}/motif_inference/motif_scoring.log
fi

## Check if all motifs are ready
bash CODE/check_inference_status.sh ${OUTPUT_DIR}/motif_inference/motif_scoring.log $REGULATORS $FLAG