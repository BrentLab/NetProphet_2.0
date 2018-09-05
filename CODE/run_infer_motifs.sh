#!/bin/bash
function infer_motifs {
	FN_REGULATORS=$1	# a list of tf names
	FN_FASTA=$2 		# promoter sequence file
	DIR_BINNED_EXPR=$3 	# directory of binned expression files
	LOG_FILE=$4
	while read regulator; do
		if [[ ! -z ${regulator} ]]; then
		perl ${FIREDIR}/fire.pl --expfiles=${DIR_BINNED_EXPR}/$regulator --exptype=discrete --fastafile_dna=${FN_FASTA} --k=7 --jn=20 --jn_t=16 --nodups=1 --dorna=0 --dodnarna=0
		echo $regulator >> $LOG_FILE
		fi
	done < $FN_REGULATORS
}


OUTPUT_DIR=$1
NETWORK=$2
REGULATORS=$3
GENES=$4
PROMOTER=$5
FLAG=$6
USE_SERIAL=$7

## Prepare score bins
printf "Binning promoters based on network scores ... "
python CODE/parse_network_scores.py -a $NETWORK -r $REGULATORS -t $GENES -o ${OUTPUT_DIR}/motif_inference/network_scores/
python CODE/parse_quantized_bins.py -n 20 -i ${OUTPUT_DIR}/motif_inference/network_scores/ -o ${OUTPUT_DIR}/motif_inference/network_bins/
printf "DONE\n"

## Infer FIRE motifs using SLURM array scheme
rm -f ${OUTPUT_DIR}/motif_inference/motif_inference.log
touch ${OUTPUT_DIR}/motif_inference/motif_inference.log

printf "Inferring DNA binding motifs using FIRE ... "
if $USE_SERIAL; then
	infer_motifs $REGULATORS $PROMOTER ${OUTPUT_DIR}/motif_inference/network_bins/ ${OUTPUT_DIR}/motif_inference/motif_inference.log
else
	num_regulators=$( wc -l ${REGULATORS} | cut -d" " -f1 )
	sbatch --array=1-${num_regulators}%48 CODE/infer_motifs.sh $REGULATORS $PROMOTER ${OUTPUT_DIR}/motif_inference/network_bins/ ${OUTPUT_DIR}/motif_inference/motif_inference.log
fi

## Check if all motifs are ready
bash CODE/check_inference_status.sh ${OUTPUT_DIR}/motif_inference/motif_inference.log $REGULATORS $FLAG
