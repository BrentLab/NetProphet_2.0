#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH -D ./LOG
#SBATCH -o motif_scoring_%A_%a.out
#SBATCH -e motif_scoring_%A_%a.err
#SBATCH -J score_motif

# Input variables
FN_REGULATORS=$1	# list of tf names
FN_TF_PWM=$2		# directory of tf pwm
FN_PROMOTERS=$3		# promoter sequence file (e.g. yeast_promoter_seq/s_cerevisiae.promoters.fasta, fly_promoter_seq/rsat_dmel_upstream_-2000_+200.filtered.fasta)
OUT_FIMO=$4			# directory of fimo alignment output 
LOG_FILE=$5

read regulator < <( sed -n ${SLURM_ARRAY_TASK_ID}p $FN_REGULATORS )
set -e

export MEME_BIN=/scratch/mblab/dabid/netprophet/code_netprophet2.0/SRC/meme/bin

if [[ ! -z ${regulator} ]]; then
	if [ -f $FN_TF_PWM/$regulator ]; then
		${MEME_BIN}/fimo -o $OUT_FIMO/$regulator --thresh 5e-3 $FN_TF_PWM/$regulator $FN_PROMOTERS
		sed ' 1d ' $OUT_FIMO/$regulator/fimo.txt | cut -f 1,2,7 > $OUT_FIMO/$regulator/temp.txt
		ruby /scratch/mblab/dabid/netprophet/code_netprophet3.0/code/netprophet2/estimate_affinity.rb -i $OUT_FIMO/$regulator/temp.txt > $OUT_FIMO/${regulator}.summary
	else
		printf "No motif inferred for %s\n" $regulator
	fi
	sleep $[ ($RANDOM % 20) + 1]s # sleep between 1-20 sec to prevent write lock
	echo $regulator >> $LOG_FILE
fi
