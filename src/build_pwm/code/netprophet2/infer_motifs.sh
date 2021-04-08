#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH -D ./LOG
#SBATCH -o motif_inference_%A_%a.out
#SBATCH -e motif_inference_%A_%a.err
#SBATCH -J infer_motifs

FN_REGULATORS=$1	# a list of tf names
FN_FASTA=$2 		# promoter sequence file
DIR_BINNED_EXPR=$3 	# directory of binned expression files
LOG_FILE=$4

read regulator < <( sed -n ${SLURM_ARRAY_TASK_ID}p $FN_REGULATORS )
set -e

if [[ ! -z ${regulator} ]]; then
    
	perl /scratch/mblab/dabid/netprophet/code_netprophet2.0/SRC/FIRE-1.1a/fire.pl \
        --expfiles=${DIR_BINNED_EXPR}/$regulator \
        --exptype=discrete \
        --fastafile_dna=${FN_FASTA} \
        --k=7 \
        --jn=20 \
        --jn_t=16 \
        --nodups=1 \
        --dorna=0 \
        --dodnarna=0
	echo $regulator >> $LOG_FILE
fi
