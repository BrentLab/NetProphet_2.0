#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -D ./LOG
#SBATCH -o motif_inference_%A.out
#SBATCH -e motif_inference_%A.err

fn_tf_list=$1		# a list of tf names
fn_fasta=$2 		# promoter sequence file
dir_binned_expr=$3 	# directory of binned expression files
log_file=$4


while read -a line
do
	tf=${line[0]}
	echo "__@__PROCESSING TF: $tf"
	perl ${FIREDIR}/fire.pl --expfiles=${dir_binned_expr}/$tf --exptype=discrete --fastafile_dna=${fn_fasta} --k=7 --jn=20 --jn_t=16 --nodups=1 --dorna=0 --dodnarna=0
	echo $tf >> ${log_file}
done < ${fn_tf_list}
