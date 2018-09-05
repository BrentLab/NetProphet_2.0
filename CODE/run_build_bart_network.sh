#!/bin/bash
#SBATCH -n 32
#SBATCH --mem=20G
#SBATCH -D ./LOG
#SBATCH -o map_bart_network_%A.out
#SBATCH -e map_bart_network_%A.err

data_fc_expr=${1}
pert_matrix=${2}
tf_names=${3}
output_adjmtr=${4}
use_serial=${5}

if $use_serial ; then
	Rscript --vanilla ./CODE/build_bart_network.r fcFile=${data_fc_expr} isPerturbedFile=${pert_matrix} tfNameFile=${tf_names} saveTo=${output_adjmtr}.tsv useMpi=FALSE
else 
	module load R/3.2.1
	module load openmpi
	Rscript --vanilla ../CODE/build_bart_network.r fcFile=${data_fc_expr} isPerturbedFile=${pert_matrix} tfNameFile=${tf_names} saveTo=${output_adjmtr}.tsv useMpi=TRUE mpiBlockSize=32
fi
sed '1d' ${output_adjmtr}.tsv > ${output_adjmtr}
awk -i inplace '{sub(/^\S+\s*/,"")}1' ${output_adjmtr}
