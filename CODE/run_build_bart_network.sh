#!/bin/bash

data_fc_expr=${1}
pert_matrix=${2}
tf_names=${3}
output_adjmtr=${4}
use_serial=${5}

module load R/3.4.3 

if $use_serial ; then
	Rscript --vanilla /scratch/mblab/dabid/netprophet/code_netprophet2.1/CODE/build_bart_network.r fcFile=${data_fc_expr} isPerturbedFile=${pert_matrix} tfNameFile=${tf_names} saveTo=${output_adjmtr}.tsv useMpi=FALSE
else 

	module load openmpi/1.8.8
    Rscript --vanilla /scratch/mblab/dabid/netprophet/code_netprophet2.1/CODE/build_bart_network.r fcFile=${data_fc_expr} isPerturbedFile=${pert_matrix} tfNameFile=${tf_names} saveTo=${output_adjmtr}.tsv useMpi=TRUE mpiBlockSize=32
fi
sed '1d' ${output_adjmtr}.tsv > ${output_adjmtr}
awk -i inplace '{sub(/^\S+\s*/,"")}1' ${output_adjmtr}
