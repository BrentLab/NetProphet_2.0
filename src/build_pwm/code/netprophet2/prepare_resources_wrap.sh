#!/bin/bash

#SBATCH --mem-per-cpu=10G

# get the parameters
f_genes=${1}
f_regulators=${2}
f_data_expr=${3}
f_conditions=${4}
f_tmp_rdata_expr=${5}
f_tmp_data_fc_tsv=${6}
f_tmp_allowed_adj=${7}
f_tmp_data_pert_adj=${8}
f_tmp_data_pert_tsv=${9}

# load modules
module load anaconda3/4.1.1
source activate netprophet

# run prepare data
python /scratch/mblab/dabid/netprophet/code_netprophet3.0/CODE/prepare_resources.py \
-g ${f_genes} \
-r ${f_regulators} \
-e ${f_data_expr} \
-c ${f_conditions} \
-or ${f_tmp_rdata_expr} \
-of ${f_tmp_data_fc_tsv} \
-oa ${f_tmp_allowed_adj} \
-op1 ${f_tmp_data_pert_adj} \
-op2 ${f_tmp_data_pert_tsv}
