#!/bin/bash

module load anaconda3/4.1.1

resort=${1}
npwa_nbwa=${2}
nm=${3}
p_out_net=${4}
nbwa_bnwa_nm=${5}
p_src_code=${6}
source activate analysis

python /scratch/mblab/dabid/netprophet/code_netprophet2.1/CODE/combine_networks.py \
-s resort \
-n ${npwa_nbwa} \
-b ${nm} \
-od ${p_out_net} \
-om ${nbwa_bnwa_nm}

source deactivate analysis