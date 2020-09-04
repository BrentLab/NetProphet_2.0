#!/bin/bash

module load anaconda3/4.1.1
source activate analysis

p_in_net=${1}
p_in_reg=${2}
p_dbd_pid=${3}
d=${4}
f=${5}
p_out_net=${6}
p_src_code=${7}

python ${p_src_code}CODE/weighted_avg_similar_dbds.py \
-n ${p_in_net} \
-r ${p_in_reg} \
-a ${p_dbd_pid} \
-d ${d} \
-f ${f} \
-o ${p_out_net}

source deactivate analysis