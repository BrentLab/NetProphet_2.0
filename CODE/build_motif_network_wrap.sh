#!/bin/bash

module load anaconda3/4.1.1

p_motif=${1}
p_in_reg=${2}
p_in_target=${3}
p_motifs_score=${4}
t=${5}
v=${6}
p_mn=${7}  # output file

source activate np2
python /scratch/mblab/dabid/netprophet/code_netprophet2.1/CODE/build_motif_network.py \
-i ${p_motif} \
-r ${p_in_reg} \
-g ${p_in_target} \
-f ${p_motifs_score} \
-t ${t} \
-v ${v} \
-o ${p_mn}

source deactivate np2