#!/bin/bash

p_npwa=${1}
p_bnwa=${2}
p_npwa_bnwa=${3}
p_src_code=${4}
flag_slurm=${5}
flag_singularity=${6}
p_singularity_img=${7}
p_singularity_bindpath=${8}

cmd=""

if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ] 
    then 
        source ${p_src_code}src/helper/load_modules.sh
    fi
    
fi

cmd+="Rscript ${p_src_code}src/combine_networks/code/quantile_combine_networks.r \
     ${p_npwa} \
     ${p_bnwa} \
     ${p_npwa_bnwa}"

eval ${cmd}