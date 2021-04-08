#!/bin/bash

p_net_np1_bart=${1}
p_net_pwm=${2}
p_out_dir=${3}
f_out_name_np2=${4}
p_src_code=${5}
flag_slurm=${6}
flag_singularity=${7}
p_singularity_img=${8}
p_singularity_bindpath=${9}


cmd=""
if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ] 
    then 
        source ${p_src_code}src/helper/load_modules.sh
        source activate np2
        ls -l ${SLURM_SUBMIT_DIR}np2/bin > /dev/null
    fi
    
fi

cmd+="python3 ${p_src_code}src/combine_networks/code/combine_networks.py \
          -s resort \
          -n ${p_net_np1_bart} \
          -b ${p_net_pwm} \
          -od ${p_out_dir} \
          -om ${f_out_name_np2}"
eval ${cmd}


if [ ${flag_singularity} == "OFF" ] && [ ${flag_slurm} == "ON" ]; then source deactivate np2; fi