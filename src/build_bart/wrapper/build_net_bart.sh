#!/bin/bash

# =========================================================================== #
# |                         **** PARSE ARGUMENTS ****                       | #
# =========================================================================== #
p_in_expr_target=${1}
p_in_expr_reg=${2}
p_out_dir=${3}
fname_bart=${4}
flag_slurm=${5}
p_src_code=${6}
flag_singularity=${7}
p_singularity_img=${8}
p_singularity_bindpath=${9}


# =========================================================================== #
# |                        **** BUILD BART ****                             | #
# =========================================================================== #
echo "build BART network.."
cmd=""
if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_modules.sh; fi
fi

cmd+="Rscript --no-save --vanilla ${p_src_code}src/build_bart/code/build_net_bart.R \
     --p_in_expr_target ${p_in_expr_target} \
     --p_in_expr_reg ${p_in_expr_reg} \
     --fname_bart ${fname_bart} \
     --p_out_dir ${p_out_dir} \
     --flag_slurm ${flag_slurm} \
     --p_src_code ${p_src_code}"

eval ${cmd}     