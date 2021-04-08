#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        -)
            case "${OPTARG}" in
                p_in_net)
                    p_in_net="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_in_motif_list)
                    p_in_motif_list="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_in_dir_score)
                    p_in_dir_score="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_out_pwm)
                    p_out_pwm="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                flag_slurm)
                    flag_slurm="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                flag_singularity)
                    flag_singularity="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_singularity_img)
                    p_singularity_img="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_singularity_bindpath)
                    p_singularity_bindpath="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_src_code)
                    p_src_code="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
            esac;;
    esac
done


cmd=""
if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ]; then
        source ${p_src_code}src/helper/load_modules.sh
        source activate np3
        ls -l ${SLURM_SUBMIT_DIR}np3/bin > /dev/null
    fi
fi


cmd+="python3 ${p_src_code}src/build_pwm/code/netprophet2/build_motif_network.py \
    -i ${p_in_motif_list} \
    -n ${p_in_net} \
    -f ${p_in_dir_score} \
    -t robust \
    -v 16 \
    -o ${p_out_pwm}"
    
eval ${cmd}

if [ ${flag_singularity} == "OFF" ] && [ ${flag_slurm} == "ON" ]; then source deactivate np3; fi 