#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        h)
            exit 2
            ;;
        -)
            case "${OPTARG}" in
                # Input
                p_in_net)
                    p_in_net="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_dbd_pids)
                    p_in_dbd_pids="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                # Output
                p_out_dir)
                    p_out_dir="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_smoothed)
                    f_out_name_smoothed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                # logistics
                p_src_code)
                    p_src_code="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # slurm
                flag_slurm)
                    flag_slurm="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # singularity
                flag_singularity)
                    flag_singularity="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_singularity_img)
                    p_singularity_img="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_singularity_bindpath)
                    p_singularity_bindpath="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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
    if [ ${flag_slurm} == "ON" ] 
    then 
        source ${p_src_code}src/helper/load_modules.sh
        source activate np2
        ls -l ${SLURM_SUBMIT_DIR}np2/bin > /dev/null
    fi
    
fi

# get the list of regulators
p_in_reg=${p_out_dir}tmp_np1/reg.tsv
mkdir -p ${p_out_dir}tmp_np1/
tail -n+2 ${p_in_net} | cut -f1 > ${p_in_reg}  # -n+2 to ignore the header

cmd+="python3 ${p_src_code}src/smooth_networks/code/weighted_avg_similar_dbds.py \
         -n ${p_in_net} \
         -r ${p_in_reg} \
         -a ${p_in_dbd_pids} \
         -d 50 \
         -f single_dbds \
         -o ${p_out_dir}${f_out_name_smoothed}"

eval ${cmd}

if [ ${flag_singularity} == "OFF" ] && [ ${flag_slurm} == "ON" ]; then source deactivate np2; fi