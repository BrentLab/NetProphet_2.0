#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        -)
            case "${OPTARG}" in
                p_in_dir_scores)
                    p_in_dir_scores="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                p_out_dir_bins)
                    p_out_dir_bins="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                flag_slurm)
                    flag_slurm="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                flag_singularity)
                    flag_singularity="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_singularity_img)
                    p_singularity_img="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_singularity_bindpath)
                    p_singularity_bindpath="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_src_code)
                    p_src_code="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
            esac;;
    esac
done

cmd=""
if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "; fi
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ]; then
        source ${p_src_code}src/helper/load_modules.sh
        source activate np2
        ls -l ${CONDA_PREFIX}/bin > /dev/null
    fi
fi


cmd+="python3 ${p_src_code}src/build_pwm/code/netprophet2/parse_quantized_bins.py \
    -n 20 \
    -i ${p_in_dir_scores} \
    -o ${p_out_dir_bins}"
    
eval ${cmd}    

if [ ${flag_singularity} == "OFF" ] && [ ${flag_slurm} == "ON" ]; then
    source deactivate np2
fi