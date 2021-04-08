#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        -)
            case "${OPTARG}" in
                p_in_pfm_reg)
                    p_in_pfm_reg="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_promoter)
                    p_in_promoter="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_out_score_reg)
                    p_out_score_reg="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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



cmd_fimo=""
cmd_estimate=""
if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd_fimo+="singularity exec ${p_singularity_img} "
    cmd_estimate+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ]; then
        source ${p_src_code}src/helper/load_modules.sh
    fi
fi

cmd_fimo+="fimo \
      --o ${p_out_score_reg} \
      --thresh 5e-3 \
      --verbosity 1 \
      ${p_in_pfm_reg} \
      ${p_in_promoter}"


cmd_prepare_data="sed '1d' ${p_out_score_reg}/fimo.tsv | cut -f 1,3,8 > ${p_out_score_reg}/temp.txt"

cmd_estimate+="ruby ${p_src_code}src/build_pwm/code/netprophet2/estimate_affinity.rb \
              -i ${p_out_score_reg}/temp.txt > ${p_out_score_reg}.summary"

eval ${cmd_fimo}
eval ${cmd_prepare_data}
eval ${cmd_estimate} &> /dev/null