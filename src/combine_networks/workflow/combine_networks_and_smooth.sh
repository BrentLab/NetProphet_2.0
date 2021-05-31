#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        -)
            case "${OPTARG}" in
                # Input
                p_in_net_np1_bart)
                    p_in_net_np1_bart="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_net_pwm)
                    p_in_net_pwm="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_dbd_pids)
                    p_in_dbd_pids="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                
                # Output
                p_out_dir)
                    p_out_dir="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_np2)
                    f_out_name_np2="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_np2_smoothed)
                    f_out_name_np2_smoothed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # Logistics
                p_src_code)
                    p_src_code="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # Slurm
                flag_slurm)
                    flag_slurm="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # Singularity
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
                 
cmd_combine="${p_src_code}src/combine_networks/wrapper/combine_networks.sh \
                 ${p_in_net_np1_bart} \
                 ${p_in_net_pwm} \
                 ${p_out_dir} \
                 ${f_out_name_np2} \
                 ${p_src_code} \
                 ${flag_slurm} \
                 ${flag_singularity} \
                 ${p_singularity_img} \
                 ${p_singularity_bindpath}"
                 
eval ${cmd_combine}

if [ ${p_in_dbd_pids} != "NONE" ]; then
    cmd_smooth="${p_src_code}src/smooth_networks/wrapper/weighted_avg_similar_dbds_wrap.sh \
                    --p_in_net ${p_out_dir}${f_out_name_np2} \
                    --p_in_dbd_pids ${p_in_dbd_pids} \
                    --p_out_dir ${p_out_dir} \
                    --f_out_name_smoothed ${f_out_name_np2_smoothed} \
                    --p_src_code ${p_src_code} \
                    --flag_slurm ${flag_slurm} \
                    --flag_singularity ${flag_singularity} \
                    --p_singularity_img ${p_singularity_img} \
                    --p_singularity_bindpath ${p_singularity_bindpath}"
    
    eval ${cmd_smooth}
fi