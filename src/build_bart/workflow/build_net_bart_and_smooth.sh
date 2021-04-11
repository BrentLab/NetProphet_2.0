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
                p_in_expr_target)
                    p_in_expr_target="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_expr_reg)
                    p_in_expr_reg="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                bart_nbr_rmpi_slave)
                    nbr_rmpi_slave="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                seed)
                    seed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_dbd_pids)
                    p_in_dbd_pids="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # Output
                p_out_dir)
                    p_out_dir="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_bart)
                    f_out_name_bart="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_bart_smoothed)
                    f_out_name_bart_smoothed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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



cmd_bart="${p_src_code}src/build_bart/wrapper/build_net_bart.sh \
        ${p_in_expr_target} \
        ${p_in_expr_reg} \
        ${p_out_dir} \
        ${f_out_name_bart} \
        ${flag_slurm} \
        ${p_src_code} \
        ${flag_singularity} \
        ${p_singularity_img} \
        ${p_singularity_bindpath} \
        ${nbr_rmpi_slave}"
    
eval ${cmd_bart}

if [ ${p_in_dbd_pids} != "NONE" ]
then
    cmd_dbd="${p_src_code}src/smooth_networks/wrapper/weighted_avg_similar_dbds_wrap.sh \
            --p_in_net ${p_out_dir}${f_out_name_bart} \
            --p_in_dbd_pids ${p_in_dbd_pids} \
            --p_out_dir ${p_out_dir} \
            --f_out_name_smoothed ${f_out_name_bart_smoothed} \
            --p_src_code ${p_src_code} \
            --flag_slurm ${flag_slurm} \
            --flag_singularity ${flag_singularity} \
            --p_singularity_img ${p_singularity_img} \
            --p_singularity_bindpath ${p_singularity_bindpath}"
    
    eval ${cmd_dbd}            
fi
