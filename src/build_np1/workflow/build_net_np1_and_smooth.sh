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
                flag_global_shrinkage)
                    flag_global_shrinkage="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                flag_local_shrinkage)
                    flag_local_shrinkage="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                seed)
                    seed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                lasso_nbr_fold)
                    lasso_nbr_fold="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                flag_microarray)
                    flag_microarray="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_net_de)
                    p_in_net_de="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_dbd_pids)
                    p_in_dbd_pids="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                    
                # Output
                p_out_dir)
                    p_out_dir="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_lasso)
                    f_out_name_lasso="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_np1)
                    f_out_name_np1="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name_np1_smoothed)
                    f_out_name_np1_smoothed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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



cmd_np1="${p_src_code}src/build_np1/wrapper/build_net_np1.sh \
        --p_in_expr_target ${p_in_expr_target} \
        --p_in_expr_reg ${p_in_expr_reg} \
        --flag_global_shrinkage ${flag_global_shrinkage} \
        --flag_local_shrinkage ${flag_local_shrinkage} \
        --lasso_nbr_fold ${lasso_nbr_fold} \
        --flag_microarray ${flag_microarray} \
        --seed ${seed}
        --p_in_net_de ${p_in_net_de} \
        --p_out_dir ${p_out_dir} \
        --f_out_name_lasso ${f_out_name_lasso} \
        --f_out_name_np1 ${f_out_name_np1} \
        --p_src_code ${p_src_code} \
        --flag_slurm ${flag_slurm} \
        --flag_singularity ${flag_singularity} \
        --p_singularity_img ${p_singularity_img} \
        --p_singularity_bindpath ${p_singularity_bindpath}"
    
eval ${cmd_np1}

if [ ${p_in_dbd_pids} != "NONE" ]
then
    cmd_smooth="${p_src_code}src/smooth_networks/wrapper/weighted_avg_similar_dbds_wrap.sh \
            --p_in_net ${p_out_dir}${f_out_name_np1} \
            --p_in_dbd_pids ${p_in_dbd_pids} \
            --p_out_dir ${p_out_dir} \
            --f_out_name_smoothed ${f_out_name_np1_smoothed} \
            --p_src_code ${p_src_code} \
            --flag_slurm ${flag_slurm} \
            --flag_singularity ${flag_singularity} \
            --p_singularity_img ${p_singularity_img} \
            --p_singularity_bindpath ${p_singularity_bindpath}"
            
    eval ${cmd_smooth}            
fi