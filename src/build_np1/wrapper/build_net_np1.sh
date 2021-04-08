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
                p_in_net_de)
                    p_in_net_de="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                seed)
                    seed="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                flag_microarray)
                    flag_microarray="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                lasso_nbr_fold)
                    lasso_nbr_fold="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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

cmd=""
if [ ${flag_singularity} == "ON" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_modules.sh; fi
fi

cmd+="Rscript ${p_src_code}src/build_np1/code/build_net_np1.R \
          --p_in_expr_target ${p_in_expr_target} \
          --p_in_expr_reg ${p_in_expr_reg} \
          --flag_global_shrinkage ${flag_global_shrinkage} \
          --flag_local_shrinkage ${flag_local_shrinkage} \
          --lasso_nbr_fold ${lasso_nbr_fold} \
          --flag_microarray ${flag_microarray} \
          --p_in_net_de ${p_in_net_de} \
          --seed ${seed}
          --p_out_dir ${p_out_dir} \
          --f_out_name_lasso net_lasso.tsv \
          --f_out_name_np1 net_np1.tsv \
          --p_src_code ${p_src_code} \
          --flag_slurm ${flag_slurm}"
          
eval ${cmd}