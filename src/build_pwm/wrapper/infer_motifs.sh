#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        -)
            case "${OPTARG}" in
                p_in_promoter)
                    p_in_promoter="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_expr_file)
                    p_in_expr_file="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
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
                p_src_code_fire)
                    p_src_code_fire="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
            esac;;
    esac
done

cmd=""
if [ ${flag_singularity} == "ON" ]; then
    p_src_code_fire=/home/packages_np3/FIRE-1.1a/
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_singularity.sh; fi
    export SINGULARITY_BINDPATH=${p_singularity_bindpath}
    cmd+="singularity exec ${p_singularity_img} "
elif [ ${flag_singularity} == "OFF" ]; then
    if [ ${flag_slurm} == "ON" ]; then source ${p_src_code}src/helper/load_modules.sh; fi
fi

cmd+="perl ${p_src_code_fire}fire.pl \
     --expfiles=${p_in_expr_file} \
     --exptype=discrete \
     --fastafile_dna=${p_in_promoter} \
     --k=7 \
     --jn=20 \
     --jn_t=16 \
     --nodups=1 \
     --dorna=0 \
     --dodnarna=0 >> /dev/null"
         
eval ${cmd} &> /dev/null
