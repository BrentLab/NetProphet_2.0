#!/bin/bash

while getopts ":h-:" OPTION
do
    case "${OPTION}" in
        -)
            case "${OPTARG}" in
                # input
                p_in_net)
                    p_in_net="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                p_in_promoter)
                    p_in_promoter="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                
                # output
                p_out_dir)
                    p_out_dir="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                f_out_name)
                    f_out_name="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ;;
                
                # SLURM
                flag_slurm)
                    flag_slurm="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                
                # Singularity
                flag_singularity)
                    flag_singularity="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_singularity_img)
                    p_singularity_img="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                p_singularity_bindpath)
                    p_singularity_bindpath="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                
                # Logistics
                p_src_code)
                    p_src_code="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
                nbr_job)
                    nbr_job="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                    ;;
            esac;;
    esac
done

p_out_dir_tmp=${p_out_dir}tmp_pwm/
mkdir -p ${p_out_dir_tmp}


# =========================================================================== #
# |                           *** Infer Motifs ***                          | #
# =========================================================================== #
echo "- bin promoters based on network scores.."
mkdir -p ${p_out_dir_tmp}network_scores/
cmd_parse_network="${p_src_code}src/build_pwm/wrapper/parse_network_scores.sh \
    --p_in_net ${p_in_net} \
    --p_out_dir_scores ${p_out_dir_tmp}network_scores/ \
    --flag_slurm ${flag_slurm} \
    --flag_singularity ${flag_singularity} \
    --p_singularity_img ${p_singularity_img} \
    --p_singularity_bindpath ${p_singularity_bindpath} \
    --p_src_code ${p_src_code}"
eval ${cmd_parse_network}

echo "- parse quantized bins.."
mkdir -p ${p_out_dir_tmp}network_scores/
mkdir -p ${p_out_dir_tmp}network_bins/
cmd_parse_bins="${p_src_code}src/build_pwm/wrapper/parse_quantized_bins.sh \
    --p_in_dir_scores ${p_out_dir_tmp}network_scores/ \
    --p_out_dir_bins ${p_out_dir_tmp}network_bins/ \
    --flag_slurm ${flag_slurm} \
    --flag_singularity ${flag_singularity} \
    --p_singularity_img ${p_singularity_img} \
    --p_singularity_bindpath ${p_singularity_bindpath} \
    --p_src_code ${p_src_code}"
eval ${cmd_parse_bins}    
    
while read reg; do
    echo "   + reg: ${reg}.."
    if [ -f ${p_out_dir_tmp}network_bins/${reg} ]; then
        cmd_infer_motif=""
        cmd_infer_motif+="${p_src_code}src/build_pwm/wrapper/infer_motifs.sh \
                         --p_in_promoter ${p_in_promoter} \
                         --p_in_expr_file ${p_out_dir_tmp}network_bins/${reg} \
                         --flag_slurm ${flag_slurm} \
                         --flag_singularity ${flag_singularity} \
                         --p_singularity_img ${p_singularity_img} \
                         --p_singularity_bindpath ${p_singularity_bindpath} \
                         --p_src_code ${p_src_code} &"
       eval ${cmd_infer_motif}
       
       # manage the number of running jobs
       nbr_running_jobs=$(jobs -p | wc -l)
       while (( nbr_running_jobs >= nbr_job ))
       do
           sleep 1
           nbr_running_jobs=$(jobs -p | wc -l)
       done
    fi
done < <(cut -f1 ${p_in_net})
wait



# =========================================================================== #
# |                           *** Score Motifs ***                          | #
# =========================================================================== #

# define command: parse motif
echo "- parse motiff summary.."
cmd_parse_motif="${p_src_code}src/build_pwm/wrapper/parse_motif_summary.sh \
    --p_in_dir_bins ${p_out_dir_tmp}network_bins/ \
    --p_out_motifs_list ${p_out_dir_tmp}motifs.txt \
    --flag_slurm ${flag_slurm} \
    --flag_singularity ${flag_singularity} \
    --p_singularity_img ${p_singularity_img} \
    --p_singularity_bindpath ${p_singularity_bindpath} \
    --p_src_code ${p_src_code}"
# run command: parse motif 
eval ${cmd_parse_motif}

# define command: convert fire2meme
mkdir -p ${p_out_dir_tmp}motifs_pfm/
cmd_convert="${p_src_code}src/build_pwm/wrapper/convert_fire2meme.sh \
    --p_in_motifs_list ${p_out_dir_tmp}motifs.txt \
    --p_out_dir_pfm ${p_out_dir_tmp}motifs_pfm/ \
    --flag_slurm ${flag_slurm} \
    --flag_singularity ${flag_singularity} \
    --p_singularity_img ${p_singularity_img} \
    --p_singularity_bindpath ${p_singularity_bindpath} \
    --p_src_code ${p_src_code}"
# run command: convert fire2meme 
eval ${cmd_convert}

# loop over regualators and generate pwm network
if [ -d ${p_out_dir_tmp}motifs_scores/ ]; then
    rm -r ${p_out_dir_tmp}motifs_scores/
fi
mkdir ${p_out_dir_tmp}motifs_scores/
while read reg; do
    echo "   + reg: ${reg}.."
    if [ -f ${p_out_dir_tmp}motifs_pfm/${reg} ]; then
        cmd_score_motif="${p_src_code}src/build_pwm/wrapper/score_motifs.sh \
                         --p_in_pfm_reg ${p_out_dir_tmp}motifs_pfm/${reg} \
                         --p_in_promoter ${p_in_promoter} \
                         --p_out_score_reg ${p_out_dir_tmp}motifs_scores/${reg} \
                         --flag_slurm ${flag_slurm} \
                         --flag_singularity ${flag_singularity} \
                         --p_singularity_img ${p_singularity_img} \
                         --p_singularity_bindpath ${p_singularity_bindpath} \
                         --p_src_code ${p_src_code} &"
        eval ${cmd_score_motif}

        nbr_running_jobs=$(jobs -p | wc -l)
        while (( ${nbr_running_jobs} > ${nbr_job} ))
        do
            sleep 1
            nbr_running_jobs=$(jobs -p | wc -l)
        done
    fi
done < <(cut -f1 ${p_in_net})
wait

# =========================================================================== #
# |                          *** Build Network ***                          | #
# =========================================================================== #


cmd_build_net="${p_src_code}src/build_pwm/wrapper/build_net_motif.sh \
              --p_in_net ${p_in_net} \
              --p_in_motif_list ${p_out_dir_tmp}motifs.txt \
              --p_in_dir_score ${p_out_dir_tmp}motifs_scores/ \
              --p_out_pwm ${p_out_dir}${f_out_name} \
              --flag_slurm ${flag_slurm} \
              --flag_singularity ${flag_singularity} \
              --p_singularity_img ${p_singularity_img} \
              --p_singularity_bindpath ${p_singularity_bindpath} \
              --p_src_code ${p_src_code}"

eval ${cmd_build_net}              
