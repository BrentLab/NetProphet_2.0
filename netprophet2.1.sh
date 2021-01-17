#!/bin/bash

# ========================================================================= #
# |                    **** USAGE OF NETPROPHET 2.1 ****                  | #      
# ========================================================================= #
usage(){
cat << EOF
    netprophet2.1 [options]
    
    INPUT arguments:
    --p_in_target            : input file for list of target genes
    --p_in_reg               : input file for list of regulators
    --p_in_sample            : input file for list of sample (condition) ids
    --p_in_expr_target       : input file for expression of target genes (target x sample)
    --p_in_net_de            : input file for network of differential expression
    --p_in_binding_event     : input file for binding events |REGULATOR|TARGET|
    
    NETPROPHET1 arguments:
    --fname_net_np1          : name of generated network for netprophet1 (optional)
    
    NETPROPHET2 arguments:
    --fname_net_bart         : name of generated network for bart (optional)
    --fname_net_np2          : name of generated network for netprophet2 (optional)
    --p_in_dbd_pids          : path of folder for DBD_PIDS
    --p_in_promoter          : path of file for promoters
    
    OUTPUT arguments:
    --p_out_dir              : path of output directory for results

    SLURM arguments:
    --p_out_logs             : output director for log files (.out & .err) for slurm runs
EOF
}


# ========================================================================= #
# |                  **** END USAGE OF NETPROPHET 2.1 ****                | #      
# ========================================================================= #

# initialize flags for specific components runs
flag_prep=0  # -p
flag_np1=0  # -n
flag_bart=0  # bart
flag_avg_np1=0  # average netprophet1 network
flag_avg_bart=0  # average bart network
flag_combine=0  # combine averaged networks
flag_infer_motif=0  # infer motifs
flag_score_motif=0  # score motifs
flag_build_motif_net=0  # build motif network
flag_combine_motif_net=0  # combine with motif network
flag_avg_final_net=0  # averate final network
flag_all=0  # to run all components



# ========================================================================= #
# |                         **** PARSE ARGUMENTS ****                     | #      
# ========================================================================= #

while getopts ":hpnbvwcismkya-:" OPTION
do
    case "${OPTION}" in
    p)
        flag_prep=1
        ;;
    n)
        flag_np1=1
        ;;
    b)
        flag_bart=1
        ;;
    v)
        flag_avg_np1=1
        ;;
    w)
        flag_avg_bart=1
        ;;
    c)
        flag_combine=1
        ;;
    i)
        flag_infer_motif=1
        ;;
    s)
        flag_score_motif=1
        ;;
    m)
        flag_build_motif_net=1
        ;;
    k)
        flag_combine_motif_net=1
        ;;
    y)
        flag_avg_final_net=1
        ;;
    a)
        flag_all=1
        ;;
    h)
        usage
        exit 2
        ;;
    -)
        case "${OPTARG}" in
            p_in_target)
                p_in_target="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_reg)
                p_in_reg="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_sample)
                p_in_sample="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_expr_target)
                p_in_expr_target="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_net_de)
                p_in_net_de="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_binding_event)
                p_in_binding_event="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            fname_net_bart)
                fname_net_bart="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            fname_net_np1)
                fname_net_np1="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            fname_net_np2)
                fname_net_np2="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_dbd_pids)
                p_in_dbd_pids="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_in_promoter)
                p_in_promoter="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_out_dir)
                p_out_dir="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            p_out_logs)
                p_out_logs="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
            data)
                data="${!OPTIND}"; OPTIND=$(( ${OPTIND} + 1 ))
                ;;
        esac;;
    esac

done

# ========================================================================= #
# |                       **** END PARSE ARGUMENTS ****                   | #      
# ========================================================================= #

# general parameters
p_wd=/scratch/mblab/dabid/netprophet/
p_src_code=${p_wd}code_netprophet2.1/

# OUTPUT parameters
p_out_tmp=${p_out_dir}tmp/
p_out_net=${p_out_dir}net/

# NETPROPHET1 parameters
fname_net_np1=net_np1.tsv
p_in_expr_reg=${p_out_tmp}expr_reg
p_in_fc=${p_out_tmp}fc
p_in_allowed=${p_out_tmp}allowed
p_in_pert_adj=${p_out_tmp}pert.adj
p_in_pert_tsv=${p_out_tmp}pert.tsv


# NETPROHET2 parameters
fname_net_bart=net_bart
fname_net_np1wa=net_np1wa.tsv
fname_net_bwa=net_bartwa.tsv
fname_net_np1wa_bwa=net_np1wa_bwa.tsv
fname_net_motif=net_motif.tsv
fname_net_np2=net_np2.tsv
fname_net_np2wa=net_np2wa.tsv

# initialize job ids for slurm runs
job_id_prepare=1
job_id_np1=1
job_id_bart=1
job_id_avg_np1=1
job_id_avg_bart=1
job_id_combine=1
job_id_infer_motif=1
job_id_score_motif=1
job_id_build_motif_net=1
job_id_build_final_net=1


# create output directories
mkdir -p ${p_out_tmp}
mkdir -p ${p_out_net}
mkdir -p ${p_out_tmp}motif_inference/
mkdir -p ${p_out_tmp}motif_inference/network_bins/
mkdir -p ${p_out_tmp}motif_inference/motifs_pfm/
mkdir -p ${p_out_tmp}motif_inference/motifs_score/

# ==================================================== #
# |            **** PREPARE RESOURCES ***            | #
# ==================================================== #
if [[ ${flag_prep}  == 1 || ${flag_all} == 1 ]]
then
    job_prepare=$(sbatch \
    -o ${p_out_logs}${data}_prepare_resources_%A.out \
    -e ${p_out_logs}${data}_prepare_resources_%A.err \
    -J ${data}_prepare_resources \
    ${p_src_code}CODE/prepare_resources_wrap.sh \
    ${p_in_target} \
    ${p_in_reg} \
    ${p_in_expr_target} \
    ${p_in_sample} \
    ${p_in_expr_reg} \
    ${p_in_fc} \
    ${p_in_allowed} \
    ${p_in_pert_adj} \
    ${p_in_pert_tsv})

    job_id_prepare=$(echo ${job_prepare} | awk '{split($0, a, " "); print a[4]}')
fi

# ========================================================== #
# |                    **** NETPROPHET1 ***                | #
# ========================================================== #

# --mem-per-cpu=40G \
# --cpus-per-task=2 \
    
if [[ ${flag_np1} == 1 || ${flag_all} == 1 ]]
then
    job_np1=$(sbatch \
    -o ${p_out_logs}${data}_np1_%A.out \
    -e ${p_out_logs}${data}_np1_%A.err \
    -n 11 \
    --dependency=afterany:${job_id_prepare} \
    --cpus-per-task=2 \
    --mem-per-cpu=10G \
    -D ${p_src_code}SRC/NetProphet1/ \
    -J ${data}_np1 \
    ${p_src_code}SRC/NetProphet1/netprophet \
       -m \
       -c \
       -u ${p_src_code}SRC/NetProphet1/ \
       -t ${p_in_expr_target} \
       -r ${p_in_expr_reg} \
       -a ${p_in_allowed} \
       -p ${p_in_pert_adj} \
       -d ${p_in_net_de} \
       -g ${p_in_target} \
       -f ${p_in_reg} \
       -o ${p_out_net} \
       -n ${fname_net_np1})

    job_id_np1=$(echo ${job_np1} | awk '{split($0, a, " "); print a[4]}')
fi

       
# ==================================================== #
# |                    **** BART ***                 | #
# ==================================================== #
if [[ ${flag_bart} == 1 || ${flag_all} == 1 ]]
then
    job_bart=$(sbatch \
    -o ${p_out_logs}${data}_bart_%A.out \
    -e ${p_out_logs}${data}_bart_%A.err \
    -n 32 \
    -J ${data}_bart \
    --mem=20GB \
    --dependency=afterany:${job_id_prepare} \
    ${p_src_code}CODE/run_build_bart_network.sh \
    ${p_in_fc} \
    ${p_in_pert_tsv} \
    ${p_in_reg} \
    ${p_out_net}${fname_net_bart} \
    false)

    job_id_bart=$(echo ${job_bart} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |         **** Weighted average np1 ***            | #
# ==================================================== #
if [[ ${flag_avg_np1} == 1 || ${flag_all} == 1 ]]
then
    job_avg_np1=$(sbatch \
    -o ${p_out_logs}${data}_avg_np1_%A.out \
    -e ${p_out_logs}${data}_avg_np1_%A.err \
    -J ${data}_avg_np1 \
    --dependency=afterany:${job_id_np1} \
    ${p_src_code}CODE/weighted_avg_similar_dbds_wrap.sh \
    ${p_out_net}${fname_net_np1} \
    ${p_in_reg} \
    ${p_in_dbd_pids} \
    50 \
    single_dbds \
    ${p_out_net}${fname_net_np1wa} \
    ${p_src_code})

    job_id_avg_np1=$(echo ${job_avg_np1} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |        **** Weighted average bart ***            | #
# ==================================================== #
if [[ ${flag_avg_bart} == 1 || ${flag_all} == 1 ]]
then
    job_avg_bart=$(sbatch \
    -o ${p_out_logs}${data}_avg_bart_%A.out \
    -e ${p_out_logs}${data}_avg_bart_%A.err \
    -J ${data}_avg_bart \
    --dependency=afterany:${job_id_bart} \
    ${p_src_code}CODE/weighted_avg_similar_dbds_wrap.sh \
    ${p_out_net}${fname_net_bart} \
    ${p_in_reg} \
    ${p_in_dbd_pids} \
    50 \
    single_dbds \
    ${p_out_net}${fname_net_bwa} \
    ${p_src_code})

    job_id_avg_bart=$(echo ${job_avg_bart} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |           **** Combine np & bart ***             | #
# ==================================================== #
if [[ ${flag_combine} == 1 || ${flag_all} == 1 ]]
then
    job_combine=$(sbatch \
    -o ${p_out_logs}${data}_combine_np_bart_%A.out \
    -e ${p_out_logs}${data}_combine_np_bart_%A.err \
    -J ${data}_combine_np_bart \
    --dependency=afterany:${job_id_avg_np1}:${job_id_avg_bart} \
    ${p_src_code}CODE/quantile_combine_networks_wrap.sh \
    ${p_out_net}${fname_net_np1wa} \
    ${p_out_net}${fname_net_bwa} \
    ${p_out_net}${fname_net_np1wa_bwa})

    job_id_combine=$(echo ${job_combine} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |              **** Infer Motif ***                | #
# ==================================================== #

num_regulators=$(wc -l ${p_in_reg} | cut -d" " -f1)

if [[ ${flag_infer_motif} == 1 || ${flag_all} == 1 ]]
then
    job_infer_motif=$(sbatch \
    -o ${p_out_logs}${data}_infer_motif_%A_%a.out \
    -e ${p_out_logs}${data}_infer_motif_%A_%a.err \
    -J ${data}_infer_motif \
    --array=1-${num_regulators}%48 \
    --dependency=afterany:${job_id_combine} \
    ${p_src_code}CODE/run_infer_motifs.sh \
    ${p_out_tmp} \
    ${p_out_net}${fname_net_np1wa_bwa} \
    ${p_in_reg} \
    ${p_in_target} \
    ${p_in_promoter} \
    ${p_out_tmp}flag_infer_motifs \
    false)

    job_id_infer_motif=$(echo ${job_infer_motif} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |              **** Score Motif ***                | #
# ==================================================== #
if [[ ${flag_score_motif} == 1 || ${flag_all} == 1 ]]
then
    job_score_motif=$(sbatch \
    -o ${p_out_logs}${data}_score_motif_%A_%a.out \
    -e ${p_out_logs}${data}_score_motif_%A_%a.err \
    -J ${data}_score_motif \
    --array=1-${num_regulators}%48 \
    --dependency=afterany:${job_id_infer_motif} \
    ${p_src_code}CODE/run_score_motifs.sh \
    ${p_out_tmp} \
    ${p_out_tmp}motif_inference/network_bins/ \
    ${p_in_reg} \
    ${p_in_promoter} \
    ${p_out_tmp}motif_inference/motifs.txt \
    ${p_out_tmp}flag_score_motifs \
    false)

    job_id_score_motif=$(echo ${job_score_motif} | awk '{split($0, a, " "); print a[4]}')
fi
# ==================================================== #
# |           **** Build Motif Network ***           | #
# ==================================================== #
if [[ ${flag_build_motif_net} == 1 || ${flag_all} == 1 ]]
then
    job_build_motif_net=$(sbatch \
    -o ${p_out_logs}${data}_build_motif_net_%A.out \
    -e ${p_out_logs}${data}_build_motif_net_%A.err \
    -J ${data}_build_motif_net \
    --dependency=afterany:${job_id_score_motif} \
    ${p_src_code}CODE/build_motif_network_wrap.sh \
    ${p_out_tmp}motif_inference/motifs.txt \
    ${p_in_reg} \
    ${p_in_target} \
    ${p_out_tmp}motif_inference/motifs_score/ \
    robust \
    16 \
    ${p_out_net}${fname_net_motif})

    job_id_build_motif_net=$(echo ${job_build_motif_net} | awk '{split($0, a, " "); print a[4]}')
fi


# ==================================================== #
# |         **** Assemble final network ***          | #
# ==================================================== #
if [[ ${flag_combine_motif_net} == 1 || ${flag_all} == 1 ]]
then
    job_build_final_net=$(sbatch \
    -o ${p_out_logs}${data}_build_final_net_%A.out \
    -e ${p_out_logs}${data}_build_final_net_%A.err \
    -J ${data}_build_final_net \
    --dependency=afterany:${job_id_build_motif_net} \
    ${p_src_code}CODE/combine_networks_wrap.sh \
    resort \
    ${p_out_net}${fname_net_np1wa_bwa} \
    ${p_out_net}${fname_net_motif} \
    ${p_out_net} \
    ${fname_net_np2} \
    ${p_src_code})

    job_id_build_final_net=$(echo ${job_build_final_net} | awk '{split($0, a, " "); print a[4]}')
fi

# ===================================================== #
# |        **** avg Assemble final network ***        | #
# ===================================================== #
if [[ ${flag_avg_final_net} == 1 || ${flag_all} == 1 ]]
then
    job_avg_final_net=$(sbatch \
    -o ${p_out_logs}${data}_avg_final_net_%A.out \
    -e ${p_out_logs}${data}_avg_final_net_%A.err \
    -J ${data}_avg_final_net \
    --dependency=afterany:${job_id_build_final_net} \
    ${p_src_code}CODE/weighted_avg_similar_dbds_wrap.sh \
    ${p_out_net}${fname_net_np2} \
    ${p_in_reg} \
    ${p_in_dbd_pids} \
    50 \
    single_dbds \
    ${p_out_net}${fname_net_np2wa} \
    ${p_src_code})
fi