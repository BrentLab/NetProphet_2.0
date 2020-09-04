#!/bin/bash

p_wd=/scratch/mblab/dabid/netprophet/

# input variables
data=zev
p_src_code=${p_wd}code_netprophet2.1/
p_in_reg=${p_wd}net_in/sub_in_reg_tf_151
p_in_target=${p_wd}net_in/sub_in_target_6023
p_in_sample=${p_wd}net_in/zev_in_condition_591
p_in_expr_target=${p_wd}net_in/sub2_zev_expr_6023_591
p_in_de=${p_wd}net_in/sub2_zev_de_shrunken_151_6023
p_dbd_pid=${p_wd}net_in/DBD_PIDS
p_in_promoter=${p_wd}net_in/promoter.scer.fasta

# output variables
p_out_dir=${p_wd}net_out/zev_netprophet2/
p_out_tmp=${p_out_dir}tmp/
p_out_net=${p_out_dir}networks/

p_in_expr_reg=${p_out_tmp}zev_expr_reg_151_1485
p_in_fc=${p_out_tmp}zev_fc
p_in_allowed=${p_out_tmp}allowed
p_in_pert_adj=${p_out_tmp}pert.adj
p_in_pert_tsv=${p_out_tmp}pert.tsv
f_out_net_lasso=net_np1
f_out_net_bart=net_bart
f_out_net_netprophet2=net_np2

# job ids for slurm
job_id_prepare=1
job_id_lasso=1
job_id_bart=1
job_id_avg_np=1
job_id_avg_bart=1
job_id_combine=1
job_id_infer_motif=1
job_id_score_motif=1
job_id_build_motif_net=1
job_id_build_final_net=1

# flag execution
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

while getopts "pnbvwcismkya" OPTION
do
    case ${OPTION} in
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
    esac
done

mkdir -p ${p_out_dir}
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
    -o ${p_wd}net_logs/${data}_prep_resources_%A.out \
    -e ${p_wd}net_logs/${data}_prep_resources_%A.err \
    -J ${data}_prepare_resource \
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

# ==================================================== #
# |                    **** LASSO ***                | #
# ==================================================== #

if [[ ${flag_np1} == 1 || ${flag_all} == 1 ]]
then
    job_lasso=$(sbatch \
    -o ${p_wd}net_logs/${data}_lasso_%A.out \
    -e ${p_wd}net_logs/${data}_lasso_%A.err \
    -n 11 \
    --mem-per-cpu=40G \
    --cpus-per-task=2 \
    --dependency=afterany:${job_id_prepare} \
    -D ${p_src_code}SRC/NetProphet1/ \
    -J ${data}_lasso \
    ${p_src_code}SRC/NetProphet1/netprophet \
       -m \
       -c \
       -u ${p_src_code}SRC/NetProphet1/ \
       -t ${p_in_expr_target} \
       -r ${p_in_expr_reg} \
       -a ${p_in_allowed} \
       -p ${p_in_pert_adj} \
       -d ${p_in_de} \
       -g ${p_in_target} \
       -f ${p_in_reg} \
       -o ${p_out_net} \
       -n ${f_out_net_lasso})

    job_id_lasso=$(echo ${job_lasso} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |                    **** BART ***                 | #
# ==================================================== #
if [[ ${flag_bart} == 1 || ${flag_all} == 1 ]]
then
    job_bart=$(sbatch \
    -o ${p_wd}net_logs/${data}_bart_%A.out \
    -e ${p_wd}net_logs/${data}_bart_%A.err \
    -n 32 \
    -J ${data}_bart \
    --mem=20GB \
    --dependency=afterany:${job_id_prepare} \
    ${p_src_code}CODE/run_build_bart_network.sh \
    ${p_in_fc} \
    ${p_in_pert_tsv} \
    ${p_in_reg} \
    ${p_out_net}${f_out_net_bart} \
    false)

    job_id_bart=$(echo ${job_bart} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |         **** Weighted average np ***             | #
# ==================================================== #
if [[ ${flag_avg_np1} == 1 || ${flag_all} == 1 ]]
then
    job_avg_np=$(sbatch \
    -o ${p_wd}net_logs/${data}_avg_np_%A.out \
    -e ${p_wd}net_logs/${data}_avg_np_%A.err \
    -J ${data}_avg_np \
    --dependency=afterany:${job_id_lasso} \
    ${p_src_code}CODE/weighted_avg_similar_dbds_wrap.sh \
    ${p_out_net}${f_out_net_lasso} \
    ${p_in_reg} \
    ${p_dbd_pid} \
    50 \
    single_dbds \
    ${p_out_net}npwa.adjmtr \
    ${p_src_code})

    job_id_avg_np=$(echo ${job_avg_np} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |        **** Weighted average bart ***            | #
# ==================================================== #
if [[ ${flag_avg_bart} == 1 || ${flag_all} == 1 ]]
then
    job_avg_bart=$(sbatch \
    -o ${p_wd}net_logs/${data}_avg_bart_%A.out \
    -e ${p_wd}net_logs/${data}_avg_bart_%A.err \
    -J ${data}_avg_bart \
    --dependency=afterany:${job_id_bart} \
    ${p_src_code}CODE/weighted_avg_similar_dbds_wrap.sh \
    ${p_out_net}${f_out_net_bart} \
    ${p_in_reg} \
    ${p_dbd_pid} \
    50 \
    single_dbds \
    ${p_out_net}bnwa.adjmtr \
    ${p_src_code})

    job_id_avg_bart=$(echo ${job_avg_bart} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |           **** Combine np & bart ***             | #
# ==================================================== #
if [[ ${flag_combine} == 1 || ${flag_all} == 1 ]]
then
    job_combine=$(sbatch \
    -o ${p_wd}net_logs/${data}_combine_np_bart_%A.out \
    -e ${p_wd}net_logs/${data}_combine_np_bart_%A.err \
    -J ${data}_combine_np_bart \
    --dependency=afterany:${job_id_avg_np}:${job_id_avg_bart} \
    ${p_src_code}CODE/quantile_combine_networks_wrap.sh \
    ${p_out_net}npwa.adjmtr \
    ${p_out_net}bnwa.adjmtr \
    ${p_out_net}npwa_bnwa.adjmtr)

    job_id_combine=$(echo ${job_combine} | awk '{split($0, a, " "); print a[4]}')
fi

# ==================================================== #
# |              **** Infer Motif ***                | #
# ==================================================== #
num_regulators=$(wc -l ${p_in_reg} | cut -d" " -f1)

if [[ ${flag_infer_motif} == 1 || ${flag_all} == 1 ]]
then
    job_infer_motif=$(sbatch \
    -o ${p_wd}net_logs/${data}_infer_motif_%A_%a.out \
    -e ${p_wd}net_logs/${data}_infer_motif_%A_%a.err \
    -J ${data}_infer_motif \
    --array=1-${num_regulators}%48 \
    --dependency=afterany:${job_id_combine} \
    ${p_src_code}CODE/run_infer_motifs.sh \
    ${p_out_tmp} \
    ${p_out_net}npwa_bnwa.adjmtr \
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
    -o ${p_wd}net_logs/${data}_score_motif_%A_%a.out \
    -e ${p_wd}net_logs/${data}_score_motif_%A_%a.err \
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
    -o ${p_wd}net_logs/${data}_build_motif_net_%A.out \
    -e ${p_wd}net_logs/${data}_build_motif_net_%A.err \
    -J ${data}_build_motif_net \
    --dependency=afterany:${job_id_score_motif} \
    ${p_src_code}CODE/build_motif_network_wrap.sh \
    ${p_out_tmp}motif_inference/motifs.txt \
    ${p_in_reg} \
    ${p_in_target} \
    ${p_out_tmp}motif_inference/motifs_score/ \
    robust \
    16 \
    ${p_out_net}mn.adjmtr)

    job_id_build_motif_net=$(echo ${job_build_motif_net} | awk '{split($0, a, " "); print a[4]}')
fi


# ==================================================== #
# |         **** Assemble final network ***          | #
# ==================================================== #
if [[ ${flag_combine_motif_net} == 1 || ${flag_all} == 1 ]]
then
    job_build_final_net=$(sbatch \
    -o ${p_wd}net_logs/${data}_build_final_net_%A.out \
    -e ${p_wd}net_logs/${data}_build_final_net_%A.err \
    -J ${data}_build_final_net \
    --dependency=afterany:${job_id_build_motif_net} \
    ${p_src_code}CODE/combine_networks_wrap.sh \
    resort \
    ${p_out_net}npwa_bnwa.adjmtr \
    ${p_out_net}mn.adjmtr \
    ${p_out_net} \
    npwa_bnwa_mn.adjmtr \
    ${p_src_code})

    job_id_build_final_net=$(echo ${job_build_final_net} | awk '{split($0, a, " "); print a[4]}')
fi

# ===================================================== #
# |        **** avg Assemble final network ***        | #
# ===================================================== #
if [[ ${flag_avg_final_net} == 1 || ${flag_all} == 1 ]]
then
    job_avg_final_net=$(sbatch \
    -o ${p_wd}net_logs/${data}_avg_final_net_%A.out \
    -e ${p_wd}net_logs/${data}_avg_final_net_%A.err \
    -J ${data}_avg_final_net \
    --dependency=afterany:${job_id_build_final_net} \
    ${p_src_code}CODE/weighted_avg_similar_dbds_wrap.sh \
    ${p_out_net}npwa_bnwa_mn.adjmtr \
    ${p_in_reg} \
    ${p_dbd_pid} \
    50 \
    single_dbds \
    ${p_out_net}${f_out_net_netprophet2} \
    ${p_src_code})
fi