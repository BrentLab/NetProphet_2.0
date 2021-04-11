#!/bin/bash

p_in_expr_target=NetProphet_2.0/toy_example/zev_expr_500_100_indexed
p_in_expr_reg=NetProphet_2.0/toy_example/zev_expr_reg_50_100_indexed

# lasso
flag_global_shrinkage="ON"
flag_local_shrinkage="OFF"

# smoothing
p_in_dbd_pids=NetProphet_2.0/toy_example/DBD_PIDS

# pwn
p_in_promoter=NetProphet_2.0/toy_example/promoter.scer.fasta

# logistics
p_src_code=/home/NetProphet_2.0/

# singularity
flag_singularity="ON"
p_singularity_img=/home/vagrant/s_np3
p_singularity_bindpath=/home/vagrant/

# output
p_out_dir=/home/vagrant/res_toy_example/

${p_src_code}np2 \
    --p_in_expr_target ${p_in_expr_target} \
    --p_in_expr_reg ${p_in_expr_reg} \
    --flag_global_shrinkage ${flag_global_shrinkage} \
    --flag_local_shrinkage ${flag_local_shrinkage} \
    --p_in_dbd_pids ${p_in_dbd_pids} \
    --p_in_promoter ${p_in_promoter} \
    --p_src_code ${p_src_code} \
    --flag_singularity ${flag_singularity} \
    --p_singularity_img ${p_singularity_img} \
    --p_singularity_bindpath ${p_singularity_bindpath} \
    --p_out_dir ${p_out_dir}
