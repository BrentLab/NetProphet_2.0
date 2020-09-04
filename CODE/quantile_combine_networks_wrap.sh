#!/bin/bash

module load R/3.4.3

p_npwa=${1}
p_bnwa=${2}
p_npwa_bnwa=${3}

Rscript /scratch/mblab/dabid/netprophet/code_netprophet2.1/CODE/quantile_combine_networks.r \
${p_npwa} \
${p_bnwa} \
${p_npwa_bnwa}