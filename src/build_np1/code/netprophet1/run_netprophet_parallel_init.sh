#!/bin/bash
#SBATCH -n 11
#SBATCH --cpus-per-task=2
#SBATCH -D ./SRC/NetProphet1/
#SBATCH -o ../../LOG/map_netprophet1_network_%A.out
#SBATCH -e ../../LOG/map_netprophet1_network_%A.err

module load R/3.2.1
module load openmpi/1.8.3

targetExpressionFile=${1}
regulatorExpressionFile=${2}
allowedMatrixFile=${3}
perturbationMatrixFile=${4}
differentialExpressionMatrixFile=${5}
microarrayFlag=${6}
nonGlobalShrinkageFlag=${7}
lassoAdjMtrFileName=${8}
combinedModelAdjMtrFileName=${9}
outputDirectory=${10}
combinedAdjLstFileName=${11}
regulatorGeneNamesFileName=${12}
targetGeneNamesFileName=${13}


module load R/3.2.1
module load openmpi/1.8.8

echo "calling mpirun now, SLURM_NTASKS=${SLURM_NTASKS}"

mpirun -np ${SLURM_NTASKS} R --no-save -q --args ${targetExpressionFile} ${regulatorExpressionFile} ${allowedMatrixFile} ${perturbationMatrixFile} ${differentialExpressionMatrixFile} ${microarrayFlag} ${nonGlobalShrinkageFlag} ${lassoAdjMtrFileName} ${combinedModelAdjMtrFileName} ${outputDirectory} ${combinedAdjLstFileName} ${regulatorGeneNamesFileName} ${targetGeneNamesFileName} < run_netprophet_parallel_init.r

