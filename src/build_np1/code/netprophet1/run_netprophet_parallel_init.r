args <- commandArgs(trailingOnly = TRUE)
targetExpressionFile <- toString(args[1])
regulatorExpressionFile <- toString(args[2])
allowedMatrixFile <- toString(args[3])
perturbationMatrixFile <- toString(args[4])
differentialExpressionMatrixFile <- toString(args[5])
microarrayFlag <- as.integer(args[6])
nonGlobalShrinkageFlag <- as.integer(args[7])
lassoAdjMtrFileName <- toString(args[8])
combinedAdjMtrFileName <- toString(args[9])
outputDirectory <- toString(args[10])
combinedAdjLstFileName <- toString(args[11])
regulatorGeneNamesFileName <- toString(args[12])
targetGeneNamesFileName <- toString(args[13])

source("/scratch/mblab/dabid/netprophet/code_netprophet3.0/code/netprophet1/run_netprophet_parallel.r")
#install.packages("/scratch/mblab/dabid/netprophet/code_NetProphet_2.0/SRC/R_pkgs/Rmpi_0.5-9.tar.gz", type="source", repos=NULL)

library(Rmpi)
mpi.bcast.Robj2slave(targetExpressionFile)
mpi.bcast.Robj2slave(regulatorExpressionFile)
mpi.bcast.Robj2slave(allowedMatrixFile)
mpi.bcast.Robj2slave(perturbationMatrixFile)
mpi.bcast.Robj2slave(microarrayFlag)
mpi.bcast.Robj2slave(outputDirectory)

mpi.remote.exec(source("/scratch/mblab/dabid/netprophet/code_netprophet3.0/code/netprophet1/run_netprophet_parallel_single_process.r"))
mpi.remote.exec(reportid())
uniform.solution <- lars.multi.optimize.parallel()

lasso_component <- uniform.solution[[1]]
write.table(lasso_component,file.path(outputDirectory,lassoAdjMtrFileName),row.names=FALSE,col.names=FALSE,quote=FALSE)
#save(solution,file="solution")

de_component <- as.matrix(read.table(differentialExpressionMatrixFile))

## Perform model averaging to get final NetProphet Predictions
source("combine_models.r")

# if(length(args) == 13 & file.exists(regulatorGeneNamesFileName) & file.exists(targetGeneNamesFileName)){
#   source("make_adjacency_list.r")
# }
