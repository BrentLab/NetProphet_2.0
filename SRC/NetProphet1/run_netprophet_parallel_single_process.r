source("run_netprophet_parallel.r")
f <- mpi.comm.rank()
cat(f,"\n")

#load(paste("cv.",f,sep=""))
cv.obj <- lars.multi.cv.singlefold(tdata,rdata,pert,prior,allowed,fold,seed,f)
save(cv.obj,file=file.path(outputDirectory, paste("cv.",f,sep="")))

reportid <- function()
{
	f
}

