source("global.lars.regulators.r")
## TODO: Both seed and # of cv folds (in global.lars.regulators.r) are parameters that should be exposed to the user
seed <- 747
fold <- 10
tdata <- as.matrix(read.table(targetExpressionFile))
rdata <- as.matrix(read.table(regulatorExpressionFile))
allowed <- as.matrix(read.table(allowedMatrixFile))
pert <- as.matrix(read.table(perturbationMatrixFile))
targets <- seq(dim(tdata)[1])

if(microarrayFlag == 0) {
	##RNA-Seq Data
	tdata <- log(tdata+1)/log(2)
	rdata <- log(rdata+1)/log(2)
}

## Center data
tdata <- tdata - apply(tdata,1,mean)
rdata <- rdata - apply(rdata,1,mean)

## Scale data
t.sd <- apply(tdata,1,sd)
t.sdfloor <- mean(t.sd) + sd(t.sd)
t.norm <- apply(rbind(rep(t.sdfloor,times=length(t.sd)),t.sd),2,max) / sqrt(dim(tdata)[2])
tdata <- tdata / ( t.norm * sqrt(dim(tdata)[2]-1) )
#
r.sd <- apply(rdata,1,sd)
r.sdfloor <- mean(r.sd) + sd(r.sd)
r.norm <- apply(rbind(rep(r.sdfloor,times=length(r.sd)),r.sd),2,max) / sqrt(dim(rdata)[2])
rdata <- rdata / ( r.norm * sqrt(dim(rdata)[2]-1) )

## Compute unweighted solution
prior <- matrix(1,ncol=dim(tdata)[1] ,nrow=dim(rdata)[1] )

set.seed(seed)

all.folds <- cv.folds(dim(tdata)[2], fold)

