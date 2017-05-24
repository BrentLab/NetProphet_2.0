###### multiple response LARS (single combined)
##TODO: Need to make # of folds a parameter exposed to used
library(lars)

lars.multi.optimize <- function(tdata,rdata,pert,prior,allowed)
{
	cat(as.character(Sys.time()),"\n")
	lars.paths.cv <- lars.multi.cv(tdata,rdata,pert,prior,allowed,10)
	
	lars.paths <- lars.multi(tdata,rdata,pert,prior,allowed)
	B <- lars.multi.path.step(lars.paths,lars.paths.cv[[2]])

	B.adj <- matrix(0,nrow=dim(rdata)[1],ncol=dim(tdata)[1])

	for (i in seq(dim(tdata)[1]))
	{
		#Scale B.adj according to prior
		B.adj[,i] <- B[[i]] * prior[,i]
	}
	rval <- list()
	rval[[1]] <- B.adj
	rval[[2]] <- lars.paths
	rval[[3]] <- lars.paths.cv
	cat(as.character(Sys.time()),"\n")
	rval
}

lars.multi.cv <- function(tdata,rdata,pert,prior,allowed,fold)
{	
	all.folds <- cv.folds(dim(tdata)[2], fold)
	lars.cv.paths <- list()
	cMax <- 0;
	for (f in 1:fold)
	{
		cat("Building path for cv",f,"of",fold,"\n")
		omit <- all.folds[[f]]
		lars.cv.paths[[f]] <- lars.multi(tdata[,-omit], rdata[,-omit], pert[,-omit],prior,allowed)
		if (lars.cv.paths[[f]]$R[1,1]>cMax)
		{
			cMax <- lars.cv.paths[[f]]$R[1,1]
		}
		
	}

	cSteps <- c()
	mse <- c()
	bnorm <- c()
	
	converged <- FALSE
	P <- c(cMax/10,0)
	FP <- c(lars.multi.cv.eval(P[1],fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed),lars.multi.cv.eval(P[2],fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed))

	while (!converged)
	{
		Imin <- which.min(FP)
		Imax <- which.max(FP)
		if (Imin==Imax)
			converged=TRUE

		cat("P:",P,"F(P):",FP,"Max/Min",Imax,Imin,"\n")
		Phat <- P[Imin]
		Pref <- 2*Phat - P[Imax]

		if (Pref>cMax)
			FPref <- max(FP) + (Pref-cMax)^2
		else if (Pref<0)
			FPref <- max(FP) + Pref^2
		else
			FPref <- lars.multi.cv.eval(Pref,fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed)
		
			if (FPref < FP[Imin])		
			{
				#Attempt expansion
				Pexp <- 2*Pref - Phat
				FPexp <- lars.multi.cv.eval(Pexp,fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed)
				if (FPexp < FPref) #Expand
				{
					cat("Expand\n")
					P[Imax] <- Pexp
					FP[Imax] <- FPexp
				}
				else #Reflect
				{
					cat("Reflect\n")
					P[Imax] <- Pref
					FP[Imax] <- FPref
				}
			}
			else
			{
				#Contract
				if (FP[Imax] < FPref)
				{
					cat("Contract (1)\n")
					P[Imax] <- (P[Imax] + Phat)/2
					FP[Imax] <- lars.multi.cv.eval(P[Imax],fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed)
				}
				else
				{
					cat("Contract (2)\n")
					P[Imax] <- (Pref + Phat)/2
					FP[Imax] <- lars.multi.cv.eval(P[Imax],fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed)
				}
			}
	}
	rval <- list()
	rval[[1]] <- lars.cv.paths
	rval[[2]] <- P[which.min(FP)]
	rval[[3]] <- min(FP)
	rval[[4]] <- all.folds
	rval
}

lars.multi.cv.eval <- function(cVal,fold,lars.cv.paths,all.folds,tdata,rdata,pert,prior,allowed,targets=NULL)
{
		if (is.null(targets))
			targets <- seq(dim(tdata)[1])
		mse <- rep(0,times=fold)
		for (f in 1:fold)
		{
			B <- lars.multi.path.step(lars.cv.paths[[f]],cVal)
			for (i in targets)
			{
				omit <- all.folds[[f]]
				omit <- setdiff(omit,which(pert[i,]!=0))
				
				testx <- rdata[,omit]
				testy <- tdata[i,omit]
				testx[which(allowed[,i]==0),] <- 0
				testx <- testx * prior[,i]

				testp <- scale(t(testx) , lars.cv.paths[[f]][[i]]$meanx, FALSE) %*% matrix(B[[i]]) + lars.cv.paths[[f]][[i]]$mu
				testr <- apply((testp - testy)^2,2,mean)
				mse[f] <- mse[f] + testr
			}
			mse[f] <- mse[f]/length(targets)
		}
		totalmse <- sum(mse) / fold
		totalmse
}


lars.multi <- function(tdata,rdata,pert,prior,allowed)
{
	lars.paths <- list()
	targets <- seq(dim(tdata)[1])
	

	cat("Building path for target: ")
	for (i in targets)
	{
		cat(i,"")
		mindices <- which(pert[i,]==0)

		x <- rdata[,mindices] * prior[,i] 
		x[which(allowed[,i]==0),]<-0; 

		##### YK patch 05-23-2017
		# results.lars <- lars(t(x),tdata[i,mindices],trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE,use.Gram=FALSE);
		y <- tdata[i,mindices];
		results.lars.flag <- TRUE;
		sample.indx <- length(y); 
		while (results.lars.flag) {
			results.lars <- try(lars(t(x[,1:sample.indx]),y[1:sample.indx],trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE,use.Gram=FALSE));
			results.lars.flag <- class(results.lars) =='try-error';
			if (results.lars.flag) { cat(' *no soln* '); }
			sample.indx <- sample.indx-1; 
		}
		#####
		lars.paths[[i]] <- list(beta=results.lars$beta, mu=results.lars$mu, meanx=results.lars$meanx, C.Max=c(results.lars$C.Max,0))
	}
	cat("\n")

	cat("Building combined C.Max ranking\n")
	all.C.lars <- c()
	for (i in targets)
	{
		all.C.lars <- rbind(all.C.lars, cbind(lars.paths[[i]]$C.Max,rep(i,times=length(lars.paths[[i]]$C.Max)),1:length(lars.paths[[i]]$C.Max)))
	}
	R <- all.C.lars[sort.list(all.C.lars[,1],decreasing=TRUE),]
	lars.paths$R <- R
	lars.paths$targets <- targets
	lars.paths
}

lars.multi.path.step <- function(lars.paths,stepC)
{
	k <- 1
	if (length(which(lars.paths$R[,1]>stepC))>0)
	{
		k<-max(which(lars.paths$R[,1]>stepC))
	}

	B <- list()
	A <- c()
	Ai <- unique(lars.paths$R[1:k,2])
	for (i in Ai)
	{	
		A <- rbind(A,c(i,lars.paths$R[max(which(lars.paths$R[1:k,2]==i)),3]))
	}

	for (i in lars.paths$targets)
	{
		indx <- which(A[,1]==i)
		if (length(indx)==0)
		{
			B[[i]] <- rep(0,times=dim(lars.paths[[i]]$beta)[2])
		}
		else
		{
			plen <- dim(lars.paths[[i]]$beta)[1]
			if (A[indx,2]==plen)
			{
				B[[i]] <- lars.paths[[i]]$beta[plen,]
			}
			else
			{
				if ( lars.paths$R[k,1] < lars.paths[[i]]$C.Max[A[indx,2]+1] )
				{
					if (A[indx,2] < plen)
						A[indx,2] <- A[indx,2] + 1
				}
				gamma <- (lars.paths[[i]]$C.Max[A[indx,2]] - lars.paths$R[k,1]) / (lars.paths[[i]]$C.Max[A[indx,2]] - lars.paths[[i]]$C.Max[A[indx,2]+1])
				u <- (lars.paths[[i]]$beta[A[indx,2]+1,] - lars.paths[[i]]$beta[A[indx,2],])
				B[[i]] <- lars.paths[[i]]$beta[A[indx,2],] + u * gamma
			}
		}
	}
	B
}

###

lars.multi.cv.singlefold <- function(tdata,rdata,pert,prior,allowed,fold,seed,f)
{
	set.seed(seed)
	all.folds <- cv.folds(dim(tdata)[2], fold)
	cat("Building path for cv",f,"of",fold,"\n")
	omit <- all.folds[[f]]
	lars.cv.paths <- lars.multi(tdata[,-omit], rdata[,-omit], pert[,-omit],prior,allowed)
	lars.cv.paths
}

lars.single.cv <- function(tdata,rdata,pert,prior,allowed,fold,seed)
{
	targets <- seq(dim(tdata)[1])
	lars.paths.cv <- list()
	lars.paths <- list()
	B <- matrix(0,nrow=dim(rdata)[1],ncol=dim(tdata))

	minFrac <- c()
  cat("Building path for target: ")
  for (i in targets)
  {
    cat(i,"")
    mindices <- which(pert[i,]==0)
    x <- rdata[,mindices] * prior[,i]
    x[which(allowed[,i]==0),]<-0;
    set.seed(seed)
    lars.paths.cv[[i]] <-  cv.lars(t(x),tdata[i,mindices],K=fold,plot.it=FALSE,trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE,se=FALSE);
		minFrac[i] <- lars.paths.cv[[i]]$index[which.min(lars.paths.cv[[i]]$cv)]
		lars.paths[[i]] <- lars(t(x),tdata[i,mindices],trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE);
		B[,i] <- predict(lars.paths[[i]],s=minFrac[i],type="coefficients",mode="fraction")$coefficients
  }	
	cat("\n")
	B
}

## Parallel functionality
lars.multi.cv.eval.parallel <- function(cVal)
{
	mse <- 0;
	B <- lars.multi.path.step(cv.obj,cVal)
	for (i in targets)
	{
		omit <- all.folds[[f]]
		omit <- setdiff(omit,which(pert[i,]!=0))
		
		testx <- rdata[,omit]
		testy <- tdata[i,omit]
		testx[which(allowed[,i]==0),] <- 0
		testx <- testx * prior[,i]

		testp <- scale(t(testx) , cv.obj[[i]]$meanx, FALSE) %*% matrix(B[[i]]) + cv.obj[[i]]$mu
		testr <- apply((testp - testy)^2,2,mean)
		mse <- mse + testr
	}
	mse <- mse/length(targets)
	mse
}

lars.multi.cv.cmax.parallel <- function()
{
	cat(cv.obj$R[1,1],"\n")
	cv.obj$R[1,1]
}

lars.multi.cv.parallel <- function()
{	
	cMax <- max(as.numeric(mpi.remote.exec(lars.multi.cv.cmax.parallel())))
	cat("cMax: ", cMax, "\n")
	converged <- FALSE
	P <- c(cMax/10,0)
	
	FP <- c(mean(as.numeric(mpi.remote.exec(cmd=lars.multi.cv.eval.parallel,P[1]))),mean(as.numeric(mpi.remote.exec(cmd=lars.multi.cv.eval.parallel,P[2]))))
	
	while (!converged)
	{
		Imin <- which.min(FP)
		Imax <- which.max(FP)
		if (Imin==Imax)
			converged=TRUE
	
		cat("P:",P,"F(P):",FP,"Max/Min",Imax,Imin,"\n")
		Phat <- P[Imin]
		Pref <- 2*Phat - P[Imax]
	
		if (Pref>cMax)
			FPref <- max(FP) + (Pref-cMax)^2
		else if (Pref<0)
			FPref <- max(FP) + Pref^2
		else
			FPref <- mean(as.numeric(mpi.remote.exec(cmd=lars.multi.cv.eval.parallel,Pref)))
		
			if (FPref < FP[Imin])		
			{
				#Attempt expansion
				Pexp <- 2*Pref - Phat
				FPexp <- mean(as.numeric(mpi.remote.exec(cmd=lars.multi.cv.eval.parallel,Pexp)))
				if (FPexp < FPref) #Expand
				{
					cat("Expand\n")
					P[Imax] <- Pexp
					FP[Imax] <- FPexp
				}
				else #Reflect
				{
					cat("Reflect\n")
					P[Imax] <- Pref
					FP[Imax] <- FPref
				}
			}
			else
			{
				#Contract
				if (FP[Imax] < FPref)
				{
					cat("Contract (1)\n")
					P[Imax] <- (P[Imax] + Phat)/2
					FP[Imax] <- mean(as.numeric(mpi.remote.exec(cmd=lars.multi.cv.eval.parallel,P[Imax])))
				}
				else
				{
					cat("Contract (2)\n")
					P[Imax] <- (Pref + Phat)/2
					FP[Imax] <- mean(as.numeric(mpi.remote.exec(cmd=lars.multi.cv.eval.parallel,P[Imax])))
				}
			}
	}
	
	cat(P[which.min(FP)],min(FP),"\n")

	P[which.min(FP)]
}

lars.multi.optimize.parallel <- function()
{
	cat(as.character(Sys.time()),"\n")
	c.cvmin <- lars.multi.cv.parallel()
	
	lars.paths <- lars.multi(tdata,rdata,pert,prior,allowed)
	B <- lars.multi.path.step(lars.paths,c.cvmin)
	
	B.adj <- matrix(0,nrow=dim(rdata)[1],ncol=dim(tdata)[1])
	
	for (i in seq(dim(tdata)[1]))
	{
		#Scale B.adj according to prior
		B.adj[,i] <- B[[i]] * prior[,i] 
	}
	rval <- list()
	rval[[1]] <- B.adj
	rval[[2]] <- lars.paths
	cat(as.character(Sys.time()),"\n")
	rval
}

lars.local <- function(tdata,rdata,pert,prior,allowed,skip_reg,skip_gen)
{
  cat(as.character(Sys.time()),"\n")
  B.adj <- matrix(0,nrow=dim(rdata)[1],ncol=dim(tdata)[1])
	
  cat("Working on Gene:")
	for(i in 1:dim(tdata)[1]) {
		if (skip_gen[i] == 0) {
		 	cat(i,",")
			mindices <- which(pert[i,]==0)
			x <- rdata[,mindices] * prior[,i]
			x[which(allowed[,i]==0),]<-0;
			nindices <- which(skip_reg==0)
			x <- x[nindices,]
			lars.paths.cv <-  cv.lars(t(x),tdata[i,mindices],K=3,trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE,se=FALSE, plot.it=FALSE,use.Gram=FALSE);
			minFrac <- lars.paths.cv$index[which.min(lars.paths.cv$cv)]
			lars.paths <- lars(t(x),tdata[i,mindices],trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE,use.Gram=FALSE);
			# lars.paths <- lars(t(x),tdata[i,mindices],trace=FALSE,max.steps=600,type="lasso",normalize=FALSE,intercept=TRUE,use.Gram=TRUE);
			tempVec <- predict(lars.paths,s=minFrac,type="coefficients",mode="fraction")$coefficients
			# tempVec <- lars.paths$beta
			# print(tempVec)
			nindices <- which(skip_reg==1)
      tempCol <- c(tempVec, rep(0,length(nindices)))
      tempIndices <- c(seq_along(tempVec), nindices+.5)
      B.adj[,i] <- tempCol[order(tempIndices)]
      #Scale B.adj according to prior		
			B.adj[,i] <- B.adj[,i] * prior[,i]
		}
	}
	cat("\n")

	rval <- list()
  rval[[1]] <- B.adj
  rval
}

