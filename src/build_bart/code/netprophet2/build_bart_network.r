# if (!require("R.oo")) try(install.packages("R.oo"));
# if (!require("Rmpi")) try(install.packages("Rmpi"));
# if (!require("BayesTree")) try(install.packages("BayesTree"));
# if (!require("stats")) try(install.packages("stats"));
# if (!require("matrixStats")) try(install.packages("matrixStats"));
# if (!require("Matrix")) try(install.packages("Matrix"));
# if (!require("abind")) try(install.packages("abind"));
# if (!require("restorepoint")) try(install.packages("restorepoint")); ## deprecated
getBartNetwork <- function(tgtLevel, tfLevel, regMat, unperturbedTfLevel, nBin = 3, ...) {
	# use BART to generate a network structure prediction
	# tgtLevel: training target expression levels, in foldchange with respect to wildtype
	# tfLevel: training regulator expression levels, in foldchange with respect to wildtype, the log of which will be fed to BART
	# regMat: adjacency matrix of a prior allowed interactions: ij entry represens wether regulator i is allowed to regulate target j. If missing, everything except self-regulation is allowed
	# ...: additional parameters to pass to bartExpr
	# will construct query TF level for network reconstruction
	#	each query profile contains one and only one perturbed regulator levels
	#	unperturbedTfLevel: specifies the unperturbed regulator levels; defaults are the medians of the regulators
	#	perturbed levels of a regulator taken from its quantiles in the training data, with equally distanced probs from 0 to 1
	#	nBin: specifies the number of quantiles to use for each regulator, default is 3. (giving min, median, and max)
	#	yMean (predicted transformed tgtLevel) and yVar(posterior variance of transformed tgtLevel) will be reshaped so that:
	#		first dimension corresponds to different quantiles of regulators they are responding to
	#		second dimension corresponds to the identity of the perturbed regulator they are responding to
	#		third dimenison is the idenity of the responding target gene itself
	# returns fields including those of bartExpr, with an additional field regScore which is the signed scores of regulation likelihood
    set.seed(747)
    # load libraries
    library("R.oo")
    library("Rmpi")
    library("BayesTree")
    library("stats")
    library("matrixStats")
    library("Matrix")
    library("abind")
	# require("matrixStats");
	nTf <- ncol(tfLevel);
	#
	# constructing query regulator levels
	result <- list();
	result$perturbedTfLevel <- t(colQuantiles(tfLevel, probs = seq(0, 1, length.out = nBin), drop = FALSE)); # quantiles of regulator levels
	testTfLevel <- rbind(colMedians(tfLevel))[rep(1, nBin * nTf), , drop = FALSE]; # shaping testTfLevel
	if (!missing(unperturbedTfLevel)) testTfLevel[] <- c(unperturbedTfLevel); # if unperturbedTfLevel is specified, use it; otherwise, medians of TF levels are used
	ix <- cbind(sequence(nrow(testTfLevel)), rep(sequence(nTf), each = nBin)); # indices of perturbed entries in the matrix
	colnames(result$perturbedTfLevel) <- colnames(tfLevel);
	testTfLevel[ix] <- c(result$perturbedTfLevel); # filling perturbed TF levels to testTfLevel
	#
	# calling BART
	barted <- bartExpr(tgtLevel = tgtLevel, tfLevel = tfLevel, regMat = regMat, testTfLevel = testTfLevel,...);
	result <- c(result, barted);
	#
	# reshaping results
	dim(result$yMean) <- c(dim(result$perturbedTfLevel), ncol(tgtLevel)); # reshaping yMean so that the first dimension is regulator quantiles, second dimension is regulator identities, third dimension is targe identities
	dimnames(result$yMean) <- c(dimnames(result$perturbedTfLevel), list(colnames(tgtLevel)));
	dim(result$yVar) <- dim(result$yMean); # reshaping yVar
	dimnames(result$yVar) <- dimnames(result$yMean); 
	try(dim(result$yMeanVar) <- dim(result$yMean)); # reshaping yMeanVar
	try(dimnames(result$yMeanVar) <- dimnames(result$yMean));
	#
	result$regScore <- result$yMean[nBin, , ] - result$yMean[1, , ] # calculating regulation score
	result;
}

bartExpr <- function(tgtLevel, tfLevel, testTfLevel, regMat, verbose = TRUE, noiseModel = c("lognormal", "normal", "linear", "sqrt"), ...) {
	# use BART to predict expressions
	# expression: training target expression levels, in foldchange with respect to wildtype
	# tfLevel: training regulator expression levels, in foldchange with respect to wildtype, the log of which will be fed to BART
	# testTfLevel: query regulator expression levels
	# regMat: adjacency matrix of a prior allowed interactions: ij entry represens wether regulator i is allowed to regulate target j. If missing, everything except self-regulation is allowed
	# noiseModel: "lognormal", "normal", or an object of class AhnNoiseModel
	# returned fields:
	#	yMean: predicted transformed target level
	#	yVar: posterior variance of target level in prediction
	#	predicted: predicted target level in foldchange with respect to wild type
	result <- list();
	if (is.character(noiseModel)) { # determine transform and backTransform from noiseModel
		noiseModel <- match.arg(noiseModel);
		if (noiseModel == "normal" || noiseModel == "linear") {
			transform <- function(x) x;
			backTransform <- function(x) x;
		} else {
			if (noiseModel == "lognormal") {
				transform <- function(x) log(x);
				backTransform <- function(x) exp(x);
			} else {
				if (noiseModel == "sqrt") {
					transform <- function(x) log(sqrt(x));
					backTransform <- function(x) pmax(x, 0)^2;
				}
			}
		}
	} else {
		if ("HbNoiseTransform" %in% class(noiseModel)) {
			transform <- noiseModel$toNormal;
			backTransform <- noiseModel$fromNormal;
		} else {
			if (is.list(noiseModel)) {
				transform <- noiseModel$transform;
				backTransform <- noiseModel$backTransform;
			}
		}
	}
	# calling BART
	barted <- bartMultiresponse(x.train = log(tfLevel), y.train = transform(tgtLevel), x.test = log(testTfLevel), allowed = regMat, verbose = verbose, simplify = TRUE,...); 
	result <- c(result, barted);
	result$predicted <- backTransform(barted$yMean); # transform the prediction back to the space of fold changes
	result;
}

bartMultiresponse <- function(x.train, y.train, x.test = NULL, allowed, simplify = TRUE, verbose = TRUE, mpiComm, blockSize, saveTo, ...) {
	# wrapper for using bart on multiple responses, able to invoke mpi
	# allowed: matrix whose ij entry tells if explanatory variable i should be used for response j; by default, if variables names are given in the colnames of x.train and y.train, every explanatory relation is allowed except for a variable explaining itself
	# mpiComm: missing or integer; if specified, will invoke mpi on the given comm number mpiComm (see package Rmpi for details)
	# blockSize: when using mpiComm, each step in the interation a number blockSize of response variables will be dispatched to mpi processes for bart calculation. Default is approximately square root of number of reponse variables
	# saveTo: a file to store current result during calculation; since calculation may take along, partial results are stored in case calculation gets terminated accidentally
	# simplify: same meaning as in function bartRobust, recommend simplify = TRUE other wise spatial complexity with be astronomical
	# verbose = 1: print progress only; verbose = 2: print progress and bart output
	nResponse <- ncol(y.train);
	nExplanatory <- ncol(x.train); 
	responseName <- colnames(y.train);
	explanatoryName <- colnames(x.train);
	result <- list(); # place to store results.
	bartedList <- list();
	if (missing(allowed)) {
		# in case allowed is not specified, all exlanatory relations are by default allowed except for self-explaining
		allowed <- matrix(TRUE, nExplanatory, nResponse, dimnames = list(explanatoryName, responseName));
		selfReg <- intersect(responseName, explanatoryName);
		allowed[cbind(selfReg, selfReg)] <- FALSE;
	}
	#
	# function to apply to each response variable for BART analysis
	applicand <- function(i, ...) {
# 		require("Matrix");
        library("Matrix")
		theX <- as.matrix(x.train[, allowed[, i], drop = FALSE]);
		theY <- y.train[, i];
		theTestX <- as.matrix(x.test[, allowed[, i], drop = FALSE]);
		bartRobust(
			x.train = theX,
			y.train = theY,
			x.test = theTestX,
			...
		);
	}
	#
	if (missing(mpiComm) || is.null(mpiComm)) { # calculate without invoking mpi
# 		require("foreach");
# 		require("doParallel");
        library("foreach")
        library("doParallel")
        
		## define pool
		pool <- makeCluster(min(11, detectCores()[1]-1), outfile="")
		registerDoParallel(pool)
		print(pool)

		## pre-convert list names
		gene2indx <- list()
		for (i in seq(nResponse)) {
			gene2indx[[responseName[i]]] <- i;
		}
		## parallel for loop
		pb <- txtProgressBar(0, nResponse, style=3)
		bartedList <- foreach(i=gene2indx, 
							.inorder=T, 
							.final=function(x) setNames(x, names(gene2indx)),
							.export=c("bartRobust")
							) %dopar% {
			# if (verbose) print(paste(i, "/", nResponse, responseName[i]));
			if (verbose) setTxtProgressBar(pb, i);
			try(applicand(i, simplify = simplify, verbose = pmax(verbose - 1, 0), ...));
			}
		## exit pool
		stopCluster(pool)

		## old serial processing
		# for (i in sequence(nResponse)) {
		# 	if (verbose) print(paste(i, "/", nResponse, responseName[i]));
		# 	bartedList[[responseName[i]]] <- try(applicand(i, simplify = simplify, verbose = pmax(verbose - 1, 0), ...));
		# 	if (!missing(saveTo) && !is.null(saveTo)) save(bartedList, file = saveTo);
		# }
	} else { # calculate with mpi
# 		require("Rmpi");
        library("Rmpi")
		# try(Sys.setenv(OMPI_MCA_btl_tcp_if_include="eth0"));        
		if (mpi.comm.size(comm = mpiComm) == 0) { # initialize mpi comm if have not done so yet        
			if (missing(blockSize)) blockSize <- floor(sqrt(nResponse));
			mpi.spawn.Rslaves(nslaves = blockSize, comm = mpiComm);
		}
		nBlock <- ceiling(nResponse / blockSize);
		blockLabel <- head(rep(seq(nBlock), each = blockSize), nResponse);
		blockIxList <- split(seq(nResponse), blockLabel); # a list of which each element is index set of one block
		mpi.bcast.Robj2slave(x.train, comm = mpiComm);
		mpi.bcast.Robj2slave(y.train, comm = mpiComm);
		mpi.bcast.Robj2slave(x.test, comm = mpiComm);
		mpi.bcast.Robj2slave(allowed, comm = mpiComm);
		mpi.bcast.Robj2slave(bartRobust, comm = mpiComm);
		try(mpi.bcast.Rfun2slave(comm = mpiComm));
		environment(applicand) <- globalenv(); # very important, so function applicand broadcasted to other mpi processes will read variables from their global environments instead of bringing this environment with it
		for (i in seq(nBlock)) {
			if (verbose) print(paste(i, "/", nBlock));
			bartedList <- c(bartedList, mpi.apply(blockIxList[[i]], applicand, simplify = simplify, verbose = pmax(verbose - 1, 0), comm = mpiComm, ...));
			if (!missing(saveTo) && !is.null(saveTo)) save(bartedList, file = saveTo);
		}
	}
	names(bartedList) <- responseName;
	bartedList <- lapply(bartedList, as.list); # in case one bart call returns an error message, the error message will be covnerted to a list
	if (!missing(saveTo) && !is.null(saveTo)) save(bartedList, file = saveTo);
	if (simplify) { # simplify results
		result$yMean <- cbindCountingNull(lapply(bartedList, "[[", "yMean"));
		result$yMeanVar <- cbindCountingNull(lapply(bartedList, "[[", "yMeanVar"));
		result$yVar <- cbindCountingNull(lapply(bartedList, "[[", "yVar"));
		colnames(result$yMean) <- responseName;
		rownames(result$yMean) <- rownames(if (missing(x.test) || is.null(x.test)) x.train else x.test);
		dimnames(result$yMeanVar) <- dimnames(result$yMean);
		dimnames(result$yVar) <- dimnames(result$yMean);
		result$sigest <- c(cbindCountingNull(lapply(bartedList, "[[", "sigest")));
		names(result$sigest) <- responseName;
		result;
	} else {
		bartedList;
	}
}

bartRobust <- function(x.train, y.train, x.test, simplify = FALSE, keepBartObj = TRUE, keepevery = 1, verbose = TRUE, sigest = NA, ...) {
	# wrap bart to a more robust function, circumventing some bugs there, and adding some useful stats about the query prediction in fields in the result:
	# 	$barted: the orignal bart object returned
	# 	$yMean: mean of posterior prediction, rows are samples, ncol = 1;
	#	$y: posterior predictions, rows are samples, ncol = 1, third dimension is MCMC steps
	#	$yVar: posterior variance, rows are samples, ncol = 1
	# for internal usage, called by other functions.
	# x.test: query features, default is using training data
	# keepevery: thinning number, see bart in package BayesTree
	# sigest: may be left unspecified; see bart in package BayesTree
	# keepBartObj: if FALSE, result$barted will be removed (overridden by simplify = TRUE)
	# simplify: if TRUE, return only fields yMean, yVar, yMeanVar and sigest
# 	require("BayesTree");
# 	require("matrixStats");
    library("BayesTree")
    library("matrixStats")
    
	result <- list(); # place to store results.
	#
	#if (missing(x.test)) x.test <- x.train[c(), , drop = FALSE]; # use empty test data if not specified otherwise
	#
	nonNaIx <- c(!is.na(y.train)); 		#
	x <- x.train[nonNaIx, , drop = FALSE];		#
	y <- y.train[nonNaIx];				# remove NA entries
	if (ncol(x) > nrow(x) - 1 && is.na(sigest)) sigest <- sd(y, na.rm = TRUE); # in case ncol(x) > nrow(x) - 1, sigest get manually specified here, otherwise bart would run into a bug and crash
	if (length(unique(y)) > 1) { # only call bart if there's multiple values in training reponses y, otherwise bart would run into a bug and crash
		defaultTest <- if (missing(x.test) || is.null(x.test)) matrix(0, 0, 0) else x.test;
		result$barted <- bart(
			x.train = x,
			y.train = y,	
			x.test = defaultTest,
			keepevery = keepevery,
			verbose = (verbose != 0),
			sigest = sigest, 
			...
		);
		# summarizing results
		if (missing(x.test) || is.null(x.test)) { # summarize training result if x.test is not specified
			result$yMean <- cbind(y * NA);
			result$yMean[nonNaIx, ] <- result$barted$yhat.train.mean;
			rownames(result$yMean) <- rownames(x.train);
			result$y <- array(NA, c(dim(result$yMean), nrow(result$barted$yhat.train)));
			dimnames(result$y)[1:2] <- dimnames(result$yMean);
			result$y[nonNaIx, , ] <- c(t(result$barted$yhat.train));
			result$yMeanVar <- result$yMean * NA;
			result$yMeanVar[nonNaIx, ] <- colVars(result$barted$yhat.train);
		} else { # summarize query result if x.test is given
			result$yMean <- cbind(result$barted$yhat.test.mean); 
			rownames(result$yMean) <- rownames(x.test);
			result$y <- array(NA, c(dim(result$yMean), nrow(result$barted$yhat.test)));
			dimnames(result$y)[1:2] <- dimnames(result$yMean);
			result$y[] <- c(t(result$barted$yhat.test));
			result$yMeanVar <- result$yMean * NA;
			result$yMeanVar[] <- colVars(result$barted$yhat.test);
		}
		result$yVar <- mean(result$barted$sigma^2) + result$yMeanVar;
		result$sigest <- result$barted$sigest;
	} else { # trivial case, no need to run bart
		result$barted <- list();
		if (missing(x.test)) {
			result$yMean <- matrix(y[1], nrow(x.train), 1, dimnames = list(rownames(x.train), colnames(y.train)));
		} else {
			result$yMean <- matrix(y[1], nrow(x.test), 1, dimnames = list(rownames(x.test), colnames(y.train)));
		}
		result$y <- result$yMean;
		result$yMeanVar <- result$yMean * 0;
		result$yVar <- result$yMean * 0;
		result$sigest <- 0;
	}
	#
	if (!keepBartObj) result$barted <- NULL;
	if (simplify) result[c("yMean", "yVar", "yMeanVar", "sigest")] else result;	
}

robustBart <- bartRobust;

cbindCountingNull <- function(colList, default = NA) { # column-binding elements of colList, with null elements replaced by NA columns
	nullIx <- mapply(is.null, colList);
	nonNullCombined <- do.call("cbind", colList[!nullIx]);
	combined <- matrix(default, nrow(nonNullCombined), length(colList));
	rownames(combined) <- rownames(nonNullCombined);
	colnames(combined) <- names(colList);
	combined[, !nullIx] <- nonNullCombined;
	combined;
}

# save.image();

# # formating input arguments
# argStringList <- commandArgs()[7:length(commandArgs())]
# argTable <- read.table(text = unlist(argStringList), sep = "=", header = 0, stringsAsFactors = FALSE, row.names = 1);
# write.table(t(argTable), file = textConnection("transposedArgString", open = "w"), quote= FALSE);
# argList <- as.list(read.table(text = transposedArgString, stringsAsFactors = FALSE));
# 
# # reading data
# fc <- as.matrix(read.table(argList$fcFile, check.names=FALSE, sep="\t")); # training fold change expression matrix, rows are samples, cols are genes
# isPerturbed <- as.matrix(read.table(argList$isPerturbedFile, check.names=FALSE)); # training logical matrix of perturbation, rows are samples, cols are genes
# tfName <- read.table(argList$tfNameFile, stringsAsFactors = FALSE)[[1]]; # get TF names, each row is a TF

# # preparing training data
# tfLevel <- matrix(1, nrow(fc), length(tfName), dimnames = list(rownames(fc), tfName));
# availableTfName <- intersect(tfName, colnames(fc)); # names of TFs available in the expression data
# unavailableTfName <- setdiff(tfName, colnames(fc));
# if (length(unavailableTfName) > 0) warning(paste("data unavailable for TF", paste(unavailableTfName, collapse = ", ")));
# tfLevel[, availableTfName] <- fc[, availableTfName];
# tgtLevel <- fc;
# tgtLevel[isPerturbed] <- NA;  # masking perturbed entries in response variables from training

# # preparing regMat
# nTf <- ncol(tfLevel);
# nTgt <- ncol(tgtLevel);
# tgtName <- colnames(tgtLevel);
# regMat <- matrix(FALSE, nTf, nTgt, dimnames = list(tfName, tgtName));
# regMat[availableTfName, ] <- TRUE; # only allow a target to be regulated by TFs whose expression levels are avaiable in the data
# regMat[cbind(availableTfName, availableTfName)] <- FALSE; # disallow autoregulation

# if (!is.null(argList$saveTo)) saveProcessTo <- paste(argList$saveTo, ".bartProcess.RData", sep = "") else saveProcessTo <- NULL;
# if (!is.null(argList$useMpi))
# 	if (argList$useMpi)
# 		mpiComm <- 1 else mpiComm <- NULL;
# if (!is.null(argList$mpiBlockSize))
# 	if (argList$mpiBlockSize)
# 		mpiBlockSize <- argList$mpiBlockSize else mpiBlockSize <- 1;

# #try(Sys.setenv(OMPI_MCA_btl_tcp_if_include="eth0")); # otherwise mpi freezes, not sure why. Please ask IT support for details
# result <- getBartNetwork(
# 	tgtLevel = tgtLevel,
# 	tfLevel = tfLevel,
# 	regMat = regMat,
# 	mpiComm = mpiComm,
# 	blockSize = mpiBlockSize,
# 	saveTo = saveProcessTo
# );

# if (!is.null(argList$saveTo)) {
# 	save(result, file = paste(argList$saveTo, ".bartResult.RData", sep = ""));
# 	write.table(result$regScore, file = argList$saveTo, sep = "\t");
# }
