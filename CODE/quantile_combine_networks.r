#!/usr/bin/R
# Combine mutliple networks using Quantile Combine method (ref: Hien's Thesis)

# source("/home/mblab/ykang/proj_hien_thesis/hliowThesisScript.R")
quantileCombine <- function(scoreList, plot = FALSE, verbose = FALSE, returnScoreListOfAllRepresentation = FALSE) {
	# quantile combination of scores: score can either be specified in list(...) or scoreList
	# plot: whether to plot the histogram of the distribution of each score
	# returnScoreListOfAllRepresentation: whether to return the combination using every representative distribution or only the final combined score
	require("matrixStats");
	if (missing(scoreList)) scoreList <- list(...); # put scores in scoreList
	rawScoreMat <- do.call("cbind", lapply(scoreList, c)); # each score source as a column in the matrix
	absScoreList <- lapply(scoreList, abs);
	nonNaCount <- colSums(!is.na(rawScoreMat)); # number of non NA entries of each member
	combinedScoreMat <- rawScoreMat * NA; # matrix to store combined scores, each column is a combination using one representative distribution
	if (verbose) cat("Ranking scores ... ");
	rankMat <- apply(abs(rawScoreMat), 2, rank);  # rawScoreMat converted to ranks
	rankMat[is.na(rawScoreMat)] <- NA;
	if (verbose) cat("DONE\n");
	if (length(unique(nonNaCount)) == 1) { # if every member has the same amount of non NA entries, invoke sort and rank functions to combine scores really fast
		if (verbose) cat("Sorting scores ... ");
		sorted <- apply(abs(rawScoreMat), 2, sort);
		ix1 <- head(rep(sequence(nrow(sorted)), each = 2), -1);
		ix2 <- tail(rep(sequence(nrow(sorted)), each = 2), -1);
		sortedWithHalfwayValue <- (sorted[ix1, , drop = FALSE] + sorted[ix2, , drop = FALSE]) / 2; # sorted absScores with halfway values
		for (i in sequence(ncol(combinedScoreMat))) {
			# if (verbose) cat(i, "\n");
			convertedMat <- rankMat * NA; # matrix to store scores converted to the representative quantiles
			convertedMat[!is.na(rawScoreMat)] <- sortedWithHalfwayValue[c(rankMat[!is.na(rawScoreMat)]) * 2 - 1, i]; # quantile conversion
			combinedScoreMat[, i] <- rowMeans(convertedMat * sign(rawScoreMat), na.rm = TRUE); # averaging the converted scores
		}
	} else { # if members have different amounts of non NA entries, cannot invoke sort and rank, quantiles need to be estimated 
		if (verbose) cat("Generating empirical cdf... ");
		qList <- lapply(scoreList, ecdf); # list if ecdf
		portionMat <- (rankMat - 1) / rep(colMaxs(rankMat, na.rm = TRUE) - 1, each = nrow(rankMat)); # rankMat scaled to [0, 1]
		portionMatrix[naIx] <- 0;
		convertedMat <- rankMat * NA; # matrix to store scores converted to the representative quantiles
		for (estimatedDistribution in qList) {
			if (plot) hist(estimatedDistribution);
			if (verbose) print(dim(convertedMat));
			convertedMat[!is.na(rawScoreMat)] <- quantile(estimatedDistribution, c(portionMatrix[!is.na(rawScoreMat)]), na.rm = TRUE); # quantile conversion
			combinedScoreMat[, i] <- rowMeans(convertedMat * sign(rawScoreMat), na.rm = TRUE); # averaging the converted scores
		}
	}
	if (verbose) cat("DONE\n");
	if (returnScoreListOfAllRepresentation) { # reshape combinedScoreMat to a list of combined scores, each is bases on one representative distribution
		combinedScoreList <- list();
		for (i in sequence(ncol(combinedScoreMat))) {
			combinedScoreList[[i]] <- scoreList[[1]] * NA;
			combinedScoreList[[i]][] <- combinedScoreMat[, i];
		}
		names(combinedScoreList) <- names(scoreList);
		combinedScoreList;
	} else { # combine combined scores based on different representative distributions
		combinedRankMat <- apply(abs(combinedScoreMat), 2, rank) * sign(combinedScoreMat);
		combinedRankMat[is.na(combinedScoreMat)] <- NA;
		score <- scoreList[[1]] * NA; # to store final combined score
		score[] <- rowMeans(combinedRankMat);
		score;
	}
}

# main script
args <- commandArgs(trailingOnly=TRUE)
data <- list()
cat("Loading network ... ")
for (i in 1:(length(args)-1)) {
	data[length(data)+1] <- list(as.matrix(read.table(toString(args[i]))))
}
cat("DONE\n")
combined <- quantileCombine(data, verbose=TRUE, returnScoreListOfAllRepresentation=TRUE)
write.table(combined[[1]],args[length(args)],row.names=FALSE,col.names=FALSE,quote=FALSE)
