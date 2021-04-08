rank.decreasing <- function(x) length(x) + 1 - rank(x,ties.method="first")

buildInteractions <- function(M,B,D,net_genes,net_regulators,GIDS,RIDS,cutoff)
{
	M.rank <- matrix(rank.decreasing(M),nrow=dim(M)[1],ncol=dim(M)[2])
	net_genes.indx <- match(net_genes,GIDS)
	net_regulators.indx <- match(net_regulators,RIDS)
	interactions <- c()
	for (rindx in net_regulators.indx)
	{
		tindx <- which(M[rindx,]>cutoff)
		tindx <- intersect(tindx,net_genes.indx)
		if (length(tindx)>0)
		{
			interactions <- rbind(interactions,cbind(rep(RIDS[rindx],times=length(tindx)),GIDS[tindx],M[rindx,tindx],sign(D[rindx,tindx]),sign(B[rindx,tindx]),(sign(D[rindx,tindx])+sign(B[rindx,tindx]))/2,M.rank[rindx,tindx]))
		}
	}
	
	colnames(interactions) <- c("REGULATOR","TARGET","SCORE","DSIGN","BSIGN","CSIGN","RANK")
	interactions
}

M <- as.matrix(read.table(file.path(outputDirectory, combinedAdjMtrFileName)))
B <- as.matrix(read.table(file.path(outputDirectory, lassoAdjMtrFileName)))
D <- as.matrix(read.table(differentialExpressionMatrixFile))
RIDS <- as.matrix(read.table(regulatorGeneNamesFileName))
GIDS <- as.matrix(read.table(targetGeneNamesFileName))

interactions <- buildInteractions(M,B,D,GIDS,RIDS,GIDS,RIDS,0)
write.table(as.table(interactions),file=file.path(outputDirectory,combinedAdjLstFileName),row.names=FALSE,quote=FALSE,sep='\t')

