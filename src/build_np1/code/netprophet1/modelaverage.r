compute.model.average <- function(B,P,params)
{
	region_enrichment <- c(params[3],params[2],params[4],params[1],params[5],params[2],params[6],params[1],0)
	combined.adj <- (abs(B)+params[8]) * (abs(P)+params[7])
	Q.test.I <- combined.adj; Q.test.I[which(B <= 0)] <- 0; Q.test.I[which(P <= 0)] <- 0;
	Q.test.II <- combined.adj; Q.test.II[which(B >= 0)] <- 0; Q.test.II[which(P <= 0)] <- 0;
	Q.test.III <- combined.adj; Q.test.III[which(B >= 0)] <- 0; Q.test.III[which(P >= 0)] <- 0;
	Q.test.IV <- combined.adj; Q.test.IV[which(B <= 0)] <- 0; Q.test.IV[which(P >= 0)] <- 0;
	Q.test.PDE <- combined.adj; Q.test.PDE[which(B != 0)] <- 0; Q.test.PDE[which(P <= 0)] <- 0;
	Q.test.NDE <- combined.adj; Q.test.NDE[which(B != 0)] <- 0; Q.test.NDE[which(P >= 0)] <- 0;
	Q.test.PB <- combined.adj; Q.test.PB[which(B <= 0)] <- 0; Q.test.PB[which(P != 0)] <- 0;
	Q.test.NB <- combined.adj; Q.test.NB[which(B >= 0)] <- 0; Q.test.NB[which(P != 0)] <- 0;
	Q.test <- region_enrichment[1] * Q.test.I + region_enrichment[2] * Q.test.PDE + region_enrichment[3] * Q.test.II + region_enrichment[4] * Q.test.PB + region_enrichment[5] * Q.test.III + region_enrichment[6] * Q.test.NDE + region_enrichment[7] * Q.test.IV + region_enrichment[8] * Q.test.NB
	Q.test
}

compute.model.average.new <- function(B,D,params)
{
	# params
	# (I,II,III,IV,B,D,Cb,Cd)
	Bs <- abs(B) + params[7]
	Ds <- abs(D) + params[8]
	M <- matrix(0,nrow=dim(B)[1],ncol=dim(B)[2])
	rindx <- intersect(which(B>0),which(D>0))								# I
	M[rindx] <- abs(Bs[rindx]) * abs(Ds[rindx]) * params[1] 
	rindx <- intersect(which(B<0),which(D>0))								# II
	M[rindx] <- abs(Bs[rindx]) * abs(Ds[rindx]) * params[2] 
	rindx <- intersect(which(B<0),which(D<0))								# III
	M[rindx] <- abs(Bs[rindx]) * abs(Ds[rindx]) * params[3] 
	rindx <- intersect(which(B>0),which(D<0))								# IV
	M[rindx] <- abs(Bs[rindx]) * abs(Ds[rindx]) * params[4] 
	rindx <- intersect(which(B!=0),which(D==0))							# B
	M[rindx] <- abs(Bs[rindx]) * abs(Ds[rindx]) * params[5] 
	rindx <- intersect(which(B==0),which(D!=0))							# D
	M[rindx] <- abs(Bs[rindx]) * abs(Ds[rindx]) * params[6] 
	M
}

plot.model.average.ranking <- function(B,P,params,n=10000,fname)
{
  region_enrichment <- c(params[3],params[2],params[4],params[1],params[5],params[2],params[6],params[1],0)
  combined.adj <- (abs(B)+params[8]) * (abs(P)+params[7])
  Q.test.I <- combined.adj; Q.test.I[which(B <= 0)] <- 0; Q.test.I[which(P <= 0)] <- 0;
  Q.test.II <- combined.adj; Q.test.II[which(B >= 0)] <- 0; Q.test.II[which(P <= 0)] <- 0;
  Q.test.III <- combined.adj; Q.test.III[which(B >= 0)] <- 0; Q.test.III[which(P >= 0)] <- 0;
  Q.test.IV <- combined.adj; Q.test.IV[which(B <= 0)] <- 0; Q.test.IV[which(P >= 0)] <- 0;
  Q.test.PDE <- combined.adj; Q.test.PDE[which(B != 0)] <- 0; Q.test.PDE[which(P <= 0)] <- 0;
  Q.test.NDE <- combined.adj; Q.test.NDE[which(B != 0)] <- 0; Q.test.NDE[which(P >= 0)] <- 0;
  Q.test.PB <- combined.adj; Q.test.PB[which(B <= 0)] <- 0; Q.test.PB[which(P != 0)] <- 0;
  Q.test.NB <- combined.adj; Q.test.NB[which(B >= 0)] <- 0; Q.test.NB[which(P != 0)] <- 0;
  Q.test <- region_enrichment[1] * Q.test.I + region_enrichment[2] * Q.test.PDE + region_enrichment[3] * Q.test.II + region_enrichment[4] * Q.test.PB + region_enrichment[5] * Q.test.III + region_enrichment[6] * Q.test.NDE + region_enrichment[7] * Q.test.IV + region_enrichment[8] * Q.test.NB

  regions <- matrix(0,ncol=dim(B)[2],nrow=dim(B)[1])
  regions[intersect(which(B > 0),which(P > 0))] <- 1	#IV
  regions[intersect(which(B < 0),which(P > 0))] <- 2	#III
  regions[intersect(which(B < 0),which(P < 0))] <- 3	#II
  regions[intersect(which(B > 0),which(P < 0))] <- 4	#I
  regions[intersect(which(B == 0),which(P != 0))] <- 5	#P
  regions[intersect(which(B != 0),which(P == 0))] <- 6	#B

  r.freq <- matrix(0,ncol=n,nrow=6)
  r.nzi <- which(regions!=0)
  q.nz <- Q.test[r.nzi]
  r.nz <- regions[r.nzi]
  r.rnk <- r.nz[sort.list(q.nz,decreasing=TRUE)]
  for (i in 2:n)
  {
    r.freq[,i] <- r.freq[,i-1]
    r.freq[r.rnk[i],i] <- r.freq[r.rnk[i],i] + 1
  }
  r.frac <- matrix(0,ncol=n,nrow=6)
  r.frac[1,] <- r.freq[1,]
  for (i in 2:6) 
  { 
    r.frac[i,] <- r.freq[i,] + r.frac[i-1,] 
  }
  r.frac <- t(t(r.frac) / (seq(1,n)-1))
	jpeg(filename=fname,width=1800,height=1800,quality=100,res=300)
  plot(seq(1,n)-1,rep(1,times=n),ylim=c(0,1),type='n',xlab="Interaction rank",ylab="Cumulative representation")
  xx <- c(seq(1,n),seq(n,1))-1
	
  cols <- c(rgb(0.7,0,0.7),rgb(0.8,0.3,0.3),rgb(1,0.6,0.2),rgb(0.3,0.8,0.3),rgb(0.6,1,1),rgb(0.2,0.6,1))

  yy <- c(rep(0,times=n),rev(r.frac[1,]))
  polygon(xx,yy,col=cols[1])
  for (i in 1:5)
  {
    yy <- c(r.frac[i,],rev(r.frac[i+1,]))
    polygon(xx,yy,col=cols[i+1])
  }
	dev.off()
}

