source("modelaverage.r")

lasso_component <- lasso_component / max(abs(lasso_component))

indx <- which(de_component>0)
de_component[indx] <- de_component[indx] - min(abs(de_component[indx]))
indx <- which(de_component<0)
de_component[indx] <- de_component[indx] + min(abs(de_component[indx]))
de_component <- de_component / max(abs(de_component))
## Untrained Parameters
#combinedAdjMtr <- compute.model.average.new(betas,de_component,c(1,1,1,1,1,1,0.01,0.01))
## Trained Parameters
combinedAdjMtr <- compute.model.average.new(lasso_component,de_component,c(3,1,1,1,1,2,0.1,0.01))

# write.table(combinedAdjMtr,file.path(outputDirectory,combinedAdjMtrFileName),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
write.table(combinedAdjMtr,combinedAdjMtrFileName,row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
