# The code of this module is gonna be seen by all nodes, and hence all used variables here can be called in other function that will be run in the nodes..

source(paste(p_src_code, "src/build_np1/code/netprophet1/global.lars.regulators.r", sep=""))

# setting the seed, and create the folds of cross validation has to be here
set.seed(seed)
all.folds <- cv.folds(dim(df_expr_target)[2], nbr_cv_fold)

# submit every fold of CV to a different node
fold <- mpi.comm.rank()
cat("fold number", fold, "\n")
# TODO to describe the cv.obj here
cv.obj <- lars.multi.cv.singlefold(df_expr_target
                                   , df_expr_reg
                                   , df_perturbed
                                   , df_prior
                                   , df_allowed
                                   , nbr_cv_fold
                                   , seed
                                   , fold)

# save(cv.obj,file=file.path(p_out_dir, paste("cv.",fold,sep="")))

get_fold_nbr <- function()
{
	fold
}

