build_net_np1 = function(p_in_expr_target
                         , p_in_expr_reg
                         , flag_global_shrinkage
                         , flag_local_shrinkage
                         , lasso_nbr_fold
                         , flag_microarray
                         , p_in_net_de
                         , seed
                         , p_out_dir
                         , f_out_name_lasso
                         , f_out_name_np1
                         , p_src_code
                         , flag_slurm){
    # load libraries
    source(paste(p_src_code, 'src/build_np1/code/build_net_lasso.R', sep=''))
    source(paste(p_src_code, 'src/build_np1/code/netprophet1/modelaverage.r', sep=''))
    # get lasso
    df_net_lasso = generate_lasso_net(p_in_expr_target=p_in_expr_target
                                      , p_in_expr_reg=p_in_expr_reg
                                      , flag_global_shrinkage=flag_global_shrinkage
                                      , flag_local_shrinkage=flag_local_shrinkage
                                      , p_out_dir=p_out_dir
                                      , fname_lasso=f_out_name_lasso
                                      , flag_parallel=flag_slurm
                                      , seed=seed
                                      , nbr_cv_fold=lasso_nbr_fold
                                      , p_src_code=p_src_code
                                      , flag_microarray=flag_microarray)
    # combine lasso and de
    # normalize lasso network
    df_net_lasso <- df_net_lasso / max(abs(df_net_lasso))
    
    # read de network and normalize it
    df_net_de = as.matrix(read.csv(p_in_net_de, header=TRUE, row.names=1, sep='\t'))
    indx <- which(df_net_de>0)
    df_net_de[indx] <- df_net_de[indx] - min(abs(df_net_de[indx]))
    indx <- which(df_net_de<0)
    df_net_de[indx] <- df_net_de[indx] + min(abs(df_net_de[indx]))
    df_net_de <- df_net_de / max(abs(df_net_de))
    ## Untrained Parameters
    #df_net_np1 <- compute.model.average.new(df_net_lasso,df_net_de,c(1,1,1,1,1,1,0.01,0.01))
    ## Trained Parameters
    
    # combine lasso and de
    df_net_np1 = compute.model.average.new(df_net_lasso,df_net_de,c(3,1,1,1,1,2,0.1,0.01))
    
    write.table(df_net_np1
                ,paste(p_out_dir, f_out_name_np1, sep='')
                ,row.names=rownames(df_net_lasso)
                ,col.names=colnames(df_net_lasso)
                ,quote=FALSE
                ,sep='\t')

}

if (sys.nframe() == 0){
    library("optparse")
    
    opt_parser = OptionParser(option_list=list(
        p_in_expr_target = make_option(c("--p_in_expr_target"), type="character")
        , p_in_expr_reg = make_option(c("--p_in_expr_reg"), type="character")
        , flag_global_shrinkage = make_option(c("--flag_global_shrinkage"), type="character")
        , flag_local_shrinkage = make_option(c("--flag_local_shrinkage"), type="character")
        , lasso_nbr_fold = make_option(c("--lasso_nbr_fold"), type="character")
        , flag_microarray = make_option(c("--flag_microarray"), type="character")
        , p_in_net_de = make_option(c("--p_in_net_de"), type="character")
        , seed = make_option(c("--seed"), type="integer")
        , p_out_dir = make_option(c("--p_out_dir"), type="character")
        , f_out_name_lasso = make_option(c("--f_out_name_lasso"), type="character")
        , f_out_name_np1 = make_option(c("--f_out_name_np1"), type="character")
        , p_src_code = make_option(c("--p_src_code"), type="character")
        , flag_slurm = make_option(c("--flag_slurm"), type="character")
    ))
    
    opt = parse_args(opt_parser, positional_arguments=TRUE)$options
    
    build_net_np1(p_in_expr_target=opt$p_in_expr_target
                              , p_in_expr_reg=opt$p_in_expr_reg
                              , flag_global_shrinkage=opt$flag_global_shrinkage
                              , flag_local_shrinkage=opt$flag_local_shrinkage
                              , lasso_nbr_fold=opt$lasso_nbr_fold
                              , flag_microarray=opt$flag_microarray
                              , p_in_net_de=opt$p_in_net_de
                              , seed=opt$seed
                              , p_out_dir=opt$p_out_dir
                              , f_out_name_lasso=opt$f_out_name_lasso
                              , f_out_name_np1=opt$f_out_name_np1
                              , p_src_code=opt$p_src_code
                              , flag_slurm=opt$flag_slurm)
}