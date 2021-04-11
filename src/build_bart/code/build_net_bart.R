generate_bart_net = function(p_in_expr_target
                             , p_in_expr_reg
                             , fname_bart
                             , p_out_dir
                             , flag_slurm
                             , seed
                             , p_src_code
                             , nbr_rmpi_slave){
    
    
    # p_src_code='/scratch/mblab/dabid/netprophet/code_netprophet3.0/'
    # p_in_expr_target='netprophet/net_debug/data/expr_target_indexed'
    # p_in_expr_reg='netprophet/net_debug/data/expr_reg_indexed'
    # fname_bart='net_bart.tsv'
    # p_out_dir='netprophet/net_debug/'
    # flag_slurm='ON'
    # seed=747
        
    # read data from files 
    # p_in_expr_target='netprophet/net_in/all_kem_expr_6112_1485_indexed'
    # p_in_expr_reg='netprophet/net_in/all_kem_expr_reg_313_1485_indexed'
    
    df_expr_target = read.csv(p_in_expr_target, header=TRUE, row.names=1, sep="\t")
    df_expr_reg = read.csv(p_in_expr_reg, header=TRUE, row.names=1, sep="\t")
    l_in_target = as.factor(rownames(df_expr_target))
    l_in_reg = as.factor(rownames(df_expr_reg))
    l_in_sample = as.factor(colnames(df_expr_target))
    # rownames(df_expr_target) = NULL
    # colnames(df_expr_target) = NULL
    # rownames(df_expr_reg) = NULL
    # colnames(df_expr_reg) = NULL
    
    # transform from log(fc) to fc
    df_expr_target = 2**df_expr_target
    df_expr_reg = 2**df_expr_reg
    
    # generate intermediate files (allowed and perturbed)
    df_allowed_perturbed = generate_allowed_perturbed_matrices(l_in_target, l_in_reg, l_in_sample, NULL, p_src_code)
    df_allowed = sapply(as.data.frame(as.matrix(df_allowed_perturbed[[1]])), as.logical)
    df_perturbed = sapply(as.data.frame(as.matrix(df_allowed_perturbed[[2]])), as.logical)
    
    # masking pertubed entries in response from training
    df_expr_target[df_perturbed] = NA

    df_bart_net = getBartNetwork(tgtLevel=t(as.matrix(df_expr_target))
                                     , tfLevel=t(as.matrix(df_expr_reg))
                                     , regMat=df_allowed
                                     , mpiComm=1
                                     , blockSize=nbr_rmpi_slave
                                    )
    
    # if output directory doesn't exist, create it
    ifelse(!dir.exists(file.path(p_out_dir))
           , dir.create(file.path(p_out_dir), showWarnings=FALSE)
           , FALSE)
    # save(df_bart_net, file = paste(p_out_dir, ".bartResult.RData", sep = ""));
    write.table(df_bart_net$regScore
                , file.path(p_out_dir, fname_bart)
                , row.names=l_in_reg
                , col.names=l_in_target
                , quote=FALSE
                , sep="\t")
                        
}

if (sys.nframe() == 0){
    # =========================================== #
    # |       *** Install packages ***          | #
    # =========================================== #
    library("optparse")
    
    # =========================================== #
    # |         **** Parse Arguments ****       | #
    # =========================================== #
    p_in_expr_target = make_option(c("--p_in_expr_target"), type="character", help='input - path of expression of target genes', default=NULL)
    p_in_expr_reg = make_option(c("--p_in_expr_reg"), type="character", help="input - path of expression of regulators", default=NULL)
    fname_bart = make_option(c("--fname_bart"), type="character", default=NULL, help="output - path of generated bart network")
    p_out_dir = make_option(c("--p_out_dir"), type="character", default=NULL, help="output - path of output directory for results")
    flag_slurm = make_option(c("--flag_slurm"), type="character", default="OFF", help="ON or OFF for MPI run")
    seed = make_option(c("--seed"), type="integer", default=747, help="seed for reproducibility")
    p_src_code = make_option(c("--p_src_code"), type="character", default=NULL, help="path of the source code")
    nbr_rmpi_slave = make_option(c("--nbr_rmpi_slave"), type="integer", default=2)    
    opt_parser = OptionParser(option_list=list(p_in_expr_target, p_in_expr_reg, fname_bart, p_out_dir, flag_slurm, seed, p_src_code, nbr_rmpi_slave))
    opt = parse_args(opt_parser)
    
    if (is.null(opt$p_in_expr_target) || is.null(opt$p_in_expr_reg) || is.null(opt$fname_bart) || is.null(opt$p_out_dir) || is.null(opt$flag_slurm) || is.null(opt$p_src_code)
       )
    {
        print_help(opt_parser)
        stop("Arguments p_in_expr_target, p_in_expr_reg, fname_bart, p_out_dir, flag_slurm are mandatory")
    }
    
    # load local R
    source(paste(opt$p_src_code, "src/build_bart/code/netprophet2/build_bart_network.r", sep=""))
    source(paste(opt$p_src_code, "src/build_np1/code/prepare_data_generate_allowed_perturbed_and_scale_normalize.R", sep=""))
             
    quit(status=generate_bart_net(p_in_expr_target=opt$p_in_expr_target
                                  , p_in_expr_reg=opt$p_in_expr_reg
                                  , fname_bart=opt$fname_bart
                                  , p_out_dir=opt$p_out_dir
                                  , flag_slurm=opt$flag_slurm
                                  , seed=opt$seed
                                  , p_src_code=opt$p_src_code
                                  , nbr_rmpi_slave=opt$nbr_rmpi_slave
                                 ))
}
