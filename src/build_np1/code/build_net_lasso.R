generate_lasso_net = function(p_in_expr_target
                              , p_in_expr_reg
                              , p_out_dir
                              , flag_global_shrinkage
                              , flag_local_shrinkage
                              , fname_lasso
                              , flag_parallel
                              , seed
                              , nbr_cv_fold
                              , p_src_code
                              , flag_microarray
                              , p_in_target="NONE"
                              , p_in_reg="NONE"
                              , p_in_sample="NONE"
                              , greenfield_method="OFF"
                              , p_lasso_greenfield=NULL
                              , idx_reg_greenfield=NULL
                             ){
  # ================================================================================================= #
  # |                            **** Load local R libraries ****                                   | #
  # | These libraries are developed by the lab, including:                                          | #
  # | 1. global.lars.regulators.r where functions for creating LASSO with global/local shrinkage,   | #
  # | and parallel/sequential implementation.                                                       | #
  # | 2. prepare_data_generate_allowed_perturbed_and_scale_normalize.R where functions for creating | #
  # |    allowed/perturbed matrices, and                                                            | #
  # | scaling/normalizing the expression matrices                                                   | #
  # ================================================================================================= #
  source(paste(p_src_code, "src/build_np1/code/netprophet1/global.lars.regulators.r", sep=""))  # for generating lasso network
  source(paste(p_src_code, "src/build_np1/code/prepare_data_generate_allowed_perturbed_and_scale_normalize.R", sep=""))  # for preparing data
  source(paste(p_src_code, "src/build_np1/code/helper.R", sep=""))  # for helper functions
  
  # ================================================================================================= #
  # |                                   **** Read Input data  ****                                  | #
  # ================================================================================================= #
  if (p_in_target == "NONE" || p_in_reg == "NONE" || p_in_sample == "NONE"){
      df_expr_target = read.csv(p_in_expr_target, header=TRUE, row.names=1, sep="\t")
      df_expr_reg = read.csv(p_in_expr_reg, header=TRUE, row.names=1, sep="\t")
      l_in_target = as.factor(rownames(df_expr_target))
      l_in_reg = as.factor(rownames(df_expr_reg))
      l_in_sample = as.factor(colnames(df_expr_target))
  } else{
      l_in_target = read.csv(p_in_target, header=FALSE)[[1]]
      l_in_reg = read.csv(p_in_reg, header=FALSE)[[1]]
      l_in_sample = read.csv(p_in_sample, header=FALSE)[[1]]
      df_expr_target = read.csv(p_in_expr_target, header=FALSE, sep="\t")
      df_expr_reg = read.csv(p_in_expr_reg, header=FALSE, sep="\t") 
  }
  
   
  # ================================================================================================= #
  # |                                      **** Prepare data  ****                                  | #
  # | 1. generate the allowed and perturbed matrices                                                | #
  # | 2. scale and normalize expression matrices for regulators and target genes           
  # ================================================================================================= #
   # generate intermediate files (allowed and perturbed)
   df_allowed_perturbed = generate_allowed_perturbed_matrices(l_in_target, l_in_reg, l_in_sample, NULL, p_src_code)
   df_allowed = as.matrix(df_allowed_perturbed[[1]])
   df_perturbed = as.matrix(df_allowed_perturbed[[2]])
   # scale and normalize expression matrices for regulators and target genes
   df_expr_target = scale_normalize_expr_matrix(df_expr_target, flag_microarray)
   df_expr_target = as.matrix(df_expr_target)
   df_expr_reg = scale_normalize_expr_matrix(df_expr_reg, flag_microarray)
   df_expr_reg = as.matrix(df_expr_reg)
   
   df_prior = matrix(1,ncol=dim(df_expr_target)[1] ,nrow=dim(df_expr_reg)[1] )
   
   # ================================================================================================= #
   # |                                **** Build LASSO NETWORK ****                                  | #
   # | To build lasso network, there are multiple options of doing that                              | #
   # | 1. with local shrinkage parameter, a lambda is estimated for every target gene.               | #
   # | 2. with global shrinkage parameter, one single lambda is estimated for all target genes,      | #
   # |    which can be done in a distributed faster fashion using Rmpi package or sequentially       | #
   # | 3. Using greenfield method by training every time we remove a TF                              | #
   # | 4. Using greenfiled method by using the trained parameter in the first round.                 | #
   # ================================================================================================= #
   

   
   # ------------------------------------------------------------------------------ #
   # |                     **** PARALLEL GLOBAL Shrinkage ****                    | #
   # ------------------------------------------------------------------------------ #
   # if I don't NULL the indexes, I get an error in the generation of LASSO
   if (p_in_target == "NONE" || p_in_reg == "NONE" || p_in_sample == "NONE"){
       # if I don't NULL the indexes, I get an error in the generation of LASSO
       rownames(df_expr_target) = NULL
       colnames(df_expr_target) = NULL
       rownames(df_expr_reg) = NULL
       colnames(df_expr_reg) = NULL
     } 
    
   if((flag_global_shrinkage == "ON" && flag_parallel == "ON") || greenfield_method == "ONE"){
     df_lasso_net = create_lasso_global_shrinkage_parallel(df_expr_target
                                                             , df_expr_reg
                                                             , df_allowed
                                                             , df_perturbed
                                                             , df_prior
                                                             , seed
                                                             , p_out_dir
                                                             , nbr_cv_fold
                                                             , p_src_code)
   # ------------------------------------------------------------------------------ #
   # |              **** Rerank LASSO With Greenfield method 1 ****               | #
   # ------------------------------------------------------------------------------ #
     if (greenfield_method == "ONE"){
       mse_full_reg = calculate_mse(l_in_target, l_in_sample, df_expr_target, df_expr_reg, df_lasso_net)
       net_lasso_greenfield = list()  
       for (idx_reg in seq(1, length(l_in_reg), 1)){  # for loop for removing a TF at a time
         df_lasso_net_minus_reg = df_lasso_net[-idx_reg, ]
         # calculate MSE for removing one regulator at a time for df_lasso_net_minus_reg
         mse_minus_reg = calculate_mse(l_in_target, l_in_sample, df_expr_target, df_expr_reg[-idx_reg, ], df_lasso_net_minus_reg)
         # calculate the scors 1 - variance of full/vairance of minus
         net_lasso_greenfield[[idx_reg]] = mapply("-"
                                                  , matrix(1, ncol=1, nrow=length(l_in_target))
                                                  ,(mapply("/", mse_full_reg, mse_minus_reg, SIMPLIFY=FALSE))
                                                  , SIMPLIFY=FALSE)
       }
       df_lasso_net = data.frame(matrix(unlist(net_lasso_greenfield), ncol=max(lengths(net_lasso_greenfield)), byrow=TRUE))
     }
   # ------------------------------------------------------------------------------ #
   # |              **** Rerank LASSO With Greenfield method 2 ****               | #
   # ------------------------------------------------------------------------------ #
   } else if (greenfield_method == "TWO"){
       df_lasso_net = read.csv(p_lasso_greenfield, header=FALSE, sep=" ")
       mse_full_reg = calculate_mse(l_in_target, l_in_sample, df_expr_target, df_expr_reg, df_lasso_net)
       # generate the allowed and perturbed matrices for reg minus
       df_allowed_perturbed_minus_reg = generate_allowed_perturbed_matrices(l_in_target, l_in_reg[-idx_reg_greenfield], l_in_sample, NULL, p_src_code)
       df_allowed_minus_reg = as.matrix(df_allowed_perturbed_minus_reg[[1]])
       df_perturbed_minus_reg = as.matrix(df_allowed_perturbed_minus_reg[[2]])
       # create the df_expr_reg_minus_reg
       df_expr_reg_minus_reg = df_expr_reg[-idx_reg_greenfield, ]
       df_prior_minus_reg = matrix(1,ncol=dim(df_expr_target)[1] ,nrow=dim(df_expr_reg_minus_reg)[1])
       
       df_lasso_net_minus_reg = create_lasso_global_shrinkage_parallel(df_expr_target
                                                             , df_expr_reg_minus_reg
                                                             , df_allowed_minus_reg
                                                             , df_perturbed_minus_reg
                                                             , df_prior_minus_reg
                                                             , seed
                                                             , p_out_dir
                                                             , nbr_cv_fold
                                                             , p_src_code)
       # calculate MSE for removing one regulator at a time for df_lasso_net_minus_reg
       mse_minus_reg = calculate_mse(l_in_target, l_in_sample, df_expr_target, df_expr_reg_minus_reg, df_lasso_net_minus_reg)
       # calculate the scors 1 - variance of full/vairance of minus This is network for this TF
       df_lasso_net = mapply("-"
                            , matrix(1, ncol=1, nrow=length(l_in_target))
                            ,(mapply("/", mse_full_reg, mse_minus_reg, SIMPLIFY=FALSE))
                            , SIMPLIFY=FALSE)
    
   # ------------------------------------------------------------------------------ #
   # |                     **** SEQUENTIAL GLOBAL Shrinkage ****                  | #
   # ------------------------------------------------------------------------------ #   
   } else if(flag_global_shrinkage == "ON" && flag_parallel == "OFF" )
     {
       df_lasso_net = lars.multi.optimize(df_expr_target
                                          , df_expr_reg
                                          , df_perturbed
                                          , df_prior
                                          , df_allowed)[[1]]
    # ------------------------------------------------------------------------------ #
    # |                          **** LOCAL Shrinkage ****                         | #
    # ------------------------------------------------------------------------------ #
     } else if (flag_local_shrinkage == "ON"){
        ## Skip regression on some genes
        cat("Gene skipped count:")
        skip_gen <- rep(0, dim(df_expr_target)[1])
        for (i in 1:dim(df_expr_target)[1]) {
          if (sum(df_expr_target[i,] != 0) < dim(df_expr_target)[2]/10+1) {
              skip_gen[i] = 1
          }
        }
        cat(length(which(skip_gen == 1)), "\nRegulator skipped:")
        skip_reg <- rep(0, dim(df_expr_reg)[1])
        for (i in 1:dim(df_expr_reg)[1]) {
          if (sum(df_expr_reg[i,] != 0) < 1) {
              skip_reg[i] = 1
              cat(i, ",")
          }
        }
        cat("\n")
        df_lasso_net = lars.local(df_expr_target
                                 , df_expr_reg
                                 , df_perturbed
                                 , df_prior
                                 , df_allowed
                                 , skip_reg
                                 , skip_gen)
       }
   # ================================================================================================= #
   # |                            **** END Build LASSO NETWORK ****                                  | #
   # ================================================================================================= #
   
   
   # ================================================================================================= #
   # |                              **** Write LASSO NETWORK ****                                    | #
   # | write lasso network into tsv file, if the output directory doesn't exist create it            | #
   # ================================================================================================= #
   # if output directory doesn't exist, create it
   ifelse(!dir.exists(file.path(p_out_dir))
          , dir.create(file.path(p_out_dir), showWarnings=FALSE)
          , FALSE
   )
   # write lasso network
   write.table(df_lasso_net
               , file.path(p_out_dir, fname_lasso)
               , row.names=l_in_reg
               , col.names=l_in_target
               , quote=FALSE
               , sep="\t"
               )
   # ================================================================================================= #
   # |                            **** END Write LASSO NETWORK ****                                  | #
   # ================================================================================================= #
   
   if (is.loaded("mpi_initialize")) {
     mpi.quit
   }
   rownames(df_lasso_net) = l_in_reg
   colnames(df_lasso_net) = l_in_target
   df_lasso_net
}

if (sys.nframe() == 0){
  # =========================================== #
  # |       *** Install packages ***          | #
  # =========================================== #
  library("optparse")
  
  # =========================================== #
  # |         **** Parse Arguments ****       | #
  # =========================================== #
  p_in_expr_target = make_option(c("--p_in_expr_target"), type="character", help='input - path of expression of target genes')
  p_in_expr_reg = make_option(c("--p_in_expr_reg"), type="character", help="input - path of expression of regulators")
  p_in_target = make_option(c("--p_in_target"), type="character", default="NONE", help="input - path of list of target gene ids")
  p_in_reg = make_option(c("--p_in_reg"), type="character", default="NONE", help="input - path of list of regulator ids")
  p_in_sample = make_option(c("--p_in_sample"), type="character", default="NONE", help="path of list of samples ids")
  flag_global_shrinkage = make_option(c("--flag_global_shrinkage"), type="character", default="OFF", help="ON or OFF for netprophet1.0 global shrinkage")
  flag_local_shrinkage = make_option(c("--flag_local_shrinkage"), type="character", default="ON", help="ON or OFF for netprophet1.0 local shrinkage")
  greenfield_method = make_option(c("--greenfield_method"), type="character", default="OFF", help="ONE or TWO for greenfield method 1 or 2")
  fname_lasso = make_option(c("--fname_lasso"), type="character", default=NULL, help="output - path of generated lasso network")
  p_out_dir = make_option(c("--p_out_dir"), type="character", default=NULL, help="output - path of output directory for results")
  flag_parallel = make_option(c("--flag_parallel"), type="character", default="OFF", help="ON or OFF flag for faster parallel computation")
  seed = make_option(c("--seed"), type="integer", default=747, help="seed for reproducibility")
  nbr_cv_fold = make_option(c("--nbr_cv_fold"), type="integer", default=10, help="number of folds for cross validation, default 10")
  p_src_code = make_option(c("--p_src_code"), type="character", default=NULL, help="path for netpropohet source code")
  flag_microarray = make_option(c("--flag_microarray"), type="character", default="MICROARRAY", help="MICROARRAY or RNA-SEQ: for scaling normalizing the expression data")
  p_lasso_greenfield = make_option(c("--p_lasso_greenfield"), type="character", default=NULL, help="path of full lasso network without removing any TF. This parameter is mandatory for greenfiled method 2, otherwise the program doens't run properly")
  idx_reg_greenfield = make_option(c("--idx_reg_greenfield"), type="integer", default=NULL, help="")
  
  opt_parser = OptionParser(option_list=list(p_in_target, p_in_reg, p_in_sample, p_in_expr_target, p_in_expr_reg
                                             , flag_global_shrinkage, flag_local_shrinkage
                                             , greenfield_method, p_out_dir, fname_lasso, flag_parallel
                                             , seed, nbr_cv_fold, p_src_code, flag_microarray, p_lasso_greenfield
                                             , idx_reg_greenfield
                                             ))
  opt = parse_args(opt_parser)
  
  if (is.null(opt$p_in_target) || is.null(opt$p_in_reg) || is.null(opt$p_in_sample)
      || is.null(opt$p_in_expr_target) || is.null(opt$p_in_expr_reg)
      ){
    print_help(opt_parser)
    stop("all arguments p_in_expr_target, p_in_expr_reg are mandatory")
  }
  
  
  p_in_target=opt$p_in_target
  p_in_reg=opt$p_in_reg
  p_in_sample=opt$p_in_sample
  p_in_expr_target=opt$p_in_expr_target
  p_in_expr_reg=opt$p_in_expr_reg
  flag_global_shrinkage=opt$flag_global_shrinkage
  flag_local_shrinkage=opt$flag_local_shrinkage
  greenfield_method=opt$greenfield_method
  p_out_dir=opt$p_out_dir
  fname_lasso=opt$fname_lasso
  flag_parallel=opt$flag_parallel
  seed=opt$seed
  nbr_cv_fold = opt$nbr_cv_fold
  p_src_code = opt$p_src_code
  flag_microarray = opt$flag_microarray
  p_lasso_greenfield = opt$p_lasso_greenfield
  idx_reg_greenfield = opt$idx_reg_greenfield

  generate_lasso_net(p_in_target=p_in_target
                                 , p_in_reg=p_in_reg
                                 , p_in_sample=p_in_sample
                                 , p_in_expr_target=p_in_expr_target
                                 , p_in_expr_reg=p_in_expr_reg
                                 , flag_global_shrinkage=flag_global_shrinkage
                                 , flag_local_shrinkage=flag_local_shrinkage
                                 , greenfield_method=greenfield_method
                                 , p_out_dir=p_out_dir
                                 , fname_lasso=fname_lasso
                                 , flag_parallel=flag_parallel
                                 , seed=seed
                                 , nbr_cv_fold=nbr_cv_fold
                                 , p_src_code=p_src_code
                                 , flag_microarray=flag_microarray
                                 , p_lasso_greenfield=p_lasso_greenfield
                                 , idx_reg_greenfield=idx_reg_greenfield)
}


