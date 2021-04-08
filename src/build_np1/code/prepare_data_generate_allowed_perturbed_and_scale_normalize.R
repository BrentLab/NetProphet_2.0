# ================================================================================================== #
# |                                      **** PREPARE DATA ****                                    | #
# | Implementation of functions that prepare data for running LASSO:                               | #
# | 1. generate_allowed_perturbed_matrices: allowed matrice is a network matrix (reg x target):    | #
# |    0 when the reg is the same as the target, 1 anything else.                                  | #
# | 2. scale_normalize_expr_matrix: scale and normalize expression matrix depending on RNA-Seq or  | #
# |    microarray data                                                                             | #
# ================================================================================================== #

generate_allowed_perturbed_matrices = function (l_in_target
                                                , l_in_reg
                                                , l_in_sample
                                                , p_out_dir
                                                , p_src_code=""
                                               ){
  "**
    *@param: l_in_target: list of target genes
    *@param: l_in_reg: list of regulators
    *@param: l_in_sample: list of samples ids for expression matrix
    *@param: p_out_dir: if not NULL, path of output directore where allowed and perturbed matrix can be written
    *@param: p_src_code: path of the source directory of the code
  *"

  source(paste(p_src_code, "src/build_np1/code/helper.R", sep=""))

  # generated the allowed matrix
  df_allowed = data.frame(matrix(1, length(l_in_reg), length(l_in_target)))
  rownames(df_allowed) = l_in_reg
  colnames(df_allowed) = l_in_target
  for (reg in l_in_reg){
    if (!is.null(df_allowed[reg, reg]))
    {
      df_allowed[reg, reg] = 0
    }
  }
      
  # generate the perturbed matrix
  df_perturbed = data.frame(matrix(0, length(l_in_target), length(l_in_sample)))
  rownames(df_perturbed) = l_in_target
  colnames(df_perturbed) = l_in_sample
  for (target in l_in_target){
    if (!is.null(df_perturbed[target, target]))
    {
      df_perturbed[target, target] = 1
     }
  }
  
  # write allowed and perturbed matrices, if output directory doesn't exit, create it    
  if (!is.null(p_out_dir)){
    ifelse(!dir.exists(file.path(p_out_dir, 'tmp')), dir.create(file.path(p_out_dir, 'tmp'), showWarnings = FALSE), FALSE)
    # write allowed matrix
    write.table(file=file(paste(p_out_dir, 'tmp/', 'allowed.tsv', sep=''))
                , x=df_allowed
                , row.names=FALSE
                , col.names=FALSE
                , sep="\t")
   # write perturbed matrix
   write.table(file=file(paste(p_out_dir, 'tmp/','perturbed.tsv', sep=''))
               , x=df_perturbed
               , row.names=FALSE
               , col.names=FALSE
               , sep="\t")
  }
  
  # return allowed and perturbed matrices
  data = list()
  data[[1]] = df_allowed
  data[[2]] = df_perturbed
  
  data
}

scale_normalize_expr_matrix = function(df_expr, flag_microarray){
  "**
    *@param: df_expr: expression matrix targets x samples
    *@param: flag_microarray: MICROARRAY for microarray data, RNA-SEQ for RNA-Seq data
  *"
  #RNA-Seq Data
  if (flag_microarray == "RNA-SEQ"){
      cat("rna-seq processing..")
      df_expr <- log(df_expr+1)/log(2)
  }
    
  df_expr = df_expr - apply(df_expr, 1, mean)
  sd_expr = apply(df_expr,1,sd)
  sd_floor_expr = mean(sd_expr) + sd(sd_expr)
  norm_expr = apply(rbind(rep(sd_floor_expr,times=length(sd_expr)),sd_expr),2,max) / sqrt(dim(df_expr)[2])
  df_expr = df_expr / (norm_expr* sqrt(dim(df_expr)[2]-1))
  
  df_expr
}