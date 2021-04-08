read_csv_indexed = function(p_df
                            , p_rownames
                            , p_colnames
                            , sep="\t"){
  df = read.csv(p_df, header=FALSE, sep=sep)
  rownames = read.csv(p_rownames, header=FALSE)[1]
  colnames = read.csv(p_colnames, header=FALSE)[1]
  df$row.names = rownames
  df$col.names = colnames
  df
}

calculate_mse = function(l_in_target
                         , l_in_sample
                         , df_expr_target
                         , df_expr_reg
                         , df_lasso_net){
  "**
    *@param: l_in_target: list of target genes
    *@param: l_in_sample: list of samples ids for expression matrix
    *@param: df_expr_target: expression matrix of target genes (targets x samples)
    *@param: df_expr_reg: expression matrix of target genes (reg x samples)
    *@param: df_lasso_net: matrix or dataframe of lasso network (reg x targets)
  *"
  mse =  list()
  for (idx_target in seq(1, length(l_in_target), 1)){
    mse[idx_target] = 1/length(l_in_sample) * sum((df_expr_target[idx_target, ] - (t(df_expr_reg) %*% df_lasso_net[, idx_target]))**2)
  }
  mse
}



