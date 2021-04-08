def subsample_data( # Input
                    p_in_reg
                  , p_in_target
                  , p_in_condition
                  , p_in_expr_target
                  , p_in_expr_reg
                  , p_in_de
                  , nbr_sub_reg
                  , nbr_sub_target
                  , nbr_sub_condition
                    # Output
                  , p_out_sub_reg
                  , p_out_sub_target
                  , p_out_sub_condition
                  , p_out_sub_expr_target
                  , p_out_sub_expr_reg
                  , p_out_sub_de
                  ):
    # load libraries
    from pandas import read_csv, DataFrame, Series
    from random import sample
    
    # read input data
    l_in_reg = list(read_csv(p_in_reg, header=None)[0])
    l_in_target = list(read_csv(p_in_target, header=None)[0])
    l_in_condition = list(read_csv(p_in_condition, header=None)[0])
    df_in_expr_target = read_csv(p_in_expr_target, header=0, index_col=0, sep='\t')
    df_in_expr_reg = read_csv(p_in_expr_reg, header=0, index_col=0, sep='\t')
    df_in_de = read_csv(p_in_de, header=0, index_col=0, sep='\t')
    
    # sample data
    l_out_sub_reg = sample(l_in_reg, nbr_sub_reg)
    l_out_sub_target = sample(l_in_target, nbr_sub_target)
    l_out_sub_condition = sample(l_in_condition, nbr_sub_condition)
    df_out_sub_expr_target = df_in_expr_target.loc[df_in_expr_target.index.isin(l_out_sub_target)
                                                   , df_in_expr_target.columns.isin(l_out_sub_condition)]
    df_out_sub_expr_target = df_out_sub_expr_target.reindex(l_out_sub_target, axis='index')
    df_out_sub_expr_target = df_out_sub_expr_target.reindex(l_out_sub_condition, axis='columns')
    
    df_out_sub_expr_reg = df_in_expr_reg.loc[df_in_expr_reg.index.isin(l_out_sub_reg)
                                            , df_in_expr_reg.columns.isin(l_out_sub_condition)]
    df_out_sub_expr_reg = df_out_sub_expr_reg.reindex(l_out_sub_reg, axis='index')
    df_out_sub_expr_reg = df_out_sub_expr_reg.reindex(l_out_sub_condition, axis='columns')
    
    df_out_sub_de = df_in_de.loc[df_in_de.index.isin(l_out_sub_reg), df_in_de.columns.isin(l_out_sub_target)]
    df_out_sub_de = df_out_sub_de.reindex(l_out_sub_reg, axis='index')
    df_out_sub_de = df_out_sub_de.reindex(l_out_sub_target, axis='columns')

    # write sub-sampled data
    Series(l_out_sub_reg, name='reg').to_csv(p_out_sub_reg, header=False, index=False)
    Series(l_out_sub_target, name='target').to_csv(p_out_sub_target, header=False, index=False)
    Series(l_out_sub_condition, name='condition').to_csv(p_out_sub_condition, header=False, index=False)
    df_out_sub_expr_target.to_csv(p_out_sub_expr_target, header=True, index=True, sep='\t')
    df_out_sub_expr_reg.to_csv(p_out_sub_expr_reg, header=True, index=True, sep='\t')
    df_out_sub_de.to_csv(p_out_sub_de, header=True, index=True, sep='\t')
        
        
def main():
        p_wd = '/scratch/mblab/dabid/netprophet/net_in/'
        # zev
        nbr_sub_reg = 50
        nbr_sub_target = 500
        nbr_sub_condition = 100
        subsample_data( # Input
                    p_in_reg = p_wd + 'zev_in_reg_tf_167'
                  , p_in_target = p_wd + 'zev_in_target_6175'
                  , p_in_condition = p_wd + 'zev_in_condition_591'
                  , p_in_expr_target = p_wd + 'zev_expr_6175_591_indexed'
                  , p_in_expr_reg = p_wd + 'zev_expr_reg_167_591_indexed'
                  , p_in_de = p_wd + 'zev_de_shrunken_167_6175_indexed'
                  , nbr_sub_reg = nbr_sub_reg
                  , nbr_sub_target = nbr_sub_target
                  , nbr_sub_condition = nbr_sub_condition
                    # Output
                  , p_out_sub_reg = p_wd + 'toy_example/zev_in_reg_tf_' + str(nbr_sub_reg)
                  , p_out_sub_target = p_wd + 'toy_example/zev_in_target_' + str(nbr_sub_target)
                  , p_out_sub_condition = p_wd + 'toy_example/zev_in_condition_' + str(nbr_sub_condition)
                  , p_out_sub_expr_target = p_wd + 'toy_example/zev_expr_' + str(nbr_sub_target) + '_' + str(nbr_sub_condition) + '_indexed'
                  , p_out_sub_expr_reg = p_wd + 'toy_example/zev_expr_reg_' + str(nbr_sub_reg) + '_' + str(nbr_sub_condition) + '_indexed'
                  , p_out_sub_de = p_wd + 'toy_example/zev_de_shrunken_' + str(nbr_sub_reg) + '_' + str(nbr_sub_target) + '_indexed'
                  )
        # kem
        
if __name__ == '__main__':
    main()
