def join_duplicates(cr_df, cdr_type):
    """
    Join duplicates by cdr3 aa/nt (use after V(D)J and chain selection )

    :param cr_df: dataframe with processed TCR/BCR (.txt)
    :param cdr_type: cdr3aa - for amino acids, cdr3nt for nucleotides
    :return: dataframe of joined cdr3
    """
    
    recalculated_df = cr_df[['clonotype_count',cdr_type]].groupby(by=cdr_type, as_index=False).sum()
    recalculated_df = recalculated_df.sort_values(by=cdr_type,ascending=False)

    cr_df = cr_df.groupby(cdr_type, as_index=False).first().sort_values(
        by='clonotype_count', ascending=False)
    cr_df = cr_df.sort_values(by=cdr_type,ascending=False)
    cr_df['clonotype_count'] = recalculated_df.clonotype_count
    cr_df = cr_df.sort_values(by='clonotype_count',ascending=False)
    
    return cr_df
  
def recalculate_chain_fraction(cr_df):
    """
    Recalculate fraction of each cdr3 of one chain

    :param cr_df: dataframe with processed TCR/BCR (.txt)
    :return: dataframe with recalculated fractions
    """
    
    res_df = pd.DataFrame()
    for chain in cr_df['chain'].unique():
        chain_df = cr_df[cr_df['chain']==chain]
        chain_count_sum = chain_df['clonotype_count'].sum()
        chain_df['fraction'] = chain_df['clonotype_count'].apply(lambda x: x / chain_count_sum)
        res_df = pd.concat([res_df,chain_df])
    return res_df
  
def qc_cr(label, prc_cr_df, drop_singletons=False):
    """
    Make a qc of mixcr data after its processing

    :param label: type of cdr3. Available : tcr-ab,tcr-gd, bcr.
    :param prc_cr_df: dataframe with processed TCR/BCR (txt)
    :param drop_singletons: True if you need to drop singletons
    :return: dataframe - dataframe of found cdr3
    """
    
    prc_cr_df = prc_cr_df[(prc_cr_df['v_start_nt'].notna()) & (prc_cr_df['v_end_nt'].notna()) \
                          & (prc_cr_df['j_start_nt'].notna()) & (prc_cr_df['j_end_nt'].notna())]
    prc_cr_df = prc_cr_df[(prc_cr_df['v_start_nt'] >= 0) & (prc_cr_df['v_end_nt'] > 0) \
                          & (prc_cr_df['j_start_nt'] > 0) & (prc_cr_df['j_end_nt'] > 0)]
    
    prc_cr_df = prc_cr_df[prc_cr_df['cdr3aa'].str.startswith('C')==True]
    
    if label == 'tcr':
        prc_cr_df_ab = prc_cr_df[prc_cr_df['chain'].isin(['alpha','beta'])]
        prc_cr_df_ab = prc_cr_df_ab[prc_cr_df_ab['cdr3nt_length'].isin(TCR_AB_LEN_RANGE)]
        
        prc_cr_df_gd = prc_cr_df[prc_cr_df['chain'].isin(['gamma','delta'])]
        
        res_df = pd.DataFrame()
        prc_cr_df_g = prc_cr_df_gd[prc_cr_df_gd['chain']=='gamma']
        prc_cr_df_g = prc_cr_df_g[prc_cr_df_g['cdr3nt_length'].isin(TCR_G_LEN_RANGE)]
        res_df = pd.concat([res_df,prc_cr_df_g])
        
        prc_cr_df_d = prc_cr_df_gd[prc_cr_df_gd['chain']=='delta']
        prc_cr_df_d = prc_cr_df_d[prc_cr_df_d['cdr3nt_length'].isin(TCR_D_LEN_RANGE)]
        res_df = pd.concat([res_df,prc_cr_df_d])
        
        prc_cr_df = pd.concat([prc_cr_df_ab,res_df])
        prc_cr_df = prc_cr_df[(prc_cr_df['cdr3aa'].str.endswith('F')==True)]
        
    elif label == 'bcr':
        prc_cr_df = prc_cr_df[prc_cr_df['cdr3nt_length'].isin(BCR_LEN_RANGE)]
        
        res_df = pd.DataFrame()
        prc_cr_df_h = prc_cr_df[prc_cr_df['chain']=='heavy']
        prc_cr_df_h = prc_cr_df_h[(prc_cr_df_h['cdr3aa'].str.endswith('W')==True)]
        res_df = pd.concat([res_df,prc_cr_df_h])
        
        prc_cr_df_l = prc_cr_df[prc_cr_df['chain'].isin(['lambda','kappa'])]
        prc_cr_df_l = prc_cr_df_l[(prc_cr_df_l['cdr3aa'].str.endswith('F')==True)]
        res_df = pd.concat([res_df,prc_cr_df_l])
        
        prc_cr_df = res_df
    else: 
        print('Wrong type of cdr3')
    
    prc_cr_df['vj'] = prc_cr_df[['v_region', 'j_region']].apply(lambda x: x[0] + ' ' + x[1], axis=1)
    res_df = pd.DataFrame()
    for chain in prc_cr_df.chain.unique():
        chain_df = prc_cr_df[prc_cr_df['chain']==chain]
        chain_df = join_duplicates(chain_df,'cdr3nt')
        res_vj = pd.DataFrame() 
        
        #print('start vj in {}'.format(chain))
        for vj in chain_df['vj'].unique():
            vj_df = chain_df[chain_df['vj']==vj]
            vj_df = join_duplicates(vj_df,'cdr3aa')
            res_vj = pd.concat([res_vj,vj_df])
        
        res_df = pd.concat([res_df,res_vj])
        #print('end vj in {}'.format(chain))
        
    prc_cr_df = res_df
    prc_cr_df = recalculate_chain_fraction(prc_cr_df)
    
    total_count_sum = prc_cr_df['clonotype_count'].sum()
    prc_cr_df['frequency'] = prc_cr_df['clonotype_count'].apply(lambda x: x / total_count_sum)
    
    if drop_singletons:
        prc_cr_df = prc_cr_df[prc_cr_df['clonotype_count'] > 1]
    
    return prc_cr_df
