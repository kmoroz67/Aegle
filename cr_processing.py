TCR_AB_LEN_RANGE = range(24, 61, 3)
TCR_G_LEN_RANGE = range(12, 61, 3)
TCR_D_LEN_RANGE = range(12, 91, 3)
BCR_LEN_RANGE = range(9, 91, 3)

def define_vdj_regions(clonotypes, raw_data):
    """
    Define VDJ regions and its borders
    :param clonotypes: dataframe with clonotypes to update
    :param raw_data: dataframe with raw clonotypes
    :return: updated dataframe with clonotypes
    """

    possible_offsets_list = list(map(lambda x: list(map(lambda y: y.split(':')[9], x.split(','))),
                                     raw_data['refPoints']))

    offsets = list(map(lambda x: (int(max(x)), np.argmax(x)), possible_offsets_list))

    for r in ['v', 'd', 'j']:
        clonotypes[r + '_region'] = list(map(lambda x: x.split('*')[0],
                                             raw_data['all' + str.upper(r) + 'HitsWithScore'].replace(
                                                 np.nan, '')))

        f = lambda x, offset_num: x.split(',')[offset_num].split('|')
        left_border = lambda x, offset: max(int(f(x, offset[1])[3]) - offset[0], 0)
        right_border = lambda x, offset: int(f(x, offset[1])[4]) - 1 - offset[0]

        clonotypes[[r + '_start_nt', r + '_end_nt', r + '_start_aa', r + '_end_aa']] = \
            list(map(lambda x, offset: [left_border(x, offset), right_border(x, offset),
                                        int(left_border(x, offset) / 3), int(right_border(x, offset) / 3)
                                        ] if pd.notnull(x) else ['', '', '', ''],
                     raw_data['all' + str.upper(r) + 'Alignments'], offsets))

    clonotypes.c_region = list(map(lambda x: x.split('*')[0],
                                   raw_data.allCHitsWithScore.replace(np.nan, '')))

    return clonotypes

def make_clonotypes_df(raw_data, analysis_type='BCR'):
    """
    Parses raw MiXCR file and make clonotype table for subsequent analysis.
    :param raw_data: dataframe with mixcr output
    :param analysis_type: ('TCR', 'BCR')
    :return: dataframe with clonotypes
    """

    if analysis_type not in ['TCR', 'BCR']:
        raise Exception("The analysis_type parameter must be either 'TCR' or 'BCR'")

    raw_data = raw_data.sort_values('cloneCount', ascending=False).reset_index(drop=True)

    clonotypes = pd.DataFrame('', index=raw_data.index,
                              columns=['clonotype_id', 'fraction', 'clonotype_count', 'cdr3aa',
                                       'cdr3nt', 'cdr3nt_length', 'chain',
                                       'v_region', 'd_region', 'j_region', 'c_region',
                                       'v_start_nt', 'v_end_nt', 'v_start_aa', 'v_end_aa',
                                       'd_start_nt', 'd_end_nt', 'd_start_aa', 'd_end_aa',
                                       'j_start_nt', 'j_end_nt', 'j_start_aa', 'j_end_aa'])

    clonotypes['clonotype_count'] = raw_data['cloneCount']
    clonotypes['cdr3aa'] = raw_data['aaSeqImputedCDR3']
    clonotypes['cdr3nt'] = raw_data['nSeqImputedCDR3']
    clonotypes['cdr1aa'] = raw_data['aaSeqImputedCDR1']
    clonotypes['cdr1nt'] = raw_data['nSeqImputedCDR1']
    clonotypes['cdr2aa'] = raw_data['aaSeqImputedCDR2']
    clonotypes['cdr2nt'] = raw_data['nSeqImputedCDR2']
    
    clonotypes['fr1aa'] = raw_data['aaSeqImputedFR1']
    clonotypes['fr1nt'] = raw_data['nSeqImputedFR1']
    clonotypes['fr2aa'] = raw_data['aaSeqImputedFR2']
    clonotypes['fr2nt'] = raw_data['nSeqImputedFR2']
    clonotypes['fr3aa'] = raw_data['aaSeqImputedFR3']
    clonotypes['fr3nt'] = raw_data['nSeqImputedFR3']
    clonotypes['fr4aa'] = raw_data['aaSeqImputedFR4']
    clonotypes['fr4nt'] = raw_data['nSeqImputedFR4']
    clonotypes['cdr3nt_length'] = list(map(len, raw_data['nSeqImputedCDR3']))

    clonotypes = define_vdj_regions(clonotypes, raw_data)

    if analysis_type == 'BCR':
        clonotypes['chain'] = clonotypes.v_region.apply(lambda x: {
            'H': 'heavy', 'K': 'kappa', 'L': 'lambda'}.get(x[2], np.nan))

        heavy_clonotypes = clonotypes[clonotypes['chain'] == 'heavy'].index
        clonotypes.loc[heavy_clonotypes, 'isotype'] = \
            clonotypes.loc[heavy_clonotypes].c_region.apply(
                lambda x: 'Ig' + x[3] if x is not '' else np.nan)

    else:
        clonotypes['chain'] = clonotypes.v_region.apply(lambda x: {
            'A': 'alpha', 'B': 'beta', 'G': 'gamma', 'D': 'delta'}.get(x[2], np.nan))

    clonotypes.dropna(subset=['chain'], inplace=True)
    clonotypes = merge_identical(clonotypes)
    clonotypes = check_regions_borders(clonotypes, max_seq_len=100)

    return clonotypes
  
def merge_identical(clonotypes):
    """
    Clonotypes that have identical CDR3 nucleotide sequences and V- and J-regions
    are merged into one clonotype.
    C-region difference is ignored, with the exception of BCR heavy chain.
    """
    heavy_chain = (clonotypes.chain == 'heavy')

    df = pd.concat([merge_identical_by_chain(clonotypes[heavy_chain], True),
                    merge_identical_by_chain(clonotypes[~heavy_chain], False)])

    df = df.sort_values('clonotype_count', ascending=False).reset_index(drop=True)

    return df
  
def merge_identical_by_chain(df, heavy_chain_flag):
    """
    Merges identical clonotypes in one chain. In BCR heavy chain clonotypes are
    additionally grouped by isotype
    :param df: dataframe with clonotypes
    :param heavy_chain_flag: bool, whether it's heavy chain
    :return: updated dataframe with clonotypes
    """
    if df.shape[0] == 0:
        return None

    columns = ['chain', 'cdr3nt', 'v_region', 'j_region']
    if heavy_chain_flag:
        columns = columns + ['isotype']

    clCount = df.fillna('nan').groupby(columns, as_index=False).clonotype_count.transform('sum')

    df.drop_duplicates(columns, inplace=True)
    df.update(clCount)

    return df
  
def check_regions_borders(clonotypes, max_seq_len=100):
    """
    Check borders of VDJ genes for overlapping coordinates and so on.
    :param clonotypes: dataframe with clonotypes, including borders
    :param max_seq_len: maximum lenght of cdr3nt sequence (necessary for Hippocrates)
    :param clonotypes: dataframe with updated borders
    """

    clonotypes.j_end_nt = list(map(lambda j, l: int(min(j, l)),
                                   clonotypes.j_end_nt, clonotypes.cdr3nt_length - 1))
    clonotypes.j_end_aa = list(map(lambda j, l: int(min(j, l)),
                                   clonotypes.j_end_aa, clonotypes.cdr3nt_length / 3 - 1))

    clonotypes.loc[clonotypes[(clonotypes.j_start_nt - clonotypes.cdr3nt_length) >= 0].index,
                   ['j_start_nt', 'j_end_nt', 'j_start_aa', 'j_end_aa']] = ''

    clonotypes.loc[clonotypes[clonotypes.v_end_nt <= 0].index,
                   ['v_start_nt', 'v_end_nt', 'v_start_aa', 'v_end_aa']] = ''
    if max_seq_len:
        clonotypes.loc[clonotypes[clonotypes.cdr3nt_length >= max_seq_len].index,
                       'v_start_nt':'j_end_aa'] = ''

    # Check for overlapping End and Start in adjacent regions

    for i in clonotypes.index:
        for r1, r2 in zip(['v', 'd', 'v'], ['d', 'j', 'j']):
            for t in ['nt', 'aa']:
                el1 = clonotypes.loc[i, r1 + '_end_' + t]
                el2 = clonotypes.loc[i, r2 + '_start_' + t]
                if (el1 != '') and (el2 != ''):
                    if el1 >= el2:
                        clonotypes.loc[i, r2 + '_start_' + t] = clonotypes.loc[i, r1 + '_end_' + t] + 1
                        if clonotypes.loc[i, r2 + '_start_' + t] >= clonotypes.loc[i, r2 + '_end_' + t]:
                            clonotypes.loc[i, r2 + '_start_' + t] -= 1
                            clonotypes.loc[i, r1 + '_end_' + t] -= 1

    # Check for overlappling End and Start in the same regions

    for i in clonotypes.index:
        for r in ['v', 'd', 'j']:
            for t in ['nt', 'aa']:
                if (clonotypes.loc[i, r + '_end_' + t] != '') and (clonotypes.loc[i, r + '_end_' + t] != ''):
                    if clonotypes.loc[i, r + '_start_' + t] > clonotypes.loc[i, r + '_end_' + t]:
                        clonotypes.loc[i, r + '_start_' + t] = clonotypes.loc[i, r + '_end_' + t]

    return clonotypes
  
def calc_clonotype_fraction(chain_df):
    """
    Function calculates fraction by chain
    """
    df = chain_df.copy()
    df['fraction'] = df.clonotype_count / sum(df.clonotype_count)

    return df
  
def process_chain_clonotypes(chain_clonotypes):
    """
    Function processes dataframe with clonotypes
    :param chain_clonotypes: data with clonotypes
    :param analysis_type: TCR or BCR
    :param diagnosis: str, diagnosis
    :param t_or_b_cell_cancers: list of T/B cell cancers
    :return: dataframe with processed clonotypes
    """
    chain_clonotypes = calc_clonotype_fraction(chain_clonotypes)
    chain_clonotypes['clonotype_id'] = chain_clonotypes.index + 1

    return chain_clonotypes
  
def make_cr_dataframe(raw_tcr_file_path, output_processed_file_path,
                       analysis_type='TCR'):
    """
    Function takes as input raw BCR/TCR files from mixcr and create processed dataframe file.
    :param raw_tcr_file_path: path to file with raw TCR/BCR data
    :param output_processed_file_path: path to "XCR_processed.txt"
    :param analysis_type: TCR or BCR
    :return: None
    """

    raw_clonotypes = pd.read_csv(raw_tcr_file_path, sep='\t',engine='python',error_bad_lines=False)

    concat_clonotypes = pd.DataFrame()

    if raw_clonotypes.shape[0] == 0:
        pass
    else:
        # Make clones table
        clonotypes = make_clonotypes_df(raw_clonotypes, analysis_type)

        # Work with each chain separately
        for chain_name, chain_clonotypes in clonotypes.groupby(by='chain'):
            chain_clonotypes = process_chain_clonotypes(chain_clonotypes)

            concat_clonotypes = pd.concat([concat_clonotypes, chain_clonotypes])  # save updated chain dataframe

    # save concat clonotypes for later use in notebook analysis
    concat_clonotypes = concat_clonotypes.reset_index(drop=True)
    concat_clonotypes.to_csv(output_processed_file_path, sep='\t', index=False, header=True)
