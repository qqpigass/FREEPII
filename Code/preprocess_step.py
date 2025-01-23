#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import re
import pyreadr
import pickle
import argparse
from itertools import combinations
from collections import OrderedDict, Counter
from scipy import sparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('-root',        '--root',             default='/FREEPII',   help='root',                          type=str)
parser.add_argument('-sp',          '--species',          default='Human',      help='species',                       type=str)
parser.add_argument('-e_name',      '--exp_name',         default='PXD002892',  help='exp_name',                      type=str)
parser.add_argument('-e_cond',      '--exp_cond',         default='SEC2-heavy', help='exp_condition',                 type=str)
parser.add_argument('-split_ratio', '--train_test_ratio', default='7',          help='train_test_ratio',              type=int)
parser.add_argument('-cv_fold',     '--cv_fold',          default='5',          help='cv_fold_number',                type=int)
parser.add_argument('-pn_ratio',    '--pn_ratio',         default='1',          help='pos_to_neg_ratio',              type=int)
parser.add_argument('-out_path',    '--output_path',      default='/Output',    help='output_path',                   type=str)
parser.add_argument('-pretrain_w',  '--pretrain_w',       default='null',       help='path_to_pretrain_model_weight', type=str)
args = parser.parse_args()
# args = parser.parse_args([])

if '\\' in os.getcwd():
    args.root = ['/'.join(os.getcwd().split('\\')[:(i+1)]) for i,v in enumerate(os.getcwd().split('\\')) 
                 if v=='FREEPII'][0]

args.root = ['/'.join(os.getcwd().split('/')[:(i+1)]) for i,v in enumerate(os.getcwd().split('/')) if v=='FREEPII'][0]

# if args.output_path=='/Output':
#     args.output_path = '/'.join([args.root + args.output_path, args.exp_name, args.exp_cond])

# if not os.path.exists(args.output_path):
#     os.makedirs(args.output_path)

def process_complex(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', cur_root='/FREEPII'):
    complex_raw_path  = cur_root + '/Protein complex-data/allComplexes_Corum.txt'
    complex_path      = cur_root + '/Protein complex-data/Complexes_gene_' + cur_species + '_filter.csv'
    uniprot_gene_path = cur_root + '/UniProt-data/uniprot_gene_' + cur_species + '.rds'
    
    ## Check if preprocessed data exists
    process_indicator = os.path.exists(complex_path)==False
    if process_indicator:
        print('Start preprocess protein complexes')
        
        #### Read data
        complex_df_raw = pd.read_csv(complex_raw_path, delimiter = "\t")
        complex_df_raw = complex_df_raw.loc[complex_df_raw['Organism']==cur_species].drop_duplicates()
        complex_df_raw = complex_df_raw[['ComplexName', 'subunits(UniProt IDs)', 'subunits(Gene name)', 
                                         'subunits(Entrez IDs)']].drop_duplicates()
        complex_df_raw.columns = ['ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez']
        print('Number of unique complexes in source: ', len(set(complex_df_raw['ComplexName'])))
        
        uniprot_df = list(pyreadr.read_r(uniprot_gene_path).values())[0]
        
        
        #### Release genes in complex df
        ## Match split in each row
        print('Release genes in complexes')
        
        complex_df_list = []
        for ComplexName in list(set(complex_df_raw['ComplexName'])):
            # print(ComplexName)
            uniprot_id = complex_df_raw.loc[complex_df_raw['ComplexName']==ComplexName]['UniProt_entry'].values[0].split(';')
            g_name = complex_df_raw.loc[complex_df_raw['ComplexName']==ComplexName]['Gene_name'].values[0].split(';')
            g_entrez = complex_df_raw.loc[complex_df_raw['ComplexName']==ComplexName]['Gene_entrez'].values[0].split(';')
            
            uniprot_id = [i for i in uniprot_id if i != '']
            g_name = [i for i in g_name if i != '']
            g_entrez = [i for i in g_entrez if i != '']
            
            complex_df_list += [[ComplexName, uniprot_id[i], g_name[i], g_entrez[i]] for i in range(len(uniprot_id))]
        
        complex_df = pd.DataFrame(complex_df_list, columns =['ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'])
        
        ## Filter None in UniProt_entry and Gene_name
        print('Filter None in complexes')
        complex_df = complex_df.loc[(complex_df['UniProt_entry']!='None')&(complex_df['Gene_name']!='None')].drop_duplicates()
        
        
        #### Match UniProt ID and gene name to complex
        print('Match UniProt ID and Gene Name')
        complex_df = pd.merge(complex_df, uniprot_df[['UniProt_entry', 'Gene_name', 'Gene_name_primary']], 
                              how='left', on='Gene_name').drop_duplicates()
        
        if complex_df.isnull().values.any():
            complex_df = complex_df.fillna('None')
        
        a = complex_df.loc[complex_df['UniProt_entry_x']==complex_df['UniProt_entry_y']]
        b = complex_df.loc[complex_df['UniProt_entry_x']!=complex_df['UniProt_entry_y']]
        b = pd.melt(b, id_vars=['ComplexName', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'], 
                    value_vars=['UniProt_entry_x', 'UniProt_entry_y'])
        b.columns = ['ComplexName', 'Gene_name', 'Gene_name_primary', 'Gene_entrez', 'Group', 'UniProt_entry_x']
        
        complex_df = pd.concat([a[['ComplexName', 'UniProt_entry_x', 'Gene_name', 'Gene_name_primary', 'Gene_entrez']], 
                                b[['ComplexName', 'UniProt_entry_x', 'Gene_name', 'Gene_name_primary', 'Gene_entrez']]]).drop_duplicates()
        complex_df.columns = ['ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez']
        complex_df = complex_df.loc[complex_df['UniProt_entry']!='None'].drop_duplicates()
        
        
        #### Filter complex by size (Gene_name)
        print('Filter complexes by size')
        
        freq_df = [[i, 
                    len([j for j in set(complex_df.loc[complex_df['ComplexName']==i]['UniProt_entry'].values) if j !='None']),
                    len([j for j in set(complex_df.loc[complex_df['ComplexName']==i]['Gene_name'].values) if j !='None']),
                    len([j for j in set(complex_df.loc[complex_df['ComplexName']==i]['Gene_name_primary'].values) if j !='None'])
                   ] for i in set(complex_df['ComplexName'])]
        freq_df = pd.DataFrame(freq_df, columns =['ComplexName', 'UniProt_entry_N', 'Gene_name_N', 'Gene_name_primary_N'])
        
        complex_df = pd.merge(complex_df, freq_df, how='left', on='ComplexName').drop_duplicates()
        complex_df_sub = complex_df.loc[complex_df['Gene_name_N'] > 2]
        
        
        #### Save output and show
        complex_df_sub.to_csv(complex_path, encoding='utf-8', index=False, header=True)
        print('Number of unique complexes after process: ', len(set(complex_df_sub['ComplexName'])))
        print('Done preprocess protein complexes')
        print(' ')
    else:
        print('Preprocessed protein complexes already exist')
        print(' ')


def process_epf(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', cur_root='/FREEPII', quant_type='iBAQ'):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    if not os.path.exists(input_path):
        os.makedirs(input_path)
    
    epf_raw_path = '/'.join([cur_root, 'EPF-data', cur_exp_name, cur_exp_cond, quant_type + '.tsv'])
    epf_path     = '/'.join([cur_root, 'EPF-data', cur_exp_name, cur_exp_cond, cur_exp_cond + '_epf.csv'])
    dict_path    = '/'.join([input_path, 'name_idx_dict_' + cur_exp_cond + '.pickle'])
    
    print('Start preprocess elution profile of ' + cur_exp_name + ' ' + cur_exp_cond)
    
    #### Read EPF, normalize, match protein group to gene, and return matrix with row name = gene name
    mat = pd.read_csv(epf_raw_path, sep = '\t')
    meta = pd.read_csv(re.sub(quant_type, 'metadata', epf_raw_path), sep = '\t')
    
    ## Replace missing values in EPF with zero
    mat.replace([np.inf, -np.inf], 0, inplace=True)
    mat.fillna(0, inplace=True)
    
    ## Remove row with all zero
    meta = meta.loc[(mat!=0).any(axis=1).values, :]
    mat = mat.loc[(mat!=0).any(axis=1)]
    
    ## Row-wise normalize
    mat = mat.div(mat.sum(axis=1), axis=0)
    mat.replace([np.inf, -np.inf], 0, inplace=True)
    mat.fillna(0, inplace=True)
    
    ## keep one row per gene
    gene_map = meta[['Majority protein IDs', 'Gene names']]
    gene_map['Protein_IDs'] = gene_map['Majority protein IDs']
    gene_map['Gene_name'] = gene_map['Gene names'].str.split(';')
    gene_map = gene_map.explode('Gene_name')
    gene_map = gene_map[['Protein_IDs', 'Gene_name']].dropna()
    
    Gnames = sorted(set(gene_map['Gene_name']), key=lambda x: x.lower())
    n_fractions = dict((mat.notna() & mat.notnull() & ~mat.isin([np.inf, -np.inf]) & mat!=0).sum(axis=1))
    
    gene_mat = pd.DataFrame(np.zeros([len(Gnames), mat.shape[1]]))
    gene_mat.columns = mat.columns
    gene_mat.index = Gnames
    
    out_map = []
    for gname in Gnames:
        protein_id = gene_map.loc[gene_map['Gene_name']==gname]['Protein_IDs'].values.tolist()
        
        # pick the best protein for this replicate
        n_fraction = [n_fractions[i] for i in protein_id]
        best = protein_id[n_fraction.index(max(n_fraction))]
        gene_mat.loc[gname, ] = mat.loc[best, ]
        out_map += [[best, gname]]
    out_map = pd.DataFrame(out_map, columns =['Protein_IDs', 'Gene_name'])
    
    gene_mat.columns = ['F'+ str(i) for i in range(1, gene_mat.shape[1]+1)]
    
    
    #### Save epf and map
    gene_mat.to_csv(epf_path, index=True, header=True, index_label='Gene_name')
    out_map.to_csv(re.sub('epf.csv', 'map.csv', epf_path), index=False, header=True)
    
    
    #### Save idx dict
    epf_name_idx_dict = dict(zip(gene_mat.index.values, np.arange(gene_mat.shape[0])))
    
    with open(dict_path, 'wb') as f:
        pickle.dump(epf_name_idx_dict, f, pickle.HIGHEST_PROTOCOL)
    
    print(', '.join(['Finish preprocess EPF', cur_exp_name, cur_exp_cond, 
                     'dim: ' + str(gene_mat.shape[0]) + ' ' + str(gene_mat.shape[1])]))
    print(' ')


def process_seq(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', cur_root='/FREEPII'):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    if not os.path.exists(input_path):
        os.makedirs(input_path)
    
    epf_path          = '/'.join([cur_root, 'EPF-data', cur_exp_name, cur_exp_cond, cur_exp_cond + '_epf.csv'])
    uniprot_gene_path = cur_root + '/UniProt-data/uniprot_gene_' + cur_species + '.rds'
    seq_raw_path      = cur_root + '/FCGR-data/uniprot_seq_' + cur_species + '_FCGR_16x.rds'   
    seq_path          = '/'.join([input_path, 'seq_FCGR_16x_' + cur_exp_cond + '.npz'])
    dict_path         = '/'.join([input_path, 'name_idx_dict_' + cur_exp_cond + '.pickle'])
    
    print('Start preprocess protein sequences of ' + cur_exp_name + ' ' + cur_exp_cond)
    
    #### Read data
    epf_df = pd.read_csv(epf_path).set_index('Gene_name', drop=True)
    seq_df = list(pyreadr.read_r(seq_raw_path).values())[0]
    uniprot_df = list(pyreadr.read_r(uniprot_gene_path).values())[0]
    
    with open(dict_path, 'rb') as f:
        epf_name_idx_dict = pickle.load(f)
    
    match_order = sum([list(epf_df.index)[i]==list(epf_name_idx_dict.keys())[i] for i in range(epf_df.shape[0])])
    if match_order!=epf_df.shape[0]:
        print('re-create name-idx dictionary')
        epf_name_idx_dict = dict(zip(epf_df.index.values, np.arange(epf_df.shape[0])))
        
        with open(dict_path, 'wb') as f:
            pickle.dump(epf_name_idx_dict, f, pickle.HIGHEST_PROTOCOL)
    
    
    #### Generate input of portein sequences
    match_seq_num = []
    for idx in range(uniprot_df.select_dtypes(include='object').shape[1]):
        match_seq_num.append(len( set(uniprot_df.iloc[:, idx].values)&set(seq_df.columns.values) ))
    best_match_seq_name = uniprot_df.select_dtypes(include='object').columns.values[np.argmax(match_seq_num)]
    print('best match seq name: ', best_match_seq_name)

    match_epf_num = []
    for idx in range(uniprot_df.select_dtypes(include='object').shape[1]):
        match_epf_num.append(len( set(uniprot_df.iloc[:, idx].values)&set(epf_df.index.values) ))
    best_match_epf_name = uniprot_df.select_dtypes(include='object').columns.values[np.argmax(match_epf_num)]
    print('best match epf name: ', best_match_epf_name)
    
    seq = None
    seq_in_exp_idx_list = []
    seq_off_exp_idx_list = []
    for cur_name in list(epf_df.index.values):
        cur_uname_list = uniprot_df.loc[uniprot_df[best_match_epf_name] == cur_name][best_match_seq_name].values.tolist()
        if len(cur_uname_list) > 0:
            if len(set(cur_uname_list)&set(seq_df.columns.values)) > 0:
                seq_in_exp_idx_list.append(epf_name_idx_dict[cur_name])
                cur_seq = np.nan_to_num(seq_df[cur_uname_list].values / seq_df[cur_uname_list].values.sum(0)).mean(1).reshape(1, -1)
                
                if seq is None:
                    seq = cur_seq
                else:
                    seq = np.concatenate([seq, cur_seq], 0)
            else:
                seq_off_exp_idx_list.append(epf_name_idx_dict[cur_name])
        else:
            seq_off_exp_idx_list.append(epf_name_idx_dict[cur_name])
    
    #### Save input of protein sequences
    np.savez('/'.join([input_path, 'seq_FCGR_16x_' + cur_exp_cond]),    seq)
    np.savez('/'.join([input_path, 'seq_in_exp_idx_' + cur_exp_cond]),  seq_in_exp_idx_list)
    np.savez('/'.join([input_path, 'seq_off_exp_idx_' + cur_exp_cond]), seq_off_exp_idx_list)
    
    print('Total seqs: ', epf_df.shape[0])
    print('seqs in exp: ', len(seq_in_exp_idx_list), ', seqs not in exp: ', len(seq_off_exp_idx_list))
    print(', '.join(['Finish preprocess protein sequence', cur_exp_name, cur_exp_cond, 'input dim: ' + str(seq.shape[0]) + ' ' + str(seq.shape[1])]))
    print(' ')


def generate_split(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', cur_root='/FREEPII', 
                   cur_train_test_ratio=7, cur_cv_fold=5, cur_pn_ratio=1):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    if not os.path.exists(input_path):
        os.makedirs(input_path)
    
    epf_path     = '/'.join([cur_root, 'EPF-data', cur_exp_name, cur_exp_cond, cur_exp_cond + '_epf.csv'])
    dict_path    = '/'.join([input_path, 'name_idx_dict_' + cur_exp_cond + '.pickle'])
    complex_path = cur_root + '/Protein complex-data/Complexes_gene_' + cur_species + '_filter.csv'
    
    print('Generate split data of ' + cur_exp_name + ' ' + cur_exp_cond)
    
    #### Read data
    epf_df = pd.read_csv(epf_path).set_index('Gene_name', drop=True)
    
    with open(dict_path, 'rb') as f:
        epf_name_idx_dict = pickle.load(f)
    
    match_order = sum([list(epf_df.index)[i]==list(epf_name_idx_dict.keys())[i] for i in range(epf_df.shape[0])])
    if match_order!=epf_df.shape[0]:
        print('re-create name-idx dictionary')
        epf_name_idx_dict = dict(zip(epf_df.index.values, np.arange(epf_df.shape[0])))
        
        with open(dict_path, 'wb') as f:
            pickle.dump(epf_name_idx_dict, f, pickle.HIGHEST_PROTOCOL)
    
    complex_df = pd.read_csv(complex_path)
    complex_df = complex_df[['ComplexName', 'Gene_name']].drop_duplicates()
    complex_df['ComplexName'] = [i.upper() for i in list(complex_df['ComplexName'].values)]
    Complex_gnames = sorted(set(complex_df['Gene_name']), key=lambda x: x.lower())
    print(', '.join(['train split ratio: ' + str(cur_train_test_ratio) + '0 %',  
                     'CV fold number: ' + str(cur_cv_fold), 'PN ratio: ' + str(cur_pn_ratio)]))
    
    #### All ppi in epf
    all_ppi = list(combinations(epf_df.index, 2))
    all_ppi = [sorted(i, key=lambda x: x.lower()) for i in all_ppi]
    all_ppi = pd.DataFrame(all_ppi, columns =['Gene_name_A', 'Gene_name_B'])
    
    #### filter complex based on epf
    cur_com_df = complex_df.loc[complex_df['Gene_name'].isin(epf_df.index),]
    n1_list = [k for k,v in list(dict(Counter(cur_com_df['ComplexName'])).items()) if v < 2]
    cur_com_df = cur_com_df.loc[~cur_com_df['ComplexName'].isin(n1_list),].reset_index(drop=True)
    
    #### Turn dataframe into dict and count complex subunit
    cur_com_dict = cur_com_df.groupby('ComplexName').agg(lambda x: x.unique().tolist())
    cur_com_dict = dict(zip(cur_com_dict.index, cur_com_dict['Gene_name'].values))
    cur_com_N = pd.DataFrame({'ComplexName': cur_com_dict.keys(), 'Gene_name_N': [len(v) for k, v in cur_com_dict.items()]})
    
    #### Extract pairwise interactions in complex
    cur_com_ppi = [sorted(combinations(list(cur_com_dict.items())[i][1], 2)) for i in range(len(cur_com_dict))]
    cur_com_ppi = pd.DataFrame(data={'ComplexName': cur_com_dict.keys(), 'Gene_name': cur_com_ppi})
    cur_com_ppi = cur_com_ppi.explode('Gene_name')
    cur_com_ppi['Gene_name'] = [sorted(i, key=lambda x: x.lower()) for i in cur_com_ppi['Gene_name']]
    cur_com_ppi[['Gene_name_A', 'Gene_name_B']] = pd.DataFrame(cur_com_ppi['Gene_name'].tolist(), index=cur_com_ppi.index)
    cur_com_ppi = cur_com_ppi[['ComplexName', 'Gene_name_A', 'Gene_name_B']].drop_duplicates()
    cur_com_ppi = pd.merge(cur_com_ppi, cur_com_N, how='left', on='ComplexName').drop_duplicates()
    cur_com_ppi = cur_com_ppi.sort_values(by=['Gene_name_N', 'ComplexName'], key=lambda x: x.str.lower() if 
                                          x.dtype=='object' else x).drop_duplicates(subset=['Gene_name_A', 'Gene_name_B'], keep='first').reset_index(drop=True)
    
    #### Remove complex PPI with no co-peak in epf
    cur_mat_ = epf_df.loc[sorted(set(cur_com_ppi[['Gene_name_A', 'Gene_name_B']].values.reshape(-1)), key=lambda x: x.lower())]
    cur_mat_ = cur_mat_.sort_index()
    epf_name_idx_dict_ = dict(zip(cur_mat_.index.values, np.arange(cur_mat_.shape[0])))
    
    cur_com_ppi_copeak = np.dot(np.array(cur_mat_.values > 0.01, dtype=int), np.array(cur_mat_.T.values > 0.01, dtype=int))
    np.fill_diagonal(cur_com_ppi_copeak, 0)
    
    cur_mask = np.zeros(cur_com_ppi_copeak.shape)
    ppi_idx = list(zip([epf_name_idx_dict_[i] for i in cur_com_ppi['Gene_name_A']], [epf_name_idx_dict_[i] for i in cur_com_ppi['Gene_name_B']]))
    ppi_idx = [sorted(i) for i in ppi_idx]
    cur_mask[np.array(ppi_idx)[:,0], np.array(ppi_idx)[:,1]] = 1
    cur_com_ppi_copeak = np.multiply(cur_com_ppi_copeak, cur_mask)
    
    keep_idx = np.array(list(zip(np.where(cur_com_ppi_copeak > 0)[0], np.where(cur_com_ppi_copeak > 0)[1])))
    cur_com_ppi_copeak = [[list(epf_name_idx_dict_.keys())[keep_idx[i, 0]], 
                           list(epf_name_idx_dict_.keys())[keep_idx[i, 1]]] for i in range(keep_idx.shape[0])]
    cur_com_ppi_copeak = [sorted(i, key=lambda x: x.lower()) for i in cur_com_ppi_copeak]
    cur_com_ppi_copeak = pd.DataFrame(cur_com_ppi_copeak, columns =['Gene_name_A', 'Gene_name_B']).drop_duplicates()
    
    cur_com_ppi = pd.merge(cur_com_ppi, cur_com_ppi_copeak, how='inner', on=['Gene_name_A', 'Gene_name_B'])
    
    
    #### Keep complex with size > 2
    cur_com_dict = cur_com_ppi[['ComplexName', 'Gene_name_A', 'Gene_name_B']].groupby('ComplexName').agg(lambda x: x.unique().tolist())
    cur_com_dict['Gene_name'] = [sorted(set(cur_com_dict['Gene_name_A'][i] + cur_com_dict['Gene_name_B'][i]), key=lambda x: x.lower()) 
                                 for i in range(cur_com_dict.shape[0])]
    cur_com_dict = dict(zip(cur_com_dict.index, cur_com_dict['Gene_name'].values))
    cur_com_N = pd.DataFrame({'ComplexName': cur_com_dict.keys(), 'Gene_name_N': [len(v) for k, v in cur_com_dict.items()]})

    cur_com_ppi = pd.merge(cur_com_ppi, cur_com_N, how='left', on='ComplexName').drop_duplicates()
    cur_com_ppi = cur_com_ppi.loc[cur_com_ppi['Gene_name_N_y'] > 2]
    print(', '.join([
        'Unique complex number: ' + str(len(set(cur_com_ppi['ComplexName']))),
        'Unique complex genes: ' + str(len(set(cur_com_ppi[['Gene_name_A', 'Gene_name_B']].values.reshape(-1)))),
        'All complex PPI: ' + str(cur_com_ppi.shape[0]),
        'All unique complex PPI: ' + str(cur_com_ppi.drop_duplicates(subset=['Gene_name_A', 'Gene_name_B']).shape[0])
    ]))
    
    
    #### Start to split (Reference, Held out, Exp PPIs)
    temp = cur_com_ppi[['Gene_name_A', 'Gene_name_B']]
    temp['Label'] = 1
    com_gnames = sorted(set(temp[['Gene_name_A', 'Gene_name_B']].values.reshape(-1)), key=lambda x: x.lower())

    ## Get all complex interactions
    all_com_ppi = list(combinations(com_gnames, 2))
    all_com_ppi = [sorted(i, key=lambda x: x.lower()) for i in all_com_ppi]
    all_com_ppi = pd.DataFrame(all_com_ppi, columns=['Gene_name_A', 'Gene_name_B']).drop_duplicates()

    shuff = list(range(all_com_ppi.shape[0]))
    random.seed(123)
    random.shuffle(shuff)
    all_com_ppi = all_com_ppi.iloc[shuff,:] # shuffle the interactions

    ## Get the split of interactions
    border = int(np.floor(all_com_ppi.shape[0]*cur_train_test_ratio*0.1))
    train_ppi = all_com_ppi.iloc[:border, ]
    test_ppi = all_com_ppi.iloc[border:, ]

    ## Label interactions
    reference = pd.merge(train_ppi, temp, how='left', on=['Gene_name_A', 'Gene_name_B']).drop_duplicates().reset_index(drop=True)
    held_out = pd.merge(test_ppi, temp, how='left', on=['Gene_name_A', 'Gene_name_B']).drop_duplicates().reset_index(drop=True)
    reference['Label'] = reference['Label'].fillna(0)
    held_out['Label'] = held_out['Label'].fillna(0)

    temp_ = pd.concat([reference, held_out])
    exp_ppi = pd.merge(all_ppi, temp_, how='left', on=['Gene_name_A', 'Gene_name_B']).drop_duplicates()
    exp_ppi = exp_ppi.loc[exp_ppi['Label'].isna()].reset_index(drop=True)

    exp_ppi['Type'] = 'Exp'
    reference['Type'] = 'Ref'
    held_out['Type'] = 'Held_out'
    
    
    #### Split Reference set into fold for cross validation
    ## balance class
    random.seed(123)
    tn_keep_idx = random.sample(list(reference.loc[reference['Label']==0].index), sum(reference['Label']==1)*cur_pn_ratio)
    tn_del_idx = list(set(reference.loc[reference['Label']==0].index) ^ set(tn_keep_idx))

    del_fold_ppi = reference.iloc[tn_del_idx].reset_index(drop=True)
    del_fold_ppi['Type'] = 'Del_fold'

    ## Split into cv folds
    ref_tp = reference.loc[reference['Label']==1].reset_index(drop=True)
    ref_tn = reference.iloc[tn_keep_idx].reset_index(drop=True)

    reference = pd.concat([reference.loc[reference['Label']==1], reference.iloc[tn_keep_idx]]).reset_index(drop=True)
    shuff = list(range(reference.shape[0]))
    random.seed(123)
    random.shuffle(shuff)
    reference = reference.iloc[shuff,:].reset_index(drop=True) # shuffle the interactions

    break_idx = list(range(0, ref_tp.shape[0], int(ref_tp.shape[0] / cur_cv_fold)))

    if len(break_idx)==cur_cv_fold:
        break_idx = break_idx + [ref_tp.shape[0]]
    else:
        break_idx[len(break_idx)-1] = ref_tp.shape[0]

    shuff = list(range(ref_tp.shape[0]))
    random.seed(123)
    random.shuffle(shuff)

    idx_fold_list = [shuff[break_idx[i]:break_idx[i+1]] for i in range(len(break_idx)-1)]
    df_fold_list = [pd.concat([ref_tp.iloc[idx_fold_list[i]], 
                               ref_tn.iloc[idx_fold_list[i]]]).reset_index(drop=True) for i in range(len(idx_fold_list))]
    
    for idx in range(cur_cv_fold):
        grep_idx = list(set([idx])^set(range(cur_cv_fold)))
        
        cur_train = pd.concat([df_fold_list[i] for i in grep_idx]).reset_index(drop=True).iloc[:,:3]
        cur_train.columns= ['Gene_name_A', 'Gene_name_B', 'train_' + str(idx + 1)]
        reference = pd.merge(reference, cur_train, how='left', on=['Gene_name_A', 'Gene_name_B']).drop_duplicates()
        reference['train_' + str(idx + 1)] = np.where(reference['train_' + str(idx + 1)].isna(), 0, 1)

        cur_test = df_fold_list[idx].reset_index(drop=True).iloc[:,:3]
        cur_test.columns= ['Gene_name_A', 'Gene_name_B', 'test_' + str(idx + 1)]
        reference = pd.merge(reference, cur_test, how='left', on=['Gene_name_A', 'Gene_name_B']).drop_duplicates()
        reference['test_' + str(idx + 1)] = np.where(reference['test_' + str(idx + 1)].isna(), 0, 1)
    
    
    #### Save edges, labels and masks
    reference['Edge_A'] = [epf_name_idx_dict[i] for i in reference['Gene_name_A']]
    reference['Edge_B'] = [epf_name_idx_dict[i] for i in reference['Gene_name_B']]

    held_out['Edge_A'] = [epf_name_idx_dict[i] for i in held_out['Gene_name_A']]
    held_out['Edge_B'] = [epf_name_idx_dict[i] for i in held_out['Gene_name_B']]

    del_fold_ppi['Edge_A'] = [epf_name_idx_dict[i] for i in del_fold_ppi['Gene_name_A']]
    del_fold_ppi['Edge_B'] = [epf_name_idx_dict[i] for i in del_fold_ppi['Gene_name_B']]

    exp_ppi['Edge_A'] = [epf_name_idx_dict[i] for i in exp_ppi['Gene_name_A']]
    exp_ppi['Edge_B'] = [epf_name_idx_dict[i] for i in exp_ppi['Gene_name_B']]

    np.savez('/'.join([input_path, 'ref_edge_'     + cur_exp_cond]), reference[['Edge_A', 'Edge_B']].values)
    np.savez('/'.join([input_path, 'heldout_edge_' + cur_exp_cond]), held_out[['Edge_A', 'Edge_B']].values)
    np.savez('/'.join([input_path, 'delfold_edge_' + cur_exp_cond]), del_fold_ppi[['Edge_A', 'Edge_B']].values)
    np.savez('/'.join([input_path, 'exp_edge_'     + cur_exp_cond]), exp_ppi[['Edge_A', 'Edge_B']].values)

    np.savez('/'.join([input_path, 'ref_label_'     + cur_exp_cond]), reference['Label'].values)
    np.savez('/'.join([input_path, 'heldout_label_' + cur_exp_cond]), held_out['Label'].values)
    np.savez('/'.join([input_path, 'delfold_label_' + cur_exp_cond]), del_fold_ppi['Label'].values)

    np.savez('/'.join([input_path, 'train_mask_' + cur_exp_cond]), reference.loc[:,[i for i in reference.columns if 'train' in i]])
    np.savez('/'.join([input_path, 'test_mask_'  + cur_exp_cond]), reference.loc[:,[i for i in reference.columns if 'test' in i]])

    print(', '.join(['TP in reference: ' + str(sum(reference['Label']==1)), 'TP in held-out: ' + str(sum(held_out['Label']==1))]))
    print(', '.join(['TN in reference: ' + str(sum(reference['Label']==0)), 'TN in held-out: ' + str(sum(held_out['Label']==0))]))
    print(', '.join(['All PPIs: ' + str(all_ppi.shape[0]), 'All reference PPIs: ' + str(reference.shape[0]), 
                     'All held-out PPIs: ' + str(held_out.shape[0]), 
                     'All del-fold PPIs: ' + str(del_fold_ppi.shape[0]), 
                     'All exp PPIs: ' + str(exp_ppi.shape[0])]))
    print(' '.join(['Done generate split of', cur_exp_name, cur_exp_cond]))
    print(' ')



#### 1. Pre-process protein complexes ####
process_complex(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name, 
                cur_exp_cond=args.exp_cond, cur_root=args.root)

#### 2. Pre-process elution profile ####
process_epf(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name, 
            cur_exp_cond=args.exp_cond, cur_root=args.root)

#### 3. Pre-process protein sequences ####
process_seq(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name, 
            cur_exp_cond=args.exp_cond, cur_root=args.root)

#### 4. Generate split ####
generate_split(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name,
               cur_exp_cond=args.exp_cond, cur_root=args.root, 
               cur_train_test_ratio=args.train_test_ratio, cur_cv_fold=args.cv_fold, 
               cur_pn_ratio=args.pn_ratio)




