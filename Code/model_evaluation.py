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
from sklearn.metrics import roc_auc_score
from copy import deepcopy
import requests, sys
import seaborn as sns
import matplotlib.pyplot as plt
from goatools.semantic import TermCounts, get_info_content, common_parent_go_ids, min_branch_length, semantic_distance, semantic_similarity, deepest_common_ancestor
from goatools.semantic import resnik_sim, lin_sim
from goatools.associations import read_associations, dnld_assc
from goatools.base import get_godag
from goatools.obo_parser import GODag
from Similarity_organize import *

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

if args.output_path=='/Output':
    args.output_path = '/'.join([args.root + args.output_path, args.exp_name, args.exp_cond])

if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

def classification_performance(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', 
                               cur_root='/FREEPII', output_path='/Output'):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    
    split_path = [input_path + '/' + i for i in os.listdir(input_path) if 'ref' in i or 'heldout' in i 
                  or 'delfold' in i]
    split_path = [i for i in split_path if 'label' in i]
    
    predict_path = output_path + '/predict/'
    predict_path = [predict_path + i for i in os.listdir(predict_path)]
    
    #### Combined labels and predictions
    ref_out = np.concatenate([np.load([i for i in split_path   if 'ref'   in i][0])['arr_0'].reshape(-1, 1),
                              np.load([i for i in predict_path if 'train' in i][0])['arr_0'].reshape(-1, 1)], 1)
    
    heldout_out = np.concatenate([np.load([i for i in split_path   if 'heldout' in i][0])['arr_0'].reshape(-1, 1),
                                  np.load([i for i in predict_path if 'heldout' in i][0])['arr_0'].reshape(-1, 1)], 1)
    
    delfold_out = np.concatenate([np.load([i for i in split_path   if 'delfold' in i][0])['arr_0'].reshape(-1, 1),
                                  np.load([i for i in predict_path if 'delfold' in i][0])['arr_0'].reshape(-1, 1)], 1)
    
    
    #### Confusion matrix
    ref_confumat     = [sum((ref_out[:, 0]==1)&(ref_out[:, 1]  > 0.5)), sum((ref_out[:, 0]==1)&(ref_out[:, 1] <= 0.5)),
                        sum((ref_out[:, 0]==0)&(ref_out[:, 1] <= 0.5)), sum((ref_out[:, 0]==0)&(ref_out[:, 1]  > 0.5))]
    
    heldout_confumat = [sum((heldout_out[:, 0]==1)&(heldout_out[:, 1]  > 0.5)), 
                        sum((heldout_out[:, 0]==1)&(heldout_out[:, 1] <= 0.5)),
                        sum((heldout_out[:, 0]==0)&(heldout_out[:, 1] <= 0.5)), 
                        sum((heldout_out[:, 0]==0)&(heldout_out[:, 1]  > 0.5))]
    
    delfold_confumat = [sum((delfold_out[:, 0]==1)&(delfold_out[:, 1]  > 0.5)), 
                        sum((delfold_out[:, 0]==1)&(delfold_out[:, 1] <= 0.5)),
                        sum((delfold_out[:, 0]==0)&(delfold_out[:, 1] <= 0.5)), 
                        sum((delfold_out[:, 0]==0)&(delfold_out[:, 1]  > 0.5))]
    
    
    #### Evaluate performance
    eval_df = pd.DataFrame([ref_confumat, heldout_confumat, delfold_confumat], columns=['TP', 'FN', 'TN', 'FP'])
    
    eval_df['Sensitivity'] = eval_df['TP'].values / (eval_df['TP'].values + eval_df['FN'].values)
    eval_df['Specificity'] = eval_df['TN'].values / (eval_df['TN'].values + eval_df['FP'].values)
    
    a = (eval_df['Sensitivity'].values * eval_df['Specificity'].values) - ((1-eval_df['Sensitivity'].values) * 
                                                                           (1-eval_df['Specificity'].values))
    b = (np.sqrt(eval_df['Sensitivity'].values + (1-eval_df['Specificity'].values)) * 
         np.sqrt(eval_df['Sensitivity'].values + (1-eval_df['Sensitivity'].values)) * 
         np.sqrt(eval_df['Specificity'].values + (1-eval_df['Specificity'].values)) *
         np.sqrt(eval_df['Specificity'].values + (1-eval_df['Sensitivity'].values)))
    eval_df['MCC'] = a/b
    del a
    del b
    
    eval_df['AUC_ROC'] = [roc_auc_score(ref_out[:,0], ref_out[:,1]), 
                          roc_auc_score(heldout_out[:,0], heldout_out[:,1]), np.NaN]
    eval_df['Type'] = ['Ref', 'Held_out', 'Del_fold']
    
    
    #### Save output and show results
    eval_path = output_path + '/evaluation'
    
    if not os.path.exists(eval_path):
        os.makedirs(eval_path)
    
    eval_df.to_csv(eval_path + '/classification_performance.csv', index=False, header=True)
    print('Output path: ' + eval_path)
    print('Classification performance of ', cur_exp_name + ' ' + cur_exp_cond)
    print(eval_df.iloc[:,4:])


def clustering_composite_score(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', 
                               cur_root='/FREEPII', output_path='/Output'):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    
    dict_path = input_path + '/name_idx_dict_' + cur_exp_cond + '.pickle'
    with open(dict_path, 'rb') as f:
        name_idx_dict = pickle.load(f)
    
    cluster_path = output_path + '/Cluster'
    with open(cluster_path + '/clusters.txt','r') as fh:
        clusters = []
        for oneline in fh:
            oneline = oneline.strip()
            clusters.append([int(i) for i in oneline.split(' ')])
    
    
    #### Prepross gs
    complex_path = cur_root + '/Protein complex-data/Complexes_gene_' + cur_species + '_filter.csv'
    gs_name_df = pd.read_csv(complex_path)
    
    gs_name_df_flat = gs_name_df[['ComplexName', 'Gene_name']].groupby('ComplexName').agg(lambda x: x.unique().tolist())
    gs_name_df_flat.reset_index(inplace=True)
    gs_name_df_flat['Gene_id'] = [[name_idx_dict[j] for j in i if j in name_idx_dict.keys()] 
                                  for i in list(gs_name_df_flat['Gene_name'].values)]
    gs_name_df_flat['NameSize'] = [len(i) for i in list(gs_name_df_flat['Gene_name'].values)]
    gs_name_df_flat['IDSize'] = [len(i) for i in list(gs_name_df_flat['Gene_id'].values)]
    
    gs_name_df_flat_sub = gs_name_df_flat.loc[gs_name_df_flat['IDSize'].values > 2]
    
    
    #### Overlap between gs and cluster
    ove_mat = np.array([[len(set(clusters[i])&set(gs_name_df_flat_sub['Gene_id'].values[j])) 
                         for j in range(gs_name_df_flat_sub.shape[0])] for i in range(len(clusters))])
    
    ove_score_mat = np.array([[(len(set(clusters[i])&set(gs_name_df_flat_sub['Gene_id'].values[j]))**2)/(
        len(set(clusters[i])) * len(set(gs_name_df_flat_sub['Gene_id'].values[j]))) 
                               for j in range(gs_name_df_flat_sub.shape[0])] for i in range(len(clusters))])
    
    
    #### Composite scores
    Overlapp = sum([any(ove_score_mat[i,] >= 0.25) for i in range(ove_score_mat.shape[0])]) / ove_score_mat.shape[0]
    Sensitivity = sum(ove_mat.max(0)) / sum(gs_name_df_flat_sub['IDSize'].values)
    PPV = sum(ove_mat.max(1)) / sum(ove_mat.sum(0))
    Accuracy = np.sqrt(Sensitivity*PPV)
    MMR = sum(ove_score_mat.max(0) / ove_score_mat.shape[1])
    eval_df = pd.DataFrame({'Metrics': ['Overlapp', 'Sensitivity', 'PPV', 'Accuracy', 'MMR'],
                            'Values': [Overlapp, Sensitivity, PPV, Accuracy, MMR]})
    
    
    #### Save and show output
    eval_path = output_path + '/evaluation'
    if not os.path.exists(eval_path):
        os.makedirs(eval_path)
    
    eval_df.to_csv(eval_path + '/cluster_composite_score.csv', index=False, header=True)
    print('Output path: ' + eval_path)
    print('Composite score of ', cur_exp_name + ' ' + cur_exp_cond + ' clusters')
    print(eval_df)


def clustering_colocalization_score(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', 
                                    cur_root='/FREEPII', output_path='/Output'):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    
    dict_path = input_path + '/name_idx_dict_' + cur_exp_cond + '.pickle'
    with open(dict_path, 'rb') as f:
        name_idx_dict = pickle.load(f)
    idx_name_dict = dict(zip(list(name_idx_dict.values()), list(name_idx_dict.keys())))
    
    cluster_path = output_path + '/Cluster'
    with open(cluster_path + '/clusters.txt','r') as fh:
        idx_clusters = []
        for oneline in fh:
            oneline = oneline.strip()
            idx_clusters.append([int(i) for i in oneline.split(' ')])
    name_clusters = [[idx_name_dict[j] for j in i] for i in idx_clusters]
    
    
    #### Sub-locatoin and weight
    subloc_path = cur_root + '/UniProt-data/uniprot_subloc_' + cur_species + '.rds'
    subloc_df = list(pyreadr.read_r(subloc_path).values())[0]
    
    gene_path = cur_root + '/UniProt-data/uniprot_gene_' + cur_species + '.rds'
    gene_df = list(pyreadr.read_r(gene_path).values())[0]
    
    gene_subloc_df = pd.merge(gene_df[['UniProt_entry', 'Gene_name']], subloc_df[['UniProt_entry', 'Location']], 
                              on='UniProt_entry', how='inner')
    gene_subloc_df = gene_subloc_df[['Gene_name', 'Location']].drop_duplicates()
    
    gene_subloc_df_sub = gene_subloc_df.loc[gene_subloc_df['Gene_name'].isin( list(idx_name_dict.values()) ), ]
    
    gene_subloc_df_sub = gene_subloc_df_sub.groupby('Location').agg(lambda x: x.unique().tolist())
    gene_subloc_df_sub.reset_index(inplace=True)
    
    
    #### Calculate co-localization scores
    eval_df = []
    for com_idx in range(len(name_clusters)):
        cur_com = 'Complex_' + str(com_idx+1)
        cur_com_genes = sorted(set(name_clusters[com_idx]), key=lambda x: x.lower())
        
        cur_max_subloc = 0
        cur_max_group_size = 0
        cur_assign_group = []
        
        for subloc_idx in range(gene_subloc_df_sub.shape[0]):
            cur_subloc = gene_subloc_df_sub.iloc[subloc_idx, ]['Location']
            cur_group = [i for i in gene_subloc_df_sub.iloc[subloc_idx, ]['Gene_name'] if i in cur_com_genes]
            cur_group_size = len(cur_group)
            
            if cur_group_size > cur_max_group_size and cur_group_size >= 2:
                cur_max_subloc = cur_subloc
                cur_max_group_size = cur_group_size
            
            if len(cur_group) > 0:
                cur_assign_group += cur_group
        
        cur_assign_group = sorted(set(cur_assign_group), key=lambda x: x.lower())
        
        eval_df.append([cur_com, len(cur_com_genes), cur_max_subloc, cur_max_group_size, 
                        len(cur_assign_group), cur_max_group_size / len(cur_assign_group)])
    
    eval_df = pd.DataFrame(eval_df, columns=['ComplexName', 'ComplexSize', 'Location', 'MaxGroupSize', 
                                             'AssignedGroupSize', 'Score'])
    
    
    #### Save and show output
    eval_path = output_path + '/evaluation'
    if not os.path.exists(eval_path):
        os.makedirs(eval_path)
    
    eval_df.to_csv(eval_path + '/cluster_colocalization_score.csv', index=False, header=True)
    print('Output path: ' + eval_path)
    print('Co-localization score of ', cur_exp_name + ' ' + cur_exp_cond + ' clusters: ')
    print('; '.join(['Max: ' + str(max(eval_df['Score'].values))[:5], 
                     'Min: ' + str(min(eval_df['Score'].values))[:5],
                     'Mean: ' + str(np.mean(eval_df['Score'].values))[:5],
                     'Median: ' + str(np.median(eval_df['Score'].values))[:5]]))


def clustering_GO_score(cur_species='Human', cur_exp_name='PXD002892', cur_exp_cond='SEC2-heavy', 
                        cur_root='/FREEPII', output_path='/Output'):
    input_path = '/'.join([cur_root + '/input', cur_exp_name, cur_exp_cond])
    
    dict_path = input_path + '/name_idx_dict_' + cur_exp_cond + '.pickle'
    with open(dict_path, 'rb') as f:
        name_idx_dict = pickle.load(f)
    idx_name_dict = dict(zip(list(name_idx_dict.values()), list(name_idx_dict.keys())))
    
    cluster_path = output_path + '/Cluster'
    with open(cluster_path + '/clusters.txt','r') as fh:
        idx_clusters = []
        for oneline in fh:
            oneline = oneline.strip()
            idx_clusters.append([int(i) for i in oneline.split(' ')])
    name_clusters = [[idx_name_dict[j] for j in i] for i in idx_clusters]
    
    ## GO annot
    go_path = cur_root + '/GO-data'
    go = get_godag(go_path + '/go-basic.obo', optional_attrs={'relationship'})

    assocs_symbol_BP = read_associations(go_path + '/association_symbol_BP_Human.txt')
    assocs_symbol_MF = read_associations(go_path + '/association_symbol_MF_Human.txt')
    assocs_symbol_CC = read_associations(go_path + '/association_symbol_CC_Human.txt')

    assocs_synonym_BP = read_associations(go_path + '/association_synonym_BP_Human.txt')
    assocs_synonym_MF = read_associations(go_path + '/association_synonym_MF_Human.txt')
    assocs_synonym_CC = read_associations(go_path + '/association_synonym_CC_Human.txt')

    termcounts_BP = TermCounts(go, assocs_synonym_BP)
    termcounts_MF = TermCounts(go, assocs_synonym_MF)
    termcounts_CC = TermCounts(go, assocs_synonym_CC)
    
    
    #### Calculate sementic similarity
    ## Turn complex with go id set
    go_cluster_BP = [[list(assocs_synonym_BP.get(i, {})) for i in j] for j in name_clusters]
    go_cluster_MF = [[list(assocs_synonym_MF.get(i, {})) for i in j] for j in name_clusters]
    go_cluster_CC = [[list(assocs_synonym_CC.get(i, {})) for i in j] for j in name_clusters]
    
    ## Filter complex protein without go-term
    go_cluster_BP = [[i for i in j if i!=[]] for j in go_cluster_BP]
    go_cluster_MF = [[i for i in j if i!=[]] for j in go_cluster_MF]
    go_cluster_CC = [[i for i in j if i!=[]] for j in go_cluster_CC]
    
    ## Calculate similarity of pair-wise go id in cluster
    GOGO_com_score_BP = deepcopy(go_cluster_BP)
    GOGO_com_score_MF = deepcopy(go_cluster_MF)
    GOGO_com_score_CC = deepcopy(go_cluster_CC)
    
    print('Calculating GOGO scores.....it will take a while to run')
    for idx in range(len(go_cluster_BP)):
        if len(go_cluster_BP[idx]) <2:
            cur_score_BP = 0
        elif len(go_cluster_BP[idx])==2:
            cur_score_BP = Similarity_of_Set_of_GOTerms(go_cluster_BP[idx][0], go_cluster_BP[idx][1], go, 'GOGO', silent=True)
        elif len(go_cluster_BP[idx]) >2:
            cur_BP_pair = list(combinations(go_cluster_BP[idx], 2))
            cur_score_BP = np.mean([Similarity_of_Set_of_GOTerms(i[0], i[1], go, 'GOGO', silent=True) for i in cur_BP_pair])
        
        if len(go_cluster_MF[idx]) <2:
            cur_score_MF = 0
        elif len(go_cluster_MF[idx])==2:
            cur_score_MF = Similarity_of_Set_of_GOTerms(go_cluster_MF[idx][0], go_cluster_MF[idx][1], go, 'GOGO', silent=True)
        elif len(go_cluster_MF[idx]) >2:
            cur_MF_pair = list(combinations(go_cluster_MF[idx], 2))
            cur_score_MF = np.mean([Similarity_of_Set_of_GOTerms(i[0], i[1], go, 'GOGO', silent=True) for i in cur_MF_pair])
        
        if len(go_cluster_CC[idx]) <2:
            cur_score_CC = 0
        elif len(go_cluster_CC[idx])==2:
            cur_score_CC = Similarity_of_Set_of_GOTerms(go_cluster_CC[idx][0], go_cluster_CC[idx][1], go, 'GOGO', silent=True)
        elif len(go_cluster_CC[idx]) >2:
            cur_CC_pair = list(combinations(go_cluster_CC[idx], 2))
            cur_score_CC = np.mean([Similarity_of_Set_of_GOTerms(i[0], i[1], go, 'GOGO', silent=True) for i in cur_CC_pair])
        
        GOGO_com_score_BP[idx] = cur_score_BP
        GOGO_com_score_MF[idx] = cur_score_MF
        GOGO_com_score_CC[idx] = cur_score_CC
        
        if idx+1 in [int(np.percentile(range(1, len(go_cluster_BP)+1), i)) for i in [20, 40, 60, 80, 100]]:
            print('Run ' + str(np.round((idx+1)/len(go_cluster_BP), 2)*100) + '% of data')
    
    #### Plot
    char = 'Sementic similarity (Complex, GOGO)'
    df = pd.Series(GOGO_com_score_BP, name = 'BP').to_frame()
    df = df.join(pd.Series(GOGO_com_score_MF, name = 'MF'))
    df = df.join(pd.Series(GOGO_com_score_CC, name = 'CC'))
    fig, ax = plt.subplots(figsize=(14, 7))
    sns.boxplot(data=df, width = 0.5, ax=ax).set(title=char)
    # plt.show()
    
    
    #### Save and show output
    eval_path = output_path + '/evaluation/Go_score'
    if not os.path.exists(eval_path):
        os.makedirs(eval_path)
    
    np.save(eval_path + '/GOGO_com_score_BP', GOGO_com_score_BP)
    np.save(eval_path + '/GOGO_com_score_MF', GOGO_com_score_MF)
    np.save(eval_path + '/GOGO_com_score_CC', GOGO_com_score_CC)
    fig.savefig(eval_path + '/GOGO_com_score.png')
    
    print('Output path: ' + eval_path)
    print('GOGO score of ', cur_exp_name + ' ' + cur_exp_cond + ' clusters: ')
    print('; '.join(['Max BP score: '    + str(max(GOGO_com_score_BP))[:5], 
                     'Min BP score: '    + str(min(GOGO_com_score_BP))[:5],
                     'Mean BP score: '   + str(np.mean(GOGO_com_score_BP))[:5],
                     'Median BP score: ' + str(np.median(GOGO_com_score_BP))[:5]]))
    print('; '.join(['Max MF score: '    + str(max(GOGO_com_score_MF))[:5], 
                     'Min MF score: '    + str(min(GOGO_com_score_MF))[:5],
                     'Mean MF score: '   + str(np.mean(GOGO_com_score_MF))[:5],
                     'Median MF score: ' + str(np.median(GOGO_com_score_MF))[:5]]))
    print('; '.join(['Max CC score: '    + str(max(GOGO_com_score_CC))[:5], 
                     'Min CC score: '    + str(min(GOGO_com_score_CC))[:5],
                     'Mean CC score: '   + str(np.mean(GOGO_com_score_CC))[:5],
                     'Median CC score: ' + str(np.median(GOGO_com_score_CC))[:5]]))


#### 1. Classification performance ####
classification_performance(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name, 
                           cur_exp_cond=args.exp_cond, 
                           cur_root=args.root, output_path=args.output_path)

#### 2. Clustering evaluation: composite score ####
clustering_composite_score(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name, 
                           cur_exp_cond=args.exp_cond, 
                           cur_root=args.root, output_path=args.output_path)

#### 3. Clustering evaluation: co-localization score ####
clustering_colocalization_score(cur_species=args.species.capitalize(), 
                                cur_exp_name=args.exp_name, cur_exp_cond=args.exp_cond, 
                                cur_root=args.root, output_path=args.output_path)

#### 4. Clustering evaluation: GOGO score ####
clustering_GO_score(cur_species=args.species.capitalize(), cur_exp_name=args.exp_name, 
                    cur_exp_cond=args.exp_cond, 
                    cur_root=args.root, output_path=args.output_path)



