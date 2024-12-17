#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import pickle
import math
import os
import re
import sys
import pyreadr
import random
import gc
import argparse
import shutil
from itertools import combinations, chain
from collections import OrderedDict, Counter
import matplotlib.pyplot as plt
from matplotlib.pylab import show, cm, axis
import torch
from torch import nn, einsum
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from einops import rearrange, repeat, reduce
from einops.layers.torch import Rearrange
from torchmetrics import Accuracy
from copy import deepcopy
from dynamicTreeCut import cutreeHybrid
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import networkx as nx
import seaborn as sns
np.set_printoptions(edgeitems=10, linewidth=400)
accuracy = Accuracy(task="binary")

parser = argparse.ArgumentParser()
parser.add_argument('-root',        '--root',             default='/FREEPII',   help='root',                          type=str)
parser.add_argument('-sp',          '--species',          default='Human',      help='species',                       type=str)
parser.add_argument('-e_name',      '--exp_name',         default='PXD002892',  help='exp_name',                      type=str)
parser.add_argument('-e_cond',      '--exp_cond',         default='SEC2-heavy', help='exp_condition',                 type=str)
parser.add_argument('-split_ratio', '--train_test_ratio', default='7',          help='train_test_ratio',              type=int)
parser.add_argument('-cv_fold',     '--cv_fold',          default='5',          help='cv_fold_number',                type=int)
parser.add_argument('-pn_ratio',    '--pn_ratio',         default='1',          help='pos_to_neg_ratio',              type=int)
parser.add_argument('-out_path',    '--output_path',      default='/output',    help='output_path',                   type=str)
parser.add_argument('-pretrain_w',  '--pretrain_w',       default='null',       help='path_to_pretrain_model_weight', type=str)
args = parser.parse_args()
# args = parser.parse_args([])

if '\\' in os.getcwd():
    args.root = ['/'.join(os.getcwd().split('\\')[:(i+1)]) for i,v in enumerate(os.getcwd().split('\\')) 
                 if v=='FREEPII'][0]

args.root = ['/'.join(os.getcwd().split('/')[:(i+1)]) for i,v in enumerate(os.getcwd().split('/')) if v=='FREEPII'][0]

if args.output_path=='/output':
    args.output_path = '/'.join([args.root + args.output_path, args.exp_name, args.exp_cond])

if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

print('Output path: ' + args.output_path)

DEVICE = 'cpu' # or 'cuda:0'
emb_size = [200, 256]


#### Read data

## Input
input_dir = '/'.join([args.root + '/input', args.exp_name, args.exp_cond])

train_edge    = torch.tensor( np.load(input_dir + '/ref_edge_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
heldout_edge  = torch.tensor( np.load(input_dir + '/heldout_edge_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
delfold_edge  = torch.tensor( np.load(input_dir + '/delfold_edge_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
exp_edge      = torch.tensor( np.load(input_dir + '/exp_edge_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
train_label   = torch.tensor( np.load(input_dir + '/ref_label_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
heldout_label = torch.tensor( np.load(input_dir + '/heldout_label_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
delfold_label = torch.tensor( np.load(input_dir + '/delfold_label_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
train_mask    = torch.tensor( np.load(input_dir + '/train_mask_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.bool).to(DEVICE)
test_mask     = torch.tensor( np.load(input_dir + '/test_mask_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.bool).to(DEVICE)

seq_in_exp_idx     = torch.tensor( np.load(input_dir + '/seq_in_exp_idx_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
seq_off_exp_idx    = torch.tensor( np.load(input_dir + '/seq_off_exp_idx_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.long).to(DEVICE)
seq                = torch.tensor( np.load(input_dir + '/seq_FCGR_' + str(int(emb_size[1]**0.5)) + 
                                           'x_' + args.exp_cond + '.npz')['arr_0'] ).to(torch.float32)
seq                = torch.unsqueeze(seq, -1).to(DEVICE)

## EPF
epf_path = '/'.join([args.root + '/EPF-data', args.exp_name, args.exp_cond, args.exp_cond + '_epf.csv'])
epf_df = pd.read_csv(epf_path)
epf = epf_df.select_dtypes(include=np.number).values
pad = np.zeros((epf.shape[0], emb_size[0] - epf.shape[1]))
epf = torch.tensor( np.concatenate([epf, pad], 1).reshape(-1, emb_size[0], 1)).to(torch.float32).to(DEVICE)

## Protein complex
complex_path = args.root + '/Protein complex-data/Complexes_gene_' + args.species.capitalize() + '_filter.csv'
gs = pd.read_csv(complex_path)

## Name idx dictionaray
dict_path = input_dir + '/name_idx_dict_' + args.exp_cond + '.pickle'
with open(dict_path, 'rb') as f:
    name_idx_dict = pickle.load(f)


#### Model layer
class embedding_layer(nn.Module):
    def __init__(self, dim_1, dim_2):
        super(embedding_layer, self).__init__()
        self.Emb = nn.Embedding(dim_1, dim_2, max_norm=True)
    
    def forward(self, in_idx, off_idx, in_emb):
        in_emb = in_emb.squeeze(-1)
        all_tokens = len(list(set(in_idx.tolist() + off_idx.tolist())))
        out_emb = self.Emb(torch.arange(all_tokens))
        
        past_off_idx = []
        for idx in range(all_tokens):
            if idx in in_idx:
                out_emb[idx, :] += in_emb[(idx - len(past_off_idx)), ]
            else:
                past_off_idx.append(idx)
        return out_emb.unsqueeze(-1)

class Encoder(nn.Module):
    def __init__(self, hid_dim=8, dropout=0.1, head_num=8, head_dim=32, max_size=emb_size*2):
        super().__init__()
        
        self.hid_dim = hid_dim
        
        self.conv = nn.Conv1d(1, hid_dim, 3, padding=1)
        self.fc1 = nn.Linear(max_size*(hid_dim + 1), 256)
        self.fc2 = nn.Linear(256, 64)
        self.fc3 = nn.Linear(64, 1)
        
        self.dropout = dropout
        self.dropout_layer = nn.Dropout(p=self.dropout)
        
    def forward(self, inp, e_idx, bias=None):
        x = self.conv(inp.permute(0, 2, 1)).permute(0, 2, 1) + inp.repeat(1, 1, self.hid_dim)
        
        x_e1 = inp[e_idx[:, 0]] - inp[e_idx[:, 1]]
        x_e2 = x[e_idx[:, 0]] - x[e_idx[:, 1]]
        x_e = torch.cat((x_e1 , x_e2), -1).flatten(start_dim=1)
        
        x_e = self.dropout_layer(F.relu(self.fc1(x_e)))
        x_e = self.dropout_layer(F.relu(self.fc2(x_e)))
        w = torch.sigmoid(self.fc3(x_e)).squeeze()
        feat = x.flatten(start_dim=1)
        
        if bias != None:
            feat += bias
        
        return feat, w

class Encoderr(nn.Module):
    def __init__(self, hid_dim=16, max_size=emb_size, dropout=0.1, all_tokens=20):
        super().__init__()
        self.dropout = dropout
        
        self.emb = embedding_layer(all_tokens, emb_size[1])
        self.enc = Encoder(hid_dim=hid_dim, max_size=sum(max_size), dropout=dropout)
        
    def forward(self, inp_e, inp_s, e_idx, in_list=None, off_list=None):
        
        if off_list is not None:
            inp_s = self.emb(in_list, off_list, inp_s)
        
        inp = torch.cat([inp_e, inp_s], 1)
        x, weight = self.enc(inp, e_idx, bias=None)
        
        return x, weight, inp[:, inp_e.shape[1]:, :]


#### Train
predict_weight    = 2.5
L2_weight         = 0.25
weight_decay      = 2.0e-3

if args.pretrain_w=='null':
    record_path = args.output_path + '/performance_record'
    if not os.path.exists(record_path):
        os.makedirs(record_path)
    # print(record_path)
    
    model_w_path = args.output_path + '/best_model_weight'
    if not os.path.exists(model_w_path):
        os.makedirs(model_w_path)
    # print(model_w_path)

predict_path = args.output_path + '/predict'
if not os.path.exists(predict_path):
      os.makedirs(predict_path)
# print(predict_path)

if args.pretrain_w=='null':
    print('Start training FREEPII, data: ' + args.exp_name + ' ' + args.exp_cond)
    
    best_test_loss  = float('inf')
    train_loss_list = []
    test_loss_list  = []
    train_acc_list  = []
    test_acc_list   = []
    
    for cv_idx in range(args.cv_fold):
        model = Encoderr(hid_dim=16, max_size=emb_size, dropout=0.3, all_tokens=max(seq_in_exp_idx.tolist() + 
                                                                                    seq_off_exp_idx.tolist()) + 1).to(DEVICE)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.002, weight_decay=weight_decay)
        
        print('cur cv: ', cv_idx+1)
        train_loss_list.append([])
        test_loss_list.append([])
        train_acc_list.append([])
        test_acc_list.append([])
        best_test_loss_ = float('inf')
        
        for epoch in range(100):
            #### Train =================================================================================================
            model.train()
            optimizer.zero_grad()
            
            _, train_label_predict, _ = model(epf, seq, train_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
            
            cur_train_loss = F.binary_cross_entropy(train_label_predict[train_mask[:,cv_idx]],
                                                    train_label[train_mask[:,cv_idx]].to(torch.float32))
            cur_L2_loss = sum([ (v**2).sum()/2 for v in model.parameters() ])
            
            loss = predict_weight*cur_train_loss + weight_decay*L2_weight*cur_L2_loss
            
            loss.backward()
            optimizer.step()
            
            cur_train_acc = accuracy(train_label_predict[train_mask[:,cv_idx]] > 0.5, train_label[train_mask[:,cv_idx]] > 0.5)
            
            #### Validataion ===========================================================================================
            model.eval()
            _, train_label_predict, _ = model(epf, seq, train_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
            
            cur_test_loss = F.binary_cross_entropy(train_label_predict[test_mask[:,cv_idx]],
                                                   train_label[test_mask[:,cv_idx]].to(torch.float32))
            cur_test_acc = accuracy(train_label_predict[test_mask[:,cv_idx]] > 0.5, train_label[test_mask[:,cv_idx]] > 0.5)
            
            if epoch in [0, 49, 99]:
                print("Epoch:", '%04d' % (epoch + 1), "all_loss=", "{:.4f}".format(loss))
                print("train_loss=", "{:.4f}".format(predict_weight*cur_train_loss), 
                      "test_loss=", "{:.4f}".format(predict_weight*cur_test_loss),
                      "lossL2=", "{:.4f}".format(weight_decay*L2_weight*cur_L2_loss))
            
            if cur_test_loss < best_test_loss_:
                torch.save(model.state_dict(), model_w_path + '/cv' + str(cv_idx+1) + '_stat_dict')
                best_test_loss_ = cur_test_loss
            
            if cur_test_loss < best_test_loss:
                torch.save(model.state_dict(), model_w_path + '/best_stat_dict')
                best_test_loss = cur_test_loss
            
            train_loss_list[cv_idx].append(float(cur_train_loss))
            test_loss_list[cv_idx].append(float(cur_test_loss))
            train_acc_list[cv_idx].append(cur_train_acc)
            test_acc_list[cv_idx].append(cur_test_acc)
        
        #### Evaluation
        model_dict = torch.load(model_w_path + '/cv' + str(cv_idx+1) + '_stat_dict')
        model_load = Encoderr(hid_dim=16, max_size=emb_size, dropout=0.3, 
                              all_tokens=max(seq_in_exp_idx.tolist() + seq_off_exp_idx.tolist()) + 1).to(DEVICE)
        model_load.load_state_dict(model_dict)
        model_load.eval()
        
        print('Evaluation: ')
        with torch.no_grad():
            _, cur_train_pred, _   = model_load(epf, seq, train_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
            _, cur_heldout_pred, _ = model_load(epf, seq, heldout_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
            _, cur_del_pred, _     = model_load(epf, seq, delfold_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
        
        train_tp = (((cur_train_pred > 0.5).to(torch.float32) + 
                     (train_label > 0.5).to(torch.float32))==2).sum() / (train_label > 0.5).sum()
        train_tn = (((cur_train_pred <= 0.5).to(torch.float32) + 
                     (train_label <= 0.5).to(torch.float32))==2).sum() / (train_label <= 0.5).sum()
        print('train TP: ', "{:.4f}".format(float(train_tp)), ' train TN: ', "{:.4f}".format(float(train_tn)))
        
        heldout_tp = (((cur_heldout_pred > 0.5).to(torch.float32) + 
                       (heldout_label > 0.5).to(torch.float32))==2).sum() / (heldout_label > 0.5).sum()
        heldout_tn = (((cur_heldout_pred <= 0.5).to(torch.float32) + 
                       (heldout_label <= 0.5).to(torch.float32))==2).sum() / (heldout_label <= 0.5).sum()
        print('heldout TP: ', "{:.4f}".format(float(heldout_tp)), ' heldout TN: ', "{:.4f}".format(float(heldout_tn)))
        
        del_fp = (cur_del_pred > 0.5).sum() / delfold_label.shape[0]
        del_tn = (((cur_del_pred <= 0.5).to(torch.float32) + 
                   (delfold_label <= 0.5).to(torch.float32))==2).sum() / delfold_label.shape[0]
        print('del-fold FP: ', "{:.4f}".format(float(del_fp)), ' del-fold TN: ', "{:.4f}".format(float(del_tn)))
        
        if cv_idx==0:
            train_pred   = cur_train_pred.reshape(1, -1)
            heldout_pred = cur_heldout_pred.reshape(1, -1)
            del_pred     = cur_del_pred.reshape(1, -1)
        else:
            train_pred   = torch.cat((train_pred, cur_train_pred.reshape(1, -1)), 0)
            heldout_pred = torch.cat((heldout_pred, cur_heldout_pred.reshape(1, -1)), 0)
            del_pred     = torch.cat((del_pred, cur_del_pred.reshape(1, -1)), 0)
        
        del cur_train_pred
        del cur_heldout_pred
        del cur_del_pred
        del model, model_load, model_dict
        print('*'*50)
        gc.collect()  # Collect unused memory
    
    #### Mean predictions over cv folds
    train_pred   = train_pred.mean(0)
    heldout_pred = heldout_pred.mean(0)
    del_pred     = del_pred.mean(0)
    
    np.savez(predict_path + '/train_out',   train_pred.detach().numpy())
    np.savez(predict_path + '/heldout_out', heldout_pred.detach().numpy())
    np.savez(predict_path + '/delfold_out', del_pred.detach().numpy())
    print('Done training FREEPII')
    print(' ')


#### Plot loss
if args.pretrain_w=='null':
    if torch.is_tensor( train_loss_list[0][0] ):
        train_loss_list = [[float(i) for i in train_loss_list[j]] for j in range(len(train_loss_list))]
        # print(torch.is_tensor( train_loss_list[0][0] ))
    if torch.is_tensor( test_loss_list[0][0] ):
        test_loss_list = [[float(i) for i in test_loss_list[j]] for j in range(len(test_loss_list))]
        # print(torch.is_tensor( test_loss_list[0][0] ))
    
    train_acc_list = np.array(train_acc_list)
    test_acc_list = np.array(test_acc_list)
    train_loss_list = np.array(train_loss_list)
    test_loss_list = np.array(test_loss_list)
    np.save(record_path + '/train_loss', train_loss_list)
    np.save(record_path +  '/test_loss',  test_loss_list)
    np.save(record_path +  '/train_acc',  train_acc_list)
    np.save(record_path +   '/test_acc',   test_acc_list)
    
    ## Plot
    fig, ax = plt.subplots(2,2, figsize=(13, 10))
    ax[0,0].plot(np.arange(100), train_acc_list.T)
    ax[0,0].legend(['cv'+str(i+1) for i in range(5)])
    ax[0,0].set_title('Train acc')
    ax[0,1].plot(np.arange(100), test_acc_list.T)
    ax[0,1].legend(['cv'+str(i+1) for i in range(5)])
    ax[0,1].set_title('Testing acc')
    ax[1,0].plot(np.arange(100), train_loss_list.T)
    ax[1,0].legend(['cv'+str(i+1) for i in range(5)])
    ax[1,0].set_title('Train loss')
    ax[1,1].plot(np.arange(100), test_loss_list.T)
    ax[1,1].legend(['cv'+str(i+1) for i in range(5)])
    ax[1,1].set_title('Testing loss')
    plt.show()
    fig.savefig(record_path + '/evaluation.png')
    print('Save evaluation plot')
    print(' ')


#### Prediction
print('Start predictions')
if args.pretrain_w=='null':
    model_dict = torch.load(model_w_path + '/best_stat_dict')
    print('Load the best model weights in this training')
else:
    model_dict = torch.load(args.pretrain_w)
    print('Load the best model weight from: ' + args.pretrain_w)

model_load = Encoderr(hid_dim=16, max_size=emb_size, dropout=0.3, 
                      all_tokens=max(seq_in_exp_idx.tolist() + seq_off_exp_idx.tolist()) + 1).to(DEVICE)
model_load.load_state_dict(model_dict)
model_load.eval()


all_edge = torch.cat([train_edge, heldout_edge, delfold_edge, exp_edge], 0)
label_edge = torch.cat([train_edge, heldout_edge, delfold_edge], 0)


#### Save all predictions in one process if you have enough memory
# emb, all_pred, _ = model_load(epf, seq, all_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
# np.savez(predict_path + '/all_emb_out', emb.detach().numpy())
# np.savez(predict_path + '/all_out', all_pred)


#### Save all predictions by batch
## Emb
with torch.no_grad():
    emb, label_pred, _ = model_load(epf, seq, label_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
np.savez(predict_path + '/all_emb_out', emb.detach().numpy())


## Batch predictions for unlabeled edges
tmp_dir = predict_path + '/tmp/' # Create temporary directory to store intermediate batch results
os.makedirs(tmp_dir, exist_ok=True)

batch_size = 100000  # Adjust the batch size according to your memory capacity

break_box = np.arange(0, exp_edge.shape[0], batch_size).tolist() + [exp_edge.shape[0]]
with torch.no_grad():
    for idx in range(1, len(break_box)):
        batch_edge = exp_edge[break_box[idx-1]:(break_box[idx]),]
        _, exp_pred, _ = model_load(epf, seq, batch_edge, in_list=seq_in_exp_idx, off_list=seq_off_exp_idx)
        
        # Save intermediate results for predictions
        np.savez(os.path.join(tmp_dir, f'exp_pred_batch_{idx}.npz'), exp_pred=exp_pred.detach().cpu().numpy())
        
        if idx+1 in [int(np.percentile(range(1, len(break_box)+1), i)) for i in [20, 40, 60, 80, 100]]:
            print('Processed ' + str(np.round((idx+1)/len(break_box), 2)*100) + '% of batch')
            # print(f"Processed batch {idx} / {len(break_box)-1}")

# Load all intermediate results and concatenate them
# Concatenate all_pred batches
exp_pred_list = []
for idx in range(1, len(break_box)):
    exp_pred_file = os.path.join(tmp_dir, f'exp_pred_batch_{idx}.npz')
    exp_pred_list.append(np.load(exp_pred_file)['exp_pred'])
    
    if idx+1 in [int(np.percentile(range(1, len(break_box)+1), i)) for i in [20, 40, 60, 80, 100]]:
        print('Loaded ' + str(np.round((idx+1)/len(break_box), 2)*100) + '% of batch')
        # print(f"Loaded exp_pred batch {idx} / {len(break_box)-1}")

exp_pred = np.concatenate(exp_pred_list, axis=0)
all_pred = np.concatenate([label_pred.detach().numpy(), exp_pred], axis=0)
np.savez(predict_path + '/all_out', all_pred)
del exp_pred_list
shutil.rmtree(tmp_dir)
print('Done predictions')
print(' ')


#### Clustering
class MCL_Cluster(object):
    def __init__(self, merge_threshold=0.25, min_complex_size=3, max_complex_size=100, split_depth=1):
        self.merge_threshold = merge_threshold
        self.min_complex_size = min_complex_size
        self.max_complex_size = max_complex_size
        self.split_depth = split_depth
        self.break_num = 0
        self.merge_num = 0
        self.max_break_num = 0
        self.max_merge_num = 7
    
    def create_edge_adj(self, e_idx, e_w=None, dim=None, diag=None):
        self.dim = dim
        mat = np.zeros([dim, dim])
        if e_w is not None:
            if e_idx.shape[0]==2:
                mat[e_idx[0], e_idx[1]] = e_w
            elif e_idx.shape[1]==2:
                mat[e_idx[:,0], e_idx[:,1]] = e_w
            else:
                print('wrong edge dimension')
                return
        else:
            if e_idx.shape[0]==3:
                mat[e_idx[0], e_idx[1]] = e_idx[2]
            elif e_idx.shape[1]==3:
                mat[e_idx[:,0], e_idx[:,1]] = e_idx[:,2]
            elif e_idx.shape[0]==2:
                mat[e_idx[0], e_idx[1]] = 1.0
                print('Not give edge weight, use 1.0 as weight')
            elif e_idx.shape[1]==2:
                mat[e_idx[:,0], e_idx[:,1]] = 1.0
                print('Not give edge weight, use 1.0 as weight')
        mat += mat.transpose()
        if diag is not None:
            np.fill_diagonal(mat, diag)
        return mat
    
    def create_tom_adj(self, A):
        d = A.shape[0]
        L = A.dot(A.T)
        K = A.sum(axis=1)

        A_tom = np.zeros_like(A)
        for i in range(d):
            for j in range(i+1, d):  
                numerator = L[i, j] + A[i, j]
                denominator = min(K[i], K[j]) + 1 - A[i, j]
                A_tom[i, j] = numerator / denominator

        A_tom += A_tom.T
        np.fill_diagonal(A_tom, 0)
        return np.nan_to_num(A_tom)
    
    def iterate(self, A, expand_power, inflate_power):
        A_ = np.linalg.matrix_power(A, expand_power)  # expand
        A_ = np.power(A_, inflate_power)              # inflate
        A_ = A_ / A_.sum(0)                           # col-wise normalization
        AA = A_.copy()
        A_[A_ < 0.001] = 0                            # pruning
        num_cols = AA.shape[1]
        col_indices = np.arange(num_cols)
        row_indices = A.argmax(axis=0).reshape((num_cols,)) # keep max of each col in original matrix
        A_[row_indices, col_indices] = AA[row_indices, col_indices]
        return A_
    
    def MCL_process(self, A, expand_power=2, inflate_power=2, iters=2):
        A_ = A / A.sum(0)
        for j in range(iters):
            A_ = self.iterate(A_, expand_power, inflate_power)
        np.fill_diagonal(A_, 0)
        return A_
    
    def label_nodes(self, n_idx, lab):
        sidx = lab.argsort()                                    # Get sidx (sorted indices) for label
        split_idx = np.flatnonzero(np.diff(lab[sidx]) > 0) + 1  # Get where the sorted version of label changes groups
        groups = np.split(n_idx[sidx], split_idx)               # Sort input based on the sidx and split label on split_idx
        return [tuple(sorted(i)) for i in groups]
    
    def remove_subset(self, Cs):
        return list(filter(lambda f: not any(set(f) < set(g) for g in Cs), Cs))
    
    def get_modularity(self, A, Cs, gamma=1.0):
        Lc_m = [( np.array(list(combinations(Cs[j], 2))).shape[0]  / ((np.array(A.nonzero()).T.shape[0] - A.shape[0])/2) ) 
                for j in range(len(Cs))]
        kc_2m = [(( sum([A[k, :].sum() for k in Cs[j]]) / (np.array(A.nonzero()).T.shape[0] - A.shape[0]) )**2)*gamma 
                 for j in range(len(Cs))]
        M = np.array(Lc_m)-np.array(kc_2m)
        return M
    
    def extract_clusters(self, A, keep_idx=None):
        print('extract clusters')
        
        zero_idx = np.array(list(set(np.where(A.sum(0)!=0)[0])&set(np.where(A.sum(1)!=0)[0])))
        
        I = self.create_tom_adj(np.power(A, 2.5))
        
        if keep_idx is None:
            keep_idx = zero_idx
        else:
            keep_idx = np.array(list( set(keep_idx) & set(zero_idx) ))
        
        A_ = A[keep_idx, :][:, keep_idx]
        A_ = self.MCL_process(A_, iters=3)
        I_ = I[keep_idx, :][:, keep_idx]
        
        imp_1 = A_*0.3 + I_*0.7
        imp_2 = A_*0.1 + I_*0.9
        
        dist_1 = pdist(imp_1, metric='cosine')
        dist_2 = pdist(imp_2, metric='cosine')
        
        dist_1[~np.isfinite(dist_1)] = 1.0
        dist_2[~np.isfinite(dist_2)] = 1.0
        
        link_1 = linkage(dist_1, 'ward')
        link_2 = linkage(dist_2, 'ward')
        
        Cs_lab_1 = cutreeHybrid(link_1, dist_1, minClusterSize=(self.min_complex_size), 
                                deepSplit=self.split_depth, respectSmallClusters=True, verbose=False, pamStage=True)
        Cs_lab_2 = cutreeHybrid(link_2, dist_2, minClusterSize=(self.min_complex_size), 
                                deepSplit=self.split_depth, respectSmallClusters=True, verbose=False, pamStage=True)
        Cs_lab_1 = Cs_lab_1['labels']
        Cs_lab_2 = Cs_lab_2['labels']
        
        Cs_1 = self.label_nodes(keep_idx, Cs_lab_1)
        Cs_2 = self.label_nodes(keep_idx, Cs_lab_2)
        Cs = Cs_1 + Cs_2
        
        small_Cs = [j for j in Cs if (len(j) >= self.min_complex_size)&(len(j) <= self.max_complex_size)]
        large_Cs = [j for j in Cs if len(j) > self.max_complex_size]
        
        while len(large_Cs) > 0 and self.break_num <= self.max_break_num:
            print('break clusters, break iter: ', self.break_num)
            
            idx = 0
            max_idx = len(large_Cs)
            max_idx_ori = len(large_Cs)
            while (idx < max_idx)&(idx < max_idx_ori):
                cur_large_C = large_Cs[idx]
                
                keep_idx = np.array(list( set(cur_large_C) & set(zero_idx) ))
                
                A_ = A[keep_idx, :][:, keep_idx]
                A_ = self.MCL_process(A_, iters=5)
                I_ = I[keep_idx, :][:, keep_idx]
                
                imp_1 = A_*0.3 + I_*0.7
                imp_2 = A_*0.1 + I_*0.9

                dist_1 = pdist(imp_1, metric='cosine')
                dist_2 = pdist(imp_2, metric='cosine')

                dist_1[~np.isfinite(dist_1)] = 1.0
                dist_2[~np.isfinite(dist_2)] = 1.0

                link_1 = linkage(dist_1, 'ward')
                link_2 = linkage(dist_2, 'ward')

                Cs_lab_1 = cutreeHybrid(link_1, dist_1, minClusterSize=(self.min_complex_size), 
                                        deepSplit=self.split_depth, respectSmallClusters=True, verbose=False, pamStage=True)
                Cs_lab_2 = cutreeHybrid(link_2, dist_2, minClusterSize=(self.min_complex_size), 
                                        deepSplit=self.split_depth, respectSmallClusters=True, verbose=False, pamStage=True)
                Cs_lab_1 = Cs_lab_1['labels']
                Cs_lab_2 = Cs_lab_2['labels']

                Cs_1 = self.label_nodes(keep_idx, Cs_lab_1)
                Cs_2 = self.label_nodes(keep_idx, Cs_lab_2)
                Cs = Cs_1 + Cs_2
                
                small_Cs = small_Cs + [j for j in Cs if (len(j) >= self.min_complex_size)&(len(j) <= self.max_complex_size)]
                large_Cs = large_Cs + [j for j in Cs if len(j) > self.max_complex_size]
                large_Cs = self.remove_subset(large_Cs)
                idx += 1
                max_idx = len(large_Cs)
            self.break_num += 1
        self.break_num = 0
        # small_Cs = self.remove_subset(small_Cs)
        return small_Cs
    
    def get_ove_score(self, list1, list2):
        if type(list1[0])!=list and type(list1[0])!=tuple:
            list1 = [list1]
        if type(list2[0])!=list and type(list2[0])!=tuple:
            list2 = [list2]
        
        s1 = np.array([[len(set(list1[k])&set(tuple(j)))**2 for j in list2 ] for k in range(len(list1))])
        s2 = np.array([[len(list1[k])*len(j) for j in list2 ] for k in range(len(list1))])
        ove_score = np.nan_to_num(s1 / s2)
        if ove_score.shape[0] > 1:
            np.fill_diagonal(ove_score, 0)
            return ove_score
        else:
            return ove_score.reshape(-1)
    
    def merge_clusters(self, Cs, A):
        if_merge = True
        while if_merge:
            print('merge clusters, merge iter: ', self.merge_num)
            ove_score = self.get_ove_score(Cs, Cs)
            ove_score = np.triu(ove_score, k=1) # keep only upper traingular values
            
            ## Combind cluster id with their overlap score
            ove_clusters_pair_id_w = list(zip(zip(np.arange(ove_score.shape[0]), ove_score.argmax(1)), ove_score.max(1))) 
            # [(cluster_id1, cluster_id2), ove_score]
            ove_clusters_pair_id_w = [j for j in ove_clusters_pair_id_w if j[1] >= self.merge_threshold]                  
            # keep those with overlap score >= threshold
            
            if len(ove_clusters_pair_id_w) < 1:
                if_merge = False
            else:
                ## Index of clusters with similarity highr than a threshold
                ove_clusters_id = list(set([j for k in [i[0] for i in ove_clusters_pair_id_w] for j in k]))
                
                ## Sorted by overlap score
                ove_clusters_pair_id_w = sorted(ove_clusters_pair_id_w, key=lambda x: x[1], reverse=True)
                
                ## Combind cluster id, overlap score, and overlap size
                ove_clusters_pair_id_w_size = [ove_clusters_pair_id_w[idx] + 
                                               ( len(set(Cs[ove_clusters_pair_id_w[idx][0][0]] + 
                                                         Cs[ove_clusters_pair_id_w[idx][0][1]])), ) 
                                               for idx in range(len(ove_clusters_pair_id_w))] # [(cluster_id1, cluster_id2), ove_score, ove_size]
                
                ## filter by max size
                ove_clusters_pair_id_w_size = [j for j in ove_clusters_pair_id_w_size if j[2] <= self.max_complex_size]
                
                if len(ove_clusters_pair_id_w_size) < 1:
                    if_merge = False
                else:
                    ## combind pairs to set
                    ove_clusters_pair_id_w_size_set = [[j for j in ove_clusters_pair_id_w_size if k in j[0]] for k in ove_clusters_id]
                    ove_clusters_pair_id_w_size_set = [j for j in ove_clusters_pair_id_w_size_set if len(j) > 0]
                    
                    if len(ove_clusters_pair_id_w_size_set) < 1:
                        if_merge = False
                    else:
                        ## Collect all the id of clusters that need to merge
                        ove_clusters_id_set = list(set( [tuple(sorted( 
                            set(chain.from_iterable( [j[0] for j in ove_clusters_pair_id_w_size_set[k]] )) )) 
                                                         for k in range(len(ove_clusters_pair_id_w_size_set))] ))
                        
                        ## Remove subset
                        ove_clusters_id_set_sub = self.remove_subset(ove_clusters_id_set)
                        
                        ## Index of clusters that will be merged / not be merged
                        ove_clusters_id = list(set([k for j in ove_clusters_id_set_sub for k in j]))
                        non_ove_clusters_id = list(set(range(len(Cs)))^set(ove_clusters_id))
                        
                        Cs = [tuple(sorted( set(chain.from_iterable([Cs[j] for j in ove_clusters_id_set_sub[k]])) )) 
                              for k in range(len(ove_clusters_id_set_sub))] + [Cs[j] for j in non_ove_clusters_id]
                        
                        small_Cs = [j for j in Cs if (len(j) >= self.min_complex_size)&(len(j) <= self.max_complex_size)]
                        large_Cs = [j for j in Cs if len(j) > self.max_complex_size]
                        
                        if len(large_Cs) > 0 and self.merge_num <= self.max_merge_num:
                            keep_idx = np.array(list(set(chain.from_iterable([j for j in large_Cs]))))
                            small_Cs = small_Cs + self.extract_clusters(A, keep_idx=keep_idx)
                            small_Cs = self.remove_subset(small_Cs)
                        else:
                            if_merge = False
                        
                        Cs = small_Cs
                        self.merge_num += 1
        
        self.merge_num = 0
        return Cs
    
    def run(self, e_idx, e_w=None, dim=None, diag=None):
        if dim is None:
            print('Dim must specified')
            return
        else:
            A = self.create_edge_adj(e_idx, e_w=e_w, dim=dim, diag=diag)
            Cs = self.extract_clusters(A)
            Cs = self.merge_clusters(Cs, A)
            Cs = self.remove_subset(Cs)
            Cs = list(set(Cs))
            return Cs


##### Preprocess gold standard
gs = list(gs.groupby('ComplexName')['Gene_name'].apply(list).values)

## filter gs by exp name
gs_sub = [[i for i in gs[j] if i in name_idx_dict.keys()] for j in range(len(gs))]

## match name by idx
complexes = [tuple([name_idx_dict[i] for i in gs_sub[j]]) for j in range(len(gs_sub))]

## filter by complex size
complexes = [i for i in complexes if len(i) > 2]


#### Clustering
print('Start clustering process')
mcl = MCL_Cluster(merge_threshold=0.25, min_complex_size=3, max_complex_size=100, split_depth=3)
clusters = mcl.run(all_edge, e_w=all_pred, dim=epf.shape[0], diag=1.0)
overlap_score = mcl.get_ove_score(clusters, complexes)

best_matched_cluster = overlap_score.argmax(0) # best matched cluster for each complex
best_matched_complex = overlap_score.argmax(1) # best matched complex for each cluster
best_matched_cluster_score = overlap_score.max(0) # overlap perc between best matched cluster and each complex
best_matched_complex_score = overlap_score.max(1) # overlap perc between best matched complex and each cluster

best_matched_cluster_info = [tuple([i, best_matched_cluster[i], best_matched_cluster_score[i]]) for i in range(len(best_matched_cluster))]
# 0:complex id, 1:best cluster id, 2:overlap perc
best_matched_complex_info = [tuple([i, best_matched_complex[i], best_matched_complex_score[i]]) for i in range(len(best_matched_complex))]
# 0:cluster id, 1:best complex id, 2:overlap perc
best_matched_cluster_info = sorted(best_matched_cluster_info, key=lambda x: x[2], reverse=True)
best_matched_complex_info = sorted(best_matched_complex_info, key=lambda x: x[2], reverse=True)
cluster_cluster_score = mcl.get_ove_score(clusters, clusters)

print('Number of generated cluster: ' + str(len(clusters)))
print('Minimum cluster size: ' +  str(min([len(i) for i in clusters])), ', maximum cluster size: ' +  str(max([len(i) for i in clusters])))


#### Save output
cluster_path = args.output_path + '/Cluster'
if not os.path.exists(cluster_path):
      os.makedirs(cluster_path)


#### Save clusters
with open(cluster_path + '/clusters.txt', 'w') as output:
    for line in clusters:
        s = " ".join(map(str, line))
        output.write(s+'\n')

np.save(cluster_path + '/overlap_score.npy',             overlap_score)
np.save(cluster_path + '/best_matched_cluster_info.npy', best_matched_cluster_info)
np.save(cluster_path + '/best_matched_complex_info.npy', best_matched_complex_info)
np.save(cluster_path + '/cluster_cluster_score.npy',     cluster_cluster_score)
print('Done clustering')
print(' ')



