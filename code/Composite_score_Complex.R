library(reshape2)
library(ggplot2)
library(cowplot)
library(stringr)
library(RColorBrewer)
library(reticulate)
pd <- import('pandas')

#===================================
#             Read data
#===================================
#### Set path
exp_name <- 'PXD002892'

exp_file <- list.files('/FREEPII_github/code/analysis/Cluster/')
exp_file <- exp_file[grepl('clusters.txt', exp_file)]
exp_path <- file.path('/FREEPII_github/code/analysis/Cluster/', exp_file)
exp_cond <- unique(sapply(exp_file, function(x) strsplit(x, '_')[[1]][1], USE.NAMES = F, simplify = T))
exp_file
exp_path
exp_cond


#### Read GS data
gs_name_df <- readRDS('/FREEPII_github/Protein complex-data/Complexes_gene_Human_filter.rds')
dim(gs_name_df)
head(gs_name_df)

gs_name_df_flat <- data.frame(matrix(NA, ncol = max(gs_name_df$Gene_name_N)+2, nrow = length(unique(gs_name_df$ComplexName))))
colnames(gs_name_df_flat) <- c('ComplexName', 'ComplexSize', paste0('Comp_', c(1:max(gs_name_df$Gene_name_N))))
gs_name_df_flat$ComplexName <- unique(gs_name_df$ComplexName)
for ( idx in c(1:nrow(gs_name_df_flat)) ) {
  cur_gene <- gs_name_df[gs_name_df$ComplexName==gs_name_df_flat$ComplexName[idx], ]$Gene_name
  gs_name_df_flat$ComplexSize[idx] <- length(cur_gene)
  gs_name_df_flat[idx, 3:(length(cur_gene)+2)] <- cur_gene
}
dim(gs_name_df_flat)
head(gs_name_df_flat)



#=============================================
#       Calculate Composite score
#=============================================
load_dir <- '/FREEPII_github/input/PXD002892/'
load_dir

performance <- data.frame()
for ( cur_exp_cond in exp_cond ) {
  cur_exp_path <- exp_path[grepl(cur_exp_cond, exp_path)]
  
  cur_name_idx_dict <- pd$read_pickle(file.path(load_dir, paste0('name_idx_dict_', cur_exp_cond, '.pickle')))
  cur_name_idx_dict <- data.frame(name = names(cur_name_idx_dict), idx = as.numeric(cur_name_idx_dict))
  
  cur_cluster_list <- read.table(cur_exp_path, header = F, sep = '\t')
  
  cur_cluster_idx_df_flat <- data.frame(
    matrix(NA, ncol = max(sapply(c(1:nrow(cur_cluster_list)), function(x) 
      length(strsplit(cur_cluster_list$V1[x], ' ')[[1]]))) + 2, nrow = nrow(cur_cluster_list)))
  
  colnames(cur_cluster_idx_df_flat) <- c('ComplexName', 'ComplexSize', paste0('Comp_', c(1:(ncol(cur_cluster_idx_df_flat)-2))))
  cur_cluster_idx_df_flat$ComplexName <- paste0('Complex_', c(1:nrow(cur_cluster_idx_df_flat)))
  for ( idx in c(1:nrow(cur_cluster_idx_df_flat)) ) {
    cur_gene <- as.numeric(strsplit(cur_cluster_list[idx, ], ' ')[[1]])
    cur_cluster_idx_df_flat$ComplexSize[idx] <- length(cur_gene)
    cur_cluster_idx_df_flat[idx, 3:(length(cur_gene)+2)] <- cur_gene
  }
  
  cur_gs_idx_df_flat <- gs_name_df_flat
  for ( idx in c(1:nrow(cur_gs_idx_df_flat)) ) {
    cur_gene <- cur_name_idx_dict[cur_name_idx_dict$name %in% as.character(cur_gs_idx_df_flat[idx, c(3:(cur_gs_idx_df_flat$ComplexSize[idx]+2))]), ]$idx
    if ( length(cur_gene) > 2 ) {
      cur_gs_idx_df_flat$ComplexSize[idx] <- length(cur_gene)
      cur_gs_idx_df_flat[idx, 3:ncol(cur_gs_idx_df_flat)] <- NA
      cur_gs_idx_df_flat[idx, 3:(length(cur_gene)+2)] <- cur_gene
    } else {
      cur_gs_idx_df_flat$ComplexSize[idx] <- 0
      cur_gs_idx_df_flat[idx, 3:ncol(cur_gs_idx_df_flat)] <- NA
    }
  }
  cur_gs_idx_df_flat <- cur_gs_idx_df_flat[cur_gs_idx_df_flat$ComplexSize > 0, ]
  
  
  #### Calculate composite score
  cur_ove_df <- data.frame(matrix(0, ncol = nrow(cur_gs_idx_df_flat), nrow = nrow(cur_cluster_idx_df_flat)))
  cur_ove_score_df <- data.frame(matrix(0, ncol = nrow(cur_gs_idx_df_flat), nrow = nrow(cur_cluster_idx_df_flat)))
  
  colnames(cur_ove_df) <- cur_gs_idx_df_flat$ComplexName
  colnames(cur_ove_score_df) <- cur_gs_idx_df_flat$ComplexName
  
  rownames(cur_ove_df) <- cur_cluster_idx_df_flat$ComplexName
  rownames(cur_ove_score_df) <- cur_cluster_idx_df_flat$ComplexName
  
  ## Check match for each Pred to each Ref
  for( pred in rownames(cur_ove_df) ){
    for( ref in colnames(cur_ove_df) ){
      pred_idx <- which(cur_cluster_idx_df_flat$ComplexName==pred)
      ref_idx <- which(cur_gs_idx_df_flat$ComplexName==ref)
      
      pred_comp <- cur_cluster_idx_df_flat[pred_idx, 3:(as.numeric(cur_cluster_idx_df_flat$ComplexSize[pred_idx])+2)]
      ref_comp <- cur_gs_idx_df_flat[ref_idx, 3:(as.numeric(cur_gs_idx_df_flat$ComplexSize[ref_idx])+2)]
      
      ove_num <- length(which(pred_comp %in% ref_comp))
      
      cur_ove_df[pred_idx, ref_idx] <- ove_num
      cur_ove_score_df[pred_idx, ref_idx] <- (ove_num*ove_num)/(length(pred_comp)*length(ref_comp))
    }
  }
  
  ## Calculate Overlapp score ( Fraction of pred with any overlap score >= 0.25 )
  Overlapp <- length(which(apply(cur_ove_score_df, 1, function(x) any(x >= 0.25)))) / nrow(cur_ove_score_df)
  
  ## Calculate Sensitivity score (Sn) ( sum( maximum overlap for each ref )/sum( size of each ref ) )
  max_idx <- as.numeric(apply(cur_ove_df, 2, which.max))
  Sensitivity <- sum(sapply( c(1:ncol(cur_ove_df)), function(x) as.numeric(cur_ove_df[max_idx[x], x]) ))/sum(as.numeric(cur_gs_idx_df_flat$ComplexSize))
  
  ## Calculate Positive predictive value score (PPV) ( sum( maximum overlap for each pred )/ sum( all overlap for each pred ) )
  max_idx <- as.numeric(apply(cur_ove_df, 1, which.max))
  PPV <- sum(sapply( c(1:nrow(cur_ove_df)), function(x) as.numeric(cur_ove_df[x, max_idx[x]]) ))/sum(cur_ove_df)
  
  ## Calculate Accuracy score (ACC) ( sqrt(Sn*PPV) )
  Accuracy <- sqrt(Sensitivity*PPV)
  
  ## Calculate Maximum matching ratio score (MMR) ( sum( maximum overlap score for each ref )/ number of ref )
  max_idx <- as.numeric(apply(cur_ove_score_df, 2, which.max))
  MMR <- sum(sapply( c(1:ncol(cur_ove_score_df)), function(x) as.numeric(cur_ove_score_df[max_idx[x], x]) )/ncol(cur_ove_score_df))
  
  
  #### Save Evaluate score
  cur_cluster_eval_df <- data.frame(
    ExpName = exp_name,
    Condition = cur_exp_cond,
    Metrics = c('Overlapp', 'Sensitivity', 'PPV', 'Accuracy', 'MMR'),
    Values = c(Overlapp, Sensitivity, PPV, Accuracy, MMR)
  )
  performance <- rbind(performance, cur_cluster_eval_df)
  
  print(cur_cluster_eval_df)
  print('=====================================================================')
}






