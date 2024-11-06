library(AUC)
library(PRROC)
library(ggplot2)
library(reshape2)
library(stringr)
library(cowplot)
library(reticulate)
library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)
library(plyr)
# https://cran.r-project.org/web/packages/PRROC/vignettes/PRROC.pdf
np <- import('numpy')

#===================================
#             Read data
#===================================
#### Set path
exp_name <- 'PXD002892'

exp_file <- list.files('/FREEPII_github/code/FREEPII/output')
exp_path <- file.path('/FREEPII_github/code/FREEPII/output', exp_file)
exp_cond <- unique(sapply(exp_file, function(x) strsplit(x, '_')[[1]][1], USE.NAMES = F, simplify = T))
exp_file
exp_path
exp_cond

split_dir <- paste('/FREEPII_github/Split-data', exp_name, sep = '/')
split_file <- sapply(split_dir, function(x) list.files(x)[grepl('ref', list.files(x))|grepl('heldout', list.files(x))], USE.NAMES = F, simplify = T)
split_path <- file.path(split_dir, split_file)
split_path



#===================================
#       Calculate performance
#===================================

performance <- data.frame()
for ( cur_exp_cond in exp_cond ) {
  print(cur_exp_cond)
  
  cur_exp_path <- exp_path[grepl(cur_exp_cond, exp_path)]
  cur_split_path <- split_path[grepl(cur_exp_cond, split_path)]
  
  cur_pred_df <- rbind(readRDS(cur_split_path[grepl('ref', cur_split_path)]), 
                       readRDS(cur_split_path[grepl('heldout', cur_split_path)]))
  
  cur_pred_df$Mean <- c(
    np$load(cur_exp_path[grepl('train',   cur_exp_path)])$f[['arr_0']],
    np$load(cur_exp_path[grepl('heldout', cur_exp_path)])$f[['arr_0']],
    np$load(cur_exp_path[grepl('delfold', cur_exp_path)])$f[['arr_0']]
  )
  
  idx <- grep('Mean', colnames(cur_pred_df))
  temp <- setNames(data.frame(matrix('TN', ncol = 1, nrow = nrow(cur_pred_df))), paste0('Label_Mean'))
  temp[(cur_pred_df$Label == 1)&(cur_pred_df[, idx]  > 0.5)&(complete.cases(cur_pred_df[, idx, drop = F])), 1] <- 'TP'
  temp[(cur_pred_df$Label == 1)&(cur_pred_df[, idx] <= 0.5)&(complete.cases(cur_pred_df[, idx, drop = F])), 1] <- 'FN'
  temp[(cur_pred_df$Label == 0)&(cur_pred_df[, idx]  > 0.5)&(complete.cases(cur_pred_df[, idx, drop = F])), 1] <- 'FP'
  temp[!complete.cases(cur_pred_df[, idx, drop = F]), 1] <- NA
  cur_pred_df <- cbind(cur_pred_df, temp)
  
  rm(idx)
  rm(temp)
  
  ## Separate by split
  for ( split_type in unique(cur_pred_df$Type) ){
    cur_temp <- cur_pred_df[cur_pred_df$Type==split_type, ]
    
    if( split_type=='Held_out' ) {
      cur_temp <- rbind(cur_temp[cur_temp$Label==1, ],cur_temp[sample(which(cur_temp$Label==0), sum(cur_temp$Label==1), replace = F), ])
    }
    
    col_idx <- which(grepl('cv|Mean', colnames(cur_pred_df))&grepl('Label', colnames(cur_pred_df)))
    
    for ( idx in col_idx ) {
      cur_cv_name <- gsub('Label_', '', colnames(cur_temp)[idx])
      
      cur_TP <- sum((cur_temp[,idx, drop=F]=='TP')&(cur_temp$Label==1), na.rm = T) / sum(cur_temp[complete.cases(cur_temp[,idx, drop=F]), ]$Label==1)
      cur_TN <- sum((cur_temp[,idx, drop=F]=='TN')&(cur_temp$Label==0), na.rm = T) / sum(cur_temp[complete.cases(cur_temp[,idx, drop=F]), ]$Label==0)
      cur_FN <- 1.0 - cur_TP
      cur_FP <- 1.0 - cur_TN
      cur_acc <- (cur_TP + cur_TN)/(cur_TP + cur_TN + cur_FP + cur_FN)
      cur_sn <- cur_TP/(cur_TP + cur_FN)
      cur_sp <- cur_TN/(cur_TN + cur_FP)
      cur_pr <- cur_TP/(cur_TP + cur_FP)
      cur_f1 <- 2*((cur_pr * cur_sn)/(cur_pr + cur_sn))
      cur_mcc <- ((cur_TP * cur_TN) - (cur_FP * cur_FN))/(sqrt(cur_TP + cur_FP) * sqrt(cur_TP + cur_FN) * sqrt(cur_TN + cur_FP) * sqrt(cur_TN + cur_FN))
      
      fg <- cur_temp[(cur_temp$Label==1)&complete.cases(cur_temp[,idx, drop=F]), 
                     which(grepl(cur_cv_name, colnames(cur_temp))&(!grepl('Label', colnames(cur_temp))))]
      bg <- cur_temp[(cur_temp$Label==0)&complete.cases(cur_temp[,idx, drop=F]), 
                     which(grepl(cur_cv_name, colnames(cur_temp))&(!grepl('Label', colnames(cur_temp))))]
      
      cur_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      cur_prc <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      
      cur_temp_ <- data.frame(
        ExpName = exp_name,
        Condition = cur_exp_cond,
        Split = str_to_title(split_type),
        CV = cur_cv_name,
        TP = cur_TP,
        TN = cur_TN,
        ACC = cur_acc,
        SN = cur_sn,
        SP = cur_sp,
        PR = cur_pr,
        MCC = cur_mcc,
        AUC_ROC = as.numeric(cur_roc$auc),
        AUC_PR_Integral = as.numeric(cur_prc$auc.integral),
        AUC_PR_Davis_Goadrich = as.numeric(cur_prc$auc.davis.goadrich)
      )
      performance <- rbind(performance, cur_temp_)
    }
    rm(cur_temp)
    rm(cur_TP)
    rm(cur_TN)
    rm(cur_FP)
    rm(cur_FN)
    rm(cur_acc)
    rm(cur_sn)
    rm(cur_sp)
    rm(cur_pr)
    rm(cur_f1)
    rm(cur_mcc)
    rm(bg)
    rm(fg)
    rm(cur_roc)
    rm(cur_prc)
    rm(cur_temp_)
  }
}






