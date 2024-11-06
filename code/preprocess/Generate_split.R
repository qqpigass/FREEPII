library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)
library(PrInCE)
library(plyr)

#==============================================
#      Read data and set split conditions
#==============================================
seed_idx <- 123
quant_type <- 'iBAQ'


#### Set path
exp_name <- list.files('/FREEPII_github/EPF-data')[!grepl('.R', list.files('/FREEPII_github/EPF-data'))]
exp_dir <- file.path('/FREEPII_github/EPF-data', exp_name)
exp_dir <- unlist(sapply(exp_dir, function(x) paste(x, list.files(x), sep = '/'), USE.NAMES = F))
exp_path <- sapply(exp_dir, function(x) paste(x, list.files(x)[grepl('.rds', list.files(x))], sep = '/'), USE.NAMES = F)
exp_path 

complex_path <- '/FREEPII_github/Protein complex-data/Complexes_gene_Human_filter.rds'
complex_path


#### Read complex
complex_df <- readRDS(complex_path)
nrow(complex_df)
colnames(complex_df)

complex_df <- unique(complex_df[,c('ComplexName', 'Gene_name')])
Complex_gnames <- sort(unique(complex_df$Gene_name))
c( nrow(complex_df), length(unique(complex_df$ComplexName)), length(Complex_gnames) )


#### Set conditions
train_split_ratio <- 7      # train:test = 70:30
train_split_type <- 'pair'  # if single, split by proteins; if pair, split by PPIs
cv_fold <- 5
cv_split_type <- 'pair'
pn_ratio <- 1               # positive:negative = 1:1
c( paste0('train_split_ratio: ', train_split_ratio),
   paste0('train_split_type: ', train_split_type),
   paste0('cv_fold: ', cv_fold),
   paste0('cv_split_type: ', cv_split_type),
   paste0('pn_ratio: ', pn_ratio) )



#=================================
#         Iter over epf
#=================================
out_path <- '/FREEPII_github/Split-data'
if ( !file.exists(out_path) ) {
  dir.create(out_path)
}

grep_pattern <- paste('heavy', 'tr', 'Cn', sep = '|')
grep_pattern

for ( path_idx in c(1:length(exp_path)) ) {
  path <- exp_path[path_idx]
  
  cur_exp_name <- tail(head(strsplit(path, '/')[[1]], -2), 1)
  cur_exp_cond <- gsub('.rds', '', basename(path))
  print(paste(cur_exp_name, cur_exp_cond, ':', sep = ' '))
  
  cur_mat = as.matrix(readRDS(path))
  all_ppi <- setNames(data.frame(t(combn(rownames(cur_mat), 2))), c('Gene_name_A', 'Gene_name_B'))
  
  cur_com_df <- complex_df[complex_df$Gene_name %in% rownames(cur_mat), ]
  
  
  #### Turn dataframe into list and count complex subunit
  cur_com_list <- setNames(cur_com_df$Gene_name, cur_com_df$ComplexName)
  cur_com_list <- lapply(split(cur_com_list, names(cur_com_list)), unname)
  cur_com_N <- lengths(cur_com_list)
  rm(cur_com_df)
  
  
  #### Extract pairwise interactions in complex
  cur_com_ppi <- cur_com_list %>%
    map_dfr(~ tidyr::crossing(Gname_A = ., Gname_B = .), .id = 'ComplexName') %>%
    dplyr::filter(Gname_A != Gname_B) %>% # remove self interaction
    dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>% # sort ppi
    dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
    dplyr::select(ComplexName, Gene_name_A, Gene_name_B) %>%
    dplyr::distinct() %>% # remove flip
    dplyr::mutate(Gene_name_N = cur_com_N[ComplexName] %>% unname()) %>%
    dplyr::group_by(Gene_name_A, Gene_name_B) %>%
    dplyr::arrange(Gene_name_N) %>%
    dplyr::mutate(keep = row_number() == 1) %>% # for interactions found in more than one complex, pick the smaller one
    dplyr::ungroup() %>%
    dplyr::filter(keep)
  rm(cur_com_list)
  rm(cur_com_N)
  
  
  #### Remove complex PPI with no co-peak in epf
  cur_com_ppi_copeak <- cur_mat[rownames(cur_mat) %in% unique(c(cur_com_ppi$Gene_name_A, cur_com_ppi$Gene_name_B)), ]
  cur_com_ppi_copeak <- crossprod( t(cur_com_ppi_copeak > 0.01), t(cur_com_ppi_copeak > 0.01) )
  diag(cur_com_ppi_copeak) <- 0
  cur_com_ppi_copeak[lower.tri(cur_com_ppi_copeak)] <- 0L
  
  cur_mask <- matrix(0, nrow(cur_com_ppi_copeak), ncol(cur_com_ppi_copeak))
  rownames(cur_mask) = colnames(cur_mask) = rownames(cur_com_ppi_copeak)
  cur_mask[as.matrix(cur_com_ppi[,c('Gene_name_A', 'Gene_name_B')])] <- 1
  cur_com_ppi_copeak <- cur_com_ppi_copeak * cur_mask
  
  cur_com_ppi_copeak <- as.data.frame(as.table(cur_com_ppi_copeak)) %>%
    dplyr::filter(Var1!=Var2) %>%
    dplyr::filter(Freq!=0) %>%
    dplyr::mutate(Gname_A = as.character(Var1)) %>%
    dplyr::mutate(Gname_B = as.character(Var2)) %>%
    dplyr::select(Gname_A, Gname_B, Freq) %>%
    dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>%
    dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
    dplyr::select(Gene_name_A, Gene_name_B) %>%
    dplyr::distinct()
  
  cur_com_ppi <- dplyr::inner_join(cur_com_ppi, cur_com_ppi_copeak, by = c('Gene_name_A', 'Gene_name_B'))
  rm(cur_com_ppi_copeak)
  rm(cur_mask)
  
  
  #### Keep complex with size > 2
  cur_com_ppi$Gene_name_N <- sapply( 
    c(1:nrow(cur_com_ppi)), function (x) length(unique(as.character(as.matrix(cur_com_ppi[cur_com_ppi$ComplexName==cur_com_ppi$ComplexName[x], 2:3]))))
  )
  cur_com_ppi <- cur_com_ppi[cur_com_ppi$Gene_name_N > 2, ]
  
  cat(paste(
    paste0('Unique complex number: ', length(unique(cur_com_ppi$ComplexName)), 
           '; Unique complex genes: ', length(unique(c(cur_com_ppi$Gene_name_A, cur_com_ppi$Gene_name_B)))),
    paste0('All complex PPI: ', nrow(cur_com_ppi), 
           '; All unique complex PPI: ', nrow((unique( cur_com_ppi[,c('Gene_name_A', 'Gene_name_B')] )))),
    strrep('=', 60),
    sep = '\n'
  ), '\n')
  
  temp <- cur_com_ppi %>%
    dplyr::select(Gene_name_A, Gene_name_B) %>%
    dplyr::rename(Gname_A = Gene_name_A, Gname_B = Gene_name_B) %>%
    dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>%
    dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
    dplyr::select(Gene_name_A, Gene_name_B) %>%
    dplyr::filter(Gene_name_A!=Gene_name_B) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Label = 1)
  
  com_gnames <- unique(c(temp$Gene_name_A, temp$Gene_name_B))
  
  if ( grepl(grep_pattern, cur_exp_cond) ) {
    #### Split the genes as train set and test set
    if ( train_split_type=='single' ) {
      set.seed(seed_idx)
      com_gnames <- sample(com_gnames) # shuffle the genes
      
      ## Get the split of genes
      border <- floor(length(com_gnames) * train_split_ratio*0.1)
      train_gnames <- com_gnames[c(1 : border)]
      test_gnames <- setdiff(com_gnames, train_gnames)
      
      ## All pairwise interaction in each split
      train_ppi <- setNames(as.data.frame(t(combn(train_gnames, 2))), c('Gname_A', 'Gname_B')) %>%
        dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>%
        dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
        dplyr::select(Gene_name_A, Gene_name_B) %>%
        dplyr::filter(Gene_name_A!=Gene_name_B) %>%
        dplyr::distinct()
      
      test_ppi <- setNames(as.data.frame(t(combn(test_gnames, 2))), c('Gname_A', 'Gname_B')) %>%
        dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>%
        dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
        dplyr::select(Gene_name_A, Gene_name_B) %>%
        dplyr::filter(Gene_name_A!=Gene_name_B) %>%
        dplyr::distinct()
      
      ## Label interactions
      reference <- merge(train_ppi, temp, by = c('Gene_name_A', 'Gene_name_B'), all.x = T)
      held_out <- merge(test_ppi, temp, by = c('Gene_name_A', 'Gene_name_B'), all.x = T)
      
      reference[!complete.cases(reference), ]$Label <- 0
      held_out[!complete.cases(held_out), ]$Label <- 0
      
      temp_ <- rbind(reference, held_out)
      non_com_ppi <- merge(all_ppi, temp_, by=c('Gene_name_A', 'Gene_name_B'), all.x = T)
      non_com_ppi$Type <- 'Exp'
      non_com_ppi <- non_com_ppi[!complete.cases(non_com_ppi$Label), ]
      
      cat(paste(
        'split by single',
        paste0('Train genes: ', length(train_gnames), '; Test genes: ', length(test_gnames)),
        paste0('Unique reference genes: ', length(unique(c(reference$Gene_name_A, reference$Gene_name_B))), 
               '; Unique held-out genes: ', length(unique(c(held_out$Gene_name_A, held_out$Gene_name_B))),
               '; Unique exp genes: ', length(unique(c(non_com_ppi$Gene_name_A, non_com_ppi$Gene_name_B))),
               '; Unique all genes: ', nrow(cur_mat)
        ),
        paste0('TP in reference: ', sum(reference$Label==1), '; TP in held-out: ', sum(held_out$Label==1)),
        paste0('TN in reference: ', sum(reference$Label==0), '; TN in held-out: ', sum(held_out$Label==0)),
        paste0('All complex PPIs: ', nrow(temp_), ' All reference PPIs: ', nrow(reference), '; All held-out PPIs: ', nrow(held_out)),
        paste0('All exp PPIs: ', nrow(non_com_ppi)),
        strrep('=', 60),
        sep = '\n'), '\n')
      
      rm(train_ppi)
      rm(test_ppi)
      rm(train_gnames)
      rm(test_gnames)
      rm(border)
      rm(temp_)
      
    } else if ( train_split_type=='pair' ) {
      ## Get all complex interactions
      all_com_ppi <- setNames(as.data.frame(t(combn(com_gnames, 2))), c('Gname_A', 'Gname_B')) %>%
        dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>%
        dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
        dplyr::select(Gene_name_A, Gene_name_B) %>%
        dplyr::filter(Gene_name_A!=Gene_name_B) %>%
        dplyr::distinct()
      
      set.seed(seed_idx)
      shuff = sample(seq_len(nrow(all_com_ppi)))
      all_com_ppi %<>% magrittr::extract(shuff, ) # shuffle the interactions
      
      ## Get the split of interactions
      border = floor(nrow(all_com_ppi) * train_split_ratio*0.1)
      train_ppi = all_com_ppi[c(1 : border), ]
      test_ppi = all_com_ppi[c( (border + 1) : nrow(all_com_ppi)), ]
      
      ## Label interactions
      reference <- merge(train_ppi, temp, by = c('Gene_name_A', 'Gene_name_B'), all.x = T)
      held_out <- merge(test_ppi, temp, by = c('Gene_name_A', 'Gene_name_B'), all.x = T)
      
      reference[!complete.cases(reference), ]$Label <- 0
      held_out[!complete.cases(held_out), ]$Label <- 0
      
      temp_ <- rbind(reference, held_out)
      non_com_ppi <- merge(all_ppi, temp_, by=c('Gene_name_A', 'Gene_name_B'), all.x = T)
      non_com_ppi$Type <- 'Exp'
      non_com_ppi <- non_com_ppi[!complete.cases(non_com_ppi$Label), ]
      
      cat(paste(
        'split by pair',
        paste0('Unique train genes: ', length(unique(c(train_ppi$Gene_name_A, train_ppi$Gene_name_B))), 
               '; Unique test genes: ', length(unique(c(test_ppi$Gene_name_A, test_ppi$Gene_name_B)))),
        paste0('Unique reference genes: ', length(unique(c(reference$Gene_name_A, reference$Gene_name_B))), 
               '; Unique held-out genes: ', length(unique(c(held_out$Gene_name_A, held_out$Gene_name_B))),
               '; Unique exp genes: ', length(unique(c(non_com_ppi$Gene_name_A, non_com_ppi$Gene_name_B))),
               '; Unique all genes: ', nrow(cur_mat)
        ),
        paste0('TP in reference: ', sum(reference$Label==1), '; TP in held-out: ', sum(held_out$Label==1)),
        paste0('TN in reference: ', sum(reference$Label==0), '; TN in held-out: ', sum(held_out$Label==0)),
        paste0('All complex PPIs: ', nrow(all_com_ppi), ' All reference PPIs: ', nrow(reference), '; All held-out PPIs: ', nrow(held_out)),
        paste0('All exp PPIs: ', nrow(non_com_ppi)),
        strrep('=', 60),
        sep = '\n'), '\n')
      
      rm(train_ppi)
      rm(test_ppi)
      rm(shuff)
      rm(border)
      rm(all_com_ppi)
      rm(temp_)
    }
    
    reference$Type <- 'Ref'
    
    
    #### Split training set into fold for cross validation
    if (cv_split_type == 'single') {
      train_gnames <- unique(c(reference$Gene_name_A, reference$Gene_name_B)) # split train genes into folds_keep
      
      borders = seq(1, length(train_gnames), length(train_gnames) / cv_fold) %>% floor()
      
      set.seed(seed_idx)
      fold_idxs = seq_len(length(train_gnames)) %>% sample() %>% split(borders)
      folds_keep = purrr::map(fold_idxs, ~ magrittr::extract(train_gnames, .)) %>% setNames(seq_along(.))
      folds_keep <- folds_keep %>% purrr::map(~ reference[(reference$Gene_name_A %in% .)&(reference$Gene_name_B %in% .), ]) # convert back to long
      
      cat(paste(
        'split by single',
        paste0('PPIs in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) nrow(folds_keep[[x]])))),
        paste0('Unique genes in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) length(unique(c(folds_keep[[x]]$Gene_name_A, folds_keep[[x]]$Gene_name_B)))))),
        paste0('TP in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) sum(folds_keep[[x]]$Label==1)))),
        paste0('TN in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) sum(folds_keep[[x]]$Label==0)))),
        strrep('=', 60),
        sep = '\n'
      ), '\n')
      
      rm(borders)
      rm(fold_idxs)
      
    } else if (cv_split_type == 'pair') {
      borders = seq(1, nrow(reference), nrow(reference) / cv_fold) %>% floor()
      
      set.seed(seed_idx)
      fold_idxs = seq_len(nrow(reference)) %>% sample() %>% split(borders)
      folds_keep = purrr::map(fold_idxs, ~ magrittr::extract(reference, ., )) %>% setNames(seq_along(.))
      
      cat(paste(
        'split by pair',
        paste0('PPIs in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) nrow(folds_keep[[x]])))),
        paste0('Unique genes in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) length(unique(c(folds_keep[[x]]$Gene_name_A, folds_keep[[x]]$Gene_name_B)))))),
        paste0('TP in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) sum(folds_keep[[x]]$Label==1)))),
        paste0('TN in each fold: ', toString(sapply(c(1:length(folds_keep)), function(x) sum(folds_keep[[x]]$Label==0)))),
        strrep('=', 60),
        sep = '\n'
      ), '\n')
      
      rm(borders)
      rm(fold_idxs)
    }
    
    
    #### Fix ratio for PPI in split
    cur_del_fold_TN <- data.frame()
    for ( fold_idx in c(1:length(folds_keep)) ){
      temp_ <- folds_keep[[fold_idx]]
      p_df <- temp_[temp_$Label==1, ]
      n_df <- temp_[temp_$Label==0, ]
      n_idx <- c(1:nrow(n_df))
      
      set.seed(seed_idx)
      keep_idx <- sample(n_idx, nrow(p_df)*pn_ratio)
      del_idx <- n_idx[!n_idx %in% keep_idx]
      
      folds_keep[[fold_idx]] <- sample_frac(rbind(p_df, n_df[keep_idx, ]), 1)
      folds_keep[[fold_idx]]$Type <- 'Train'
      cur_del_fold_TN <- rbind(cur_del_fold_TN, n_df[del_idx, ])
      
      rm(p_df)
      rm(n_df)
      rm(n_idx)
      rm(keep_idx)
      rm(del_idx)
      rm(temp_)
    }
    cat(paste(
      paste0('PPIs in each balance fold: ', toString(sapply(c(1:length(folds_keep)), function(x) nrow(folds_keep[[x]])))),
      paste0('Unique genes in each balance fold: ', toString(sapply(c(1:length(folds_keep)), function(x) length(unique(c(folds_keep[[x]]$Gene_name_A, folds_keep[[x]]$Gene_name_B)))))),
      paste0('TP in each balance fold: ', toString(sapply(c(1:length(folds_keep)), function(x) sum(folds_keep[[x]]$Label==1)))),
      paste0('TN in each balance fold: ', toString(sapply(c(1:length(folds_keep)), function(x) sum(folds_keep[[x]]$Label==0)))),
      paste0('Removed negative PPIs: ', nrow(cur_del_fold_TN)),
      strrep('=', 60),
      sep = '\n'
    ), '\n')
    
    held_out$Type <- 'Held_out'
    cur_del_fold_TN$Type <- 'Del_fold'
    
    reference <- merge(reference, cur_del_fold_TN, by=c('Gene_name_A', 'Gene_name_B'), all.x = T)
    reference <- reference[!complete.cases(reference$Type.y), ]
    reference <- reference[,1:4]
    colnames(reference) <- gsub('.x$', '', colnames(reference))
    
    held_out <- rbind(held_out, cur_del_fold_TN)
    rm(cur_del_fold_TN)
    
    
    #### Create cross splits
    split_cv = list()
    for (test_idx in seq_along(folds_keep)) {
      test = folds_keep[[test_idx]]
      train = folds_keep[-test_idx] %>% bind_rows()
      test$Type <- 'Test'
      split_cv[[test_idx]] = list(train = train, test = test)
      
      rm(test)
      rm(train)
    }
    split_cv %<>% setNames(seq_along(.))
    rm(folds_keep)
    
    cat(paste(
      paste0('Train PPIs in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) nrow(split_cv[[x]]$train)))),
      paste0('Test PPIs in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) nrow(split_cv[[x]]$test)))),
      paste0('Train genes in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) length(unique(c(split_cv[[x]]$train$Gene_name_A, split_cv[[x]]$train$Gene_name_B)))))),
      paste0('Test genes in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) length(unique(c(split_cv[[x]]$test$Gene_name_A, split_cv[[x]]$test$Gene_name_B)))))),
      paste0('Train TP in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) sum(split_cv[[x]]$train$Label==1)))),
      paste0('Test TP in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) sum(split_cv[[x]]$test$Label==1)))),
      paste0('Train TN in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) sum(split_cv[[x]]$train$Label==0)))),
      paste0('Test TN in each balance fold: ', toString(sapply(c(1:length(split_cv)), function(x) sum(split_cv[[x]]$test$Label==0)))),
      strrep('=', 100),
      sep = '\n'
    ), '\n')
    
    if ( !file.exists( file.path(out_path, cur_exp_name) ) ) {
      dir.create( file.path(out_path, cur_exp_name) )
    }
    
    saveRDS(cur_com_ppi, file.path(out_path, cur_exp_name, paste0('complex_ppi_', cur_exp_cond, '.rds')))
    saveRDS(reference, file.path(out_path, cur_exp_name, paste0('ref_ppi_', train_split_type, '_', cur_exp_cond, '.rds')))
    saveRDS(held_out, file.path(out_path, cur_exp_name, paste0('heldout_delfold_ppi_', train_split_type, '_', cur_exp_cond, '.rds')))
    saveRDS(non_com_ppi, file.path(out_path, cur_exp_name, paste0('exp_ppi_', train_split_type, '_', cur_exp_cond, '.rds')))
    saveRDS(split_cv, file.path(out_path, cur_exp_name, paste0('cv_split_', cv_split_type, '_', cur_exp_cond, '.rds')))
    
    rm(reference)
    rm(held_out)
    rm(non_com_ppi)
    rm(split_cv)
  } else {
    all_com_ppi <- setNames(as.data.frame(t(combn(com_gnames, 2))), c('Gname_A', 'Gname_B')) %>%
      dplyr::mutate(Gene_name_A = ifelse(Gname_A < Gname_B, Gname_A, Gname_B)) %>%
      dplyr::mutate(Gene_name_B = ifelse(Gname_B < Gname_A, Gname_A, Gname_B)) %>%
      dplyr::select(Gene_name_A, Gene_name_B) %>%
      dplyr::filter(Gene_name_A!=Gene_name_B) %>%
      dplyr::distinct()
    
    temp_ <- merge(all_com_ppi, temp, by=c('Gene_name_A', 'Gene_name_B'), all.x = T)
    temp_[!complete.cases(temp_$Label), ]$Label <- 0
    
    temp_ <- merge(all_ppi, temp_, by=c('Gene_name_A', 'Gene_name_B'), all.x = T)
    temp_$Type <- NA
    temp_[!complete.cases(temp_$Label), ]$Type <- 'Exp'
    
    cat(paste(
      paste0('TP in exp data: ', sum(temp_$Label==1, na.rm = T)),
      paste0('TN in exp data: ', sum(temp_$Label==0, na.rm = T)),
      paste0('Remain PPI in exp data: ', sum(!complete.cases(temp_$Label))),
      strrep('=', 100),
      sep = '\n'
    ), '\n')
    
    if ( !file.exists( file.path(out_path, cur_exp_name) ) ) {
      dir.create( file.path(out_path, cur_exp_name) )
    }
    
    saveRDS(cur_com_ppi, file.path(out_path, cur_exp_name, paste0('complex_ppi_', cur_exp_cond, '.rds')))
    saveRDS(temp_, file.path(out_path, cur_exp_name, paste0('all_ppi_', cur_exp_cond, '.rds')))
    
    rm(temp_)
    rm(all_com_ppi)
  }
  
  rm(cur_mat)
  rm(cur_com_ppi)
  rm(all_ppi)
  rm(temp)
}







