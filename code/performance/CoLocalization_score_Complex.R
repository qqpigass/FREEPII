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

exp_file <- list.files('/FREEPII_github/code/FREEPII/Cluster/')
exp_file <- exp_file[grepl('clusters.txt', exp_file)]
exp_path <- file.path('/FREEPII_github/code/FREEPII/Cluster/', exp_file)
exp_cond <- unique(sapply(exp_file, function(x) strsplit(x, '_')[[1]][1], USE.NAMES = F, simplify = T))
exp_file
exp_path
exp_cond

Uname_subloc_path <- '/FREEPII_github/UniProt-data/uniprot_subloc_Human.rds'
name_path <- '/FREEPII_github/UniProt-data/uniprot_gene_Human.rds'
co_subloc_path <- '/FREEPII_github/CoLocalization_weight/ppi_CoLocalization_Human.rds'
Uname_subloc_path
name_path
co_subloc_path


#### Read data
## Sub-cellular location
Uname_subloc_df <- readRDS(Uname_subloc_path)
nrow(Uname_subloc_df)
colnames(Uname_subloc_df)
c(length(unique(Uname_subloc_df$UniProt_entry)), length(unique(Uname_subloc_df$Location)))

## Gene name
name_df <- readRDS(name_path)
nrow(name_df)
colnames(name_df)

## Co-localization
co_subloc_df <- readRDS(co_subloc_path)
nrow(co_subloc_df)
colnames(co_subloc_df)
length(unique(c(co_subloc_df$Term, co_subloc_df$CoocTerm)))



#====================================================
#     Combined sub-cellular location to exp data
#====================================================
## Combined gene name with sub-cellular location information
Gname_subloc_df <- unique( merge(
  name_df[,c('UniProt_entry', 'Gene_name')], 
  Uname_subloc_df[,c('UniProt_entry', 'Location')], 
  by = 'UniProt_entry')[, c('Gene_name', 'Location')] )
c( nrow(Gname_subloc_df), sum(complete.cases(Gname_subloc_df)) )
colnames(Gname_subloc_df)



#=============================================
#       Calculate co-localization score
#=============================================
load_dir <- paste('/FREEPII_github/input', exp_name, sep = '/')
load_dir

performance <- data.frame()
for ( cur_exp_cond in exp_cond ) {
  cur_exp_path <- exp_path[grepl(cur_exp_cond, exp_path)]
  
  cur_name_idx_dict <- pd$read_pickle(file.path(load_dir, paste0('name_idx_dict_', cur_exp_cond, '.pickle')))
  cur_name_idx_dict <- data.frame(name = names(cur_name_idx_dict), idx = as.numeric(cur_name_idx_dict))
  
  cur_cluster_list <- read.table(cur_exp_path, header = F, sep = '\t')
  
  cur_cluster_df <- data.frame(matrix(0, ncol = 1, nrow = 0))
  for ( idx in c(1:nrow(cur_cluster_list)) ) {
    cur_name <- cur_name_idx_dict$name[cur_name_idx_dict$idx %in% as.numeric(strsplit(cur_cluster_list[idx, ], ' ')[[1]])]
    cur_df <- unique(data.frame( ComplexName = paste0('Complex_', idx), ComplexSize = length(cur_name), Gene_name = cur_name ))
    cur_cluster_df <- rbind(cur_cluster_df, cur_df)
  }
  rm(cur_name)
  rm(cur_df)
  
  cur_Gname_subloc_df <- unique( Gname_subloc_df[Gname_subloc_df$Gene_name %in% cur_name_idx_dict$name,] )
  
  cur_subloc_subloc_list <- list()
  cur_subloc_Gname_list <- list()
  for ( cur_subloc in unique(cur_Gname_subloc_df$Location) ){
    ## Co-localization list
    cur_co_subloc <- as.character(unique(c( co_subloc_df$CoocTerm[co_subloc_df$Term==cur_subloc], co_subloc_df$Term[co_subloc_df$CoocTerm==cur_subloc] )))
    if ( length(cur_co_subloc) > 0 ) {
      cur_subloc_subloc_list[[cur_subloc]] <- cur_co_subloc
    }
    
    ## Co-localization group list
    cur_subloc_Gname_list[[cur_subloc]] <- as.character( unique(cur_Gname_subloc_df$Gene_name[cur_Gname_subloc_df$Location==cur_subloc]) )
  }
  rm(cur_subloc)
  rm(cur_co_subloc)
  
  cur_df <- data.frame(matrix(0, ncol = 1, nrow = 0))
  for ( cur_com in unique(cur_cluster_df$ComplexName) ) {
    cur_gname <- unique( cur_cluster_df[cur_cluster_df$ComplexName==cur_com, 'Gene_name'] )
    
    cur_dff <- data.frame(matrix(0, ncol = 1, nrow = 0))
    cur_max_same_subloc <- 0
    cur_max_diff_subloc1 <- 0
    cur_max_diff_subloc2 <- 0
    cur_max_same_group_size <- 0
    cur_max_diff_group_size <- 0
    cur_same_unassign_group <- cur_gname
    cur_diff_unassign_group <- cur_gname
    
    for ( idx in c(1:length(cur_subloc_Gname_list)) ) {
      cur_same_subloc <- names(cur_subloc_Gname_list)[idx]
      cur_same_group <- cur_gname[cur_gname %in% cur_subloc_Gname_list[[cur_same_subloc]]]
      cur_same_group_size <- length(cur_same_group)
      
      if ( (cur_same_group_size > cur_max_same_group_size)&(cur_same_group_size >= 2) ) {
        cur_max_same_subloc <- cur_same_subloc
        cur_max_same_group_size <- cur_same_group_size
      }
      if ( length(cur_same_unassign_group) > 0 ) {
        cur_same_unassign_group <- setdiff(cur_same_unassign_group, cur_same_group)
      }
      
      if ( cur_same_subloc %in% names(cur_subloc_subloc_list) ) {
        for ( cur_diff_subloc in cur_subloc_subloc_list[[cur_same_subloc]] ) {
          cur_diff_group <- cur_gname[cur_gname %in% cur_subloc_Gname_list[[cur_diff_subloc]]]
          cur_diff_group <- union(cur_same_group, cur_diff_group)
          cur_diff_group_size <- length(cur_diff_group)
          
          if ( (cur_diff_group_size > cur_max_diff_group_size)&(cur_diff_group_size >= 2) ) {
            cur_max_diff_subloc1 <- cur_same_subloc
            cur_max_diff_subloc2 <- cur_diff_subloc
            cur_max_diff_group_size <- cur_diff_group_size
          }
          if ( length(cur_diff_unassign_group) > 0 ) {
            cur_diff_unassign_group <- setdiff(cur_diff_unassign_group, cur_diff_group)
          }
        }
      }
    }
    
    cur_df <- rbind(cur_df, data.frame(
      ComplexName = cur_com,
      ComplexSize = length(cur_gname),
      Main_subloc = cur_max_same_subloc,
      Pair_subloc = cur_max_same_subloc,
      CoSublocType = 'Same',
      MaxGroupSize = cur_max_same_group_size,
      AssignedGroupSize = length(cur_gname)-length(cur_same_unassign_group),
      Score = cur_max_same_group_size / ( length(cur_gname)-length(cur_same_unassign_group) )
    ))
    
    cur_df <- rbind(cur_df, data.frame(
      ComplexName = cur_com,
      ComplexSize = length(cur_gname),
      Main_subloc = cur_max_diff_subloc1,
      Pair_subloc = cur_max_diff_subloc2,
      CoSublocType = 'Significant',
      MaxGroupSize = cur_max_diff_group_size,
      AssignedGroupSize = length(cur_gname)-length(cur_diff_unassign_group),
      Score = cur_max_diff_group_size / ( length(cur_gname)-length(cur_diff_unassign_group) )
    ))
    
  }
  
  cur_df$ExpName <- exp_name
  cur_df$ExpCondition <- cur_exp_cond
  
  performance <- rbind(performance, cur_df)
  
  print(paste(cur_exp_cond, length(unique(cur_df$ComplexName)), nrow(cur_df), sep = ' '))
  print(head(cur_df))
  print('=====================================================================')
}




