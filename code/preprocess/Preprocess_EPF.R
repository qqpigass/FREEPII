library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)

#===================================
#      Read data and normalize
#===================================
quant_type <- 'iBAQ'

#### Set path
exp_dir <- '/FREEPII_github/EPF-data'
exp_name <- list.files(exp_dir)[!grepl('.R', list.files(exp_dir))]
exp_file <- sapply( file.path(exp_dir, exp_name), function(x) list.files(x)[!grepl('.rds', list.files(x))], simplify = F, USE.NAMES = F )
exp_path <- unlist(sapply( c(1:length(exp_name)), 
                           function(y) sapply( 
                             file.path(exp_dir, exp_name[y], exp_file[[y]]), 
                             function(x) paste(x, list.files(x)[grepl(quant_type, list.files(x))&grepl('tsv', list.files(x))], sep = '/'),
                             simplify = F, USE.NAMES = F
                           ),
                           simplify = F, USE.NAMES = F
))

complex_path <- '/FREEPII_github/Protein complex-data/Complexes_gene_Human_filter.rds'


#### Read complex
complex_df <- readRDS(complex_path)
nrow(complex_df)
colnames(complex_df)

complex_df <- unique(complex_df[,c('ComplexName', 'Gene_name')])
Complex_gnames <- sort(unique(complex_df$Gene_name))


#### Read EPF, normalize, match protein group to gene, and return matrix with row name = gene name
for ( path in exp_path ) {
  #### Pre-processing mat
  mat = as.matrix(read.delim(path))
  meta =  read.delim(paste(c(head(strsplit(path, '/')[[1]], -1), 'metadata.tsv'), collapse = '/'))
  
  ## Replace missing values in EPF with zero
  mat[!is.finite(mat)] = NA
  mat[is.nan(mat)] = NA
  mat[is.na(mat)] = 0
  
  ## Remove row with all zero
  meta =  meta[rowSums(mat)!=0, ]
  mat = mat[rowSums(mat)!=0, ]
  
  ## Row-wise normalize
  mat = (mat / rowSums(mat))
  mat[!is.finite(mat)] = NA
  mat[is.nan(mat)] = NA
  mat[is.na(mat)] = 0
  
  ## keep one row per gene
  gene_map = meta %>%
    select(Majority.protein.IDs, Gene.names) %>%
    mutate(Protein_IDs = Majority.protein.IDs) %>%
    mutate(Gene_name = strsplit(Gene.names, ';')) %>%
    unnest(Gene_name) %>%
    select(Protein_IDs, Gene_name) %>%
    drop_na()
  
  Gnames = unique(gene_map$Gene_name)
  gene_mat = matrix(NA, nrow = length(Gnames), ncol = ncol(mat), dimnames = list(Gnames, colnames(mat)))
  n_fractions = rowSums(!is.na(mat) & is.finite(mat) & mat != 0)
  out_map = data.frame()
  
  for (gname in Gnames) {
    protein_id = gene_map$Protein_IDs[gene_map$Gene_name == gname]
    
    # pick the best protein for this replicate
    n_fraction = n_fractions[protein_id]
    best = names(which(n_fraction == max(n_fraction))) %>% first()
    gene_mat[gname, ] = mat[best, ]
    
    # save the mapping
    out_map %<>% bind_rows(data.frame(Protein_IDs = best, Gene_name = gname))
  }
  
  
  #### Remove duplicates
  ## First ordered by name
  order = rownames(gene_mat) %>% order()
  gene_mat = gene_mat[order, ]
  out_map = out_map[order, ]
  
  ## Second ordered by 'if in complexes'. True values will be gather ahead, so later part will be those False
  order = rownames(gene_mat) %in% Complex_gnames %>% order(decreasing = TRUE)
  gene_mat = gene_mat[order, ]
  out_map = out_map[order, ]
  
  # ## Third, drop duplicated rows
  # drop = which(duplicated(gene_mat))
  # gene_mat = gene_mat[-drop, ]
  # out_map = out_map[-drop, ]
  
  ## Four, re-name columns and rows
  colnames(gene_mat) <- paste0('F', c(1:ncol(gene_mat)))
  rownames(out_map) <- NULL
  
  ## set attribute on the gene matrix
  attr(gene_mat, 'out_map') = out_map
  
  cur_exp_name <- strsplit(path, '/')[[1]][which(grepl('PXD', strsplit(path, '/')[[1]]))]
  cur_name <- strsplit(path, '/')[[1]][which(grepl('PXD', strsplit(path, '/')[[1]]))+1]
  
  saveRDS(gene_mat, gsub(paste0(quant_type, '.tsv'), paste0(cur_name, '.rds'), path))
  
  cat(paste(
    paste0(cur_name, ': '),
    paste0(nrow(gene_mat), ' ', ncol(gene_mat)),
    paste0(nrow(out_map), ' ', ncol(out_map)),
    strrep('=', 60), 
    sep = '\n'
  ))
}




