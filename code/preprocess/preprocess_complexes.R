# library(tidyverse)
library(magrittr)
library(reshape2)
library(plyr)

#===================================
#          Prepare data
#===================================
#### Set path
## Complex path
complex_path <- '/FREEPII_github/Protein complex-data/allComplexes_Corum.txt'
complex_path

## UniProt path
uniprot_gene_path <- '/FREEPII_github/UniProt-data/uniprot_gene_Human.rds'
uniprot_gene_path


#### Read complex data
complex_df_raw <- read.table(complex_path, sep='\t', quote='', header = T)
nrow(complex_df_raw)
colnames(complex_df_raw)
unique(complex_df_raw$Organism)

## Extract data of specified organism
complex_df_raw <- unique(complex_df_raw[complex_df_raw$Organism=='Human', ])
nrow(complex_df_raw)
unique(complex_df_raw$Organism)

## Select needed columns
complex_df_raw <- unique( complex_df_raw[, c('ComplexName', 'subunits.UniProt.IDs.', 'subunits.Gene.name.', 'subunits.Entrez.IDs.')] )
nrow(complex_df_raw)
colnames(complex_df_raw)

## Re-name columns
colnames(complex_df_raw)[2:ncol(complex_df_raw)] <- c('UniProt_entry', 'Gene_name', 'Gene_entrez')
colnames(complex_df_raw)

## Unique items and empty in data
rbind( sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) length(unique(complex_df_raw[,x])) ), 
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df_raw[,x]=='') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df_raw[,x]==0) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df_raw[,x]=='-') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df_raw[,x]=='None') ) )


#### Read UniProt data
uniprot_df <- readRDS(uniprot_gene_path)
nrow(uniprot_df)
colnames(uniprot_df)

## Unique items and empty in data
rbind( sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'), function(x) length(unique(uniprot_df[,x])) ),
       sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'), function(x) sum(uniprot_df[,x]=='') ),
       sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'), function(x) sum(uniprot_df[,x]==0) ),
       sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'), function(x) sum(uniprot_df[,x]=='-') ),
       sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'), function(x) sum(uniprot_df[,x]=='None') ))



#===================================================
#           Release genes in complex df
#===================================================
sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(unique(unlist(strsplit(complex_df_raw[,'UniProt_entry'], ';')))!='') ) 

## Match split in each row
complex_df <- data.frame(matrix(0, ncol = 4, nrow = 0))
for ( ComplexName in unique(complex_df_raw$ComplexName) ) {
  uniprot_id <- strsplit(complex_df_raw[complex_df_raw$ComplexName==ComplexName, 'UniProt_entry'], ';')[[1]]
  g_name <- strsplit(complex_df_raw[complex_df_raw$ComplexName==ComplexName, 'Gene_name'], ';')[[1]]
  g_entrez <- strsplit(complex_df_raw[complex_df_raw$ComplexName==ComplexName, 'Gene_entrez'], ';')[[1]]
  
  uniprot_id <- uniprot_id[uniprot_id!='']
  g_name <- g_name[g_name!='']
  g_entrez <- g_entrez[g_entrez!='']
  
  temp <- data.frame(ComplexName = ComplexName,
                     UniProt_entry = uniprot_id,
                     Gene_name = g_name,
                     Gene_entrez = g_entrez)
  complex_df <- rbind(complex_df, temp)
}
nrow(complex_df)
colnames(complex_df)
rbind( sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) length(unique(complex_df[,x])) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]=='') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]==0) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]=='-') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]=='None') ))

## Filter None in UniProt_entry and Gene_name
complex_df <- unique( complex_df[(complex_df$UniProt_entry!='None')&(complex_df$Gene_name!='None'), ] )
nrow(complex_df) # 14169
rbind( sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) length(unique(complex_df[,x])) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]=='') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]==0) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]=='-') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_entrez'), function(x) sum(complex_df[,x]=='None') ))



#===================================================
#     Match UniProt ID and gene name to complex
#===================================================
#### Combined Uniprot info to complex
## Overlap between names in complex and names in Uniprot
cbind( rbind( sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'),
                      function(x) sum(unique(complex_df$UniProt_entry) %in% unique(as.character(as.matrix(uniprot_df[,x])))) ),
              sapply( c('UniProt_entry', 'Gene_name', 'Gene_name_primary'),
                      function(x) sum(unique(complex_df$Gene_name) %in% unique(as.character(as.matrix(uniprot_df[,x])))) ) ), 
       c(length(unique(complex_df$UniProt_entry)), length(unique(complex_df$Gene_name))) )

## Merged by Gene_name
complex_df <- unique( merge(complex_df, uniprot_df[ ,c('UniProt_entry', 'Gene_name', 'Gene_name_primary')], by = 'Gene_name', all.x = T) )
complex_df[is.na(complex_df)] <- 'None'
nrow(complex_df)

## Combined UniProt_entry
a <- complex_df[complex_df$UniProt_entry.x==complex_df$UniProt_entry.y, ]
b <- complex_df[complex_df$UniProt_entry.x!=complex_df$UniProt_entry.y, ]
b <- melt(b, id=c('ComplexName', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), variable.name = 'Group', value.name = 'UniProt_entry.x')
unique(b$Group)

complex_df <- unique( rbind(a[,c('ComplexName', 'UniProt_entry.x', 'Gene_name', 'Gene_name_primary', 'Gene_entrez')], 
                            b[,c('ComplexName', 'UniProt_entry.x', 'Gene_name', 'Gene_name_primary', 'Gene_entrez')]) )
nrow(complex_df)
colnames(complex_df)[2] <- 'UniProt_entry'
colnames(complex_df)

## Unique items and empty in data
rbind( sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) length(unique(complex_df[,x])) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]=='') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]==0) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]=='-') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]=='None') ))

## Remove UniProt_entry==None
complex_df <- unique(complex_df[complex_df$UniProt_entry!='None', ])
nrow(complex_df) # 15613
rbind( sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) length(unique(complex_df[,x])) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]=='') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]==0) ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]=='-') ),
       sapply( c('ComplexName', 'UniProt_entry', 'Gene_name', 'Gene_name_primary', 'Gene_entrez'), function(x) sum(complex_df[,x]=='None') ))



#==========================================================
#     Count complex sub-units based on different names
#==========================================================
#### Subunit number
length(unique(complex_df$ComplexName))
freq_uniprot_entry <- setNames(aggregate(data = complex_df, UniProt_entry ~ ComplexName,     function(x) sum(unique(x)!='None')), c('ComplexName', 'UniProt_entry_N'))
freq_gname         <- setNames(aggregate(data = complex_df, Gene_name ~ ComplexName,         function(x) sum(unique(x)!='None')), c('ComplexName', 'Gene_name_N'))
freq_gname_primary <- setNames(aggregate(data = complex_df, Gene_name_primary ~ ComplexName, function(x) sum(unique(x)!='None')), c('ComplexName', 'Gene_name_primary_N'))
print(c(nrow(freq_uniprot_entry), nrow(freq_gname), nrow(freq_gname_primary)))

complex_df <- unique( join(join(join(complex_df, freq_uniprot_entry), freq_gname), freq_gname_primary) )
nrow(complex_df)
summary(complex_df)
summary(complex_df)[,6:ncol(complex_df)]


#### Filter complex by size (Gene_name)
complex_df_sub <- unique(complex_df[complex_df$Gene_name_N > 2,])
nrow(complex_df_sub)
summary(complex_df_sub)
summary(complex_df_sub)[,6:ncol(complex_df_sub)]



#===================
#     Save data
#===================
saveRDS(complex_df_sub, '/FREEPII_github/Protein complex-data/Complexes_gene_Human_filter.rds')




