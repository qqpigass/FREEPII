library(kaos)

seq_path <- '/FREEPII_github/UniProt-data/uniprot_seq_Human.rds'

res <- 16

temp <- readRDS(seq_path)
out <- setNames(data.frame(sapply( c(1:nrow(temp)), 
                                   function(x) vectorize(cgr( strsplit(temp$Sequence[x], '')[[1]], sf = 0.863271, res = res)) )), 
                temp$UniProt_entry)

dir.create('/FREEPII_github/FCGR-data', showWarnings = FALSE)
saveRDS(out, '/FREEPII_github/FCGR-data/uniprot_seq_Human_FCGR_16x.rds')

