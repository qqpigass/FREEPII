

#### Settings
split_dir <- '/FREEPII_github/Split-data'
exp_name <- list.files(split_dir)[!grepl('R', list.files(split_dir))]
split_file <- file.path(split_dir, exp_name)
split_path <- unlist(sapply(split_file, function(x) file.path(x, list.files(x)), USE.NAMES = F))
split_path <- split_path[grepl('cv_split', split_path)]
split_path <- split_path[grepl('.rds', split_path)]
split_path


#### Read split
for ( path_idx in c(1:length(split_path)) ) {
  path <- split_path[path_idx]
  cur_exp_name <- tail(strsplit(path, '/')[[1]],2)[1]
  cur_exp_cond <- gsub('.rds', '', tail(strsplit(tail(strsplit(path, '/')[[1]],2)[2], '_')[[1]], 1))
  
  cur_cv <- readRDS(path)
  
  for (cv_idx in c(1:length(cur_cv))) {
    if (cv_idx == 1) {
      temp <- rbind(data.frame(cur_cv[[cv_idx]]$train, CV_fold = cv_idx), 
                    data.frame(cur_cv[[cv_idx]]$test, CV_fold = cv_idx))
    } else {
      temp <- rbind(temp, 
                    rbind(data.frame(cur_cv[[cv_idx]]$train, CV_fold = cv_idx), 
                          data.frame(cur_cv[[cv_idx]]$test, CV_fold = cv_idx)))
    }
    
    cat(paste(
      paste(cur_exp_name, cur_exp_cond, paste0('cv', cv_idx, ':'), sep = ' '),
      paste0('Train TP: ', sum(temp[temp$Type=='Train', ]$Label==1), ' ', 'Train TN: ', sum(temp[temp$Type=='Train', ]$Label==0)),
      paste0('Test TP: ', sum(temp[temp$Type=='Test', ]$Label==1), ' ', 'Test TN: ', sum(temp[temp$Type=='Test', ]$Label==0)),
      strrep('=', 100),
      sep = '\n'
    ), '\n')
  }
  
  write.csv(temp, paste(split_dir, cur_exp_name, gsub('.rds$', '.csv', basename(path)), sep = '/'), row.names = F)
  rm(temp)
}





