options(stringsAsFactors=F)

#========= Paths =========#
base_dir <- list(
   hpc='/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
   mnt='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
   local='/Users/lnguyen/Documents/P0013_WGS_patterns_Diagn/'
)

for(i in base_dir){
   if(dir.exists(i)){
      base_dir <- i
      break
   }
}
rm(i)

#========= Args =========#
args <- commandArgs(trailingOnly=T)

training_data_path <- args[1]
seed <- args[2]
out_dir <- args[3]
sel_sample_cv_dir <- args[4]
test_data_path <- args[5]

core_functions_path <- if(is.na(args[6])){ 
   paste0(base_dir,'/CHORDv2/processed/training/scripts/train_rf/trainCoreFunctions.R')
} else {
   args[6] 
}

# training_data_path = paste0(base_dir,'/CHORDv2/training/scripts/train_rf/training_data.txt.gz')
# core_functions_path = paste0(base_dir,'/CHORDv2/training/scripts/train_rf/trainCoreFunctions.R')
# seed = 1
# sel_sample_cv_dir = paste0(base_dir,'/CHORDv2/training/models/2.00_expandedIndels/cv_sel_samples/')
# out_dir = paste0(base_dir,'/CHORDv2/training/models/2.00_expandedIndels/')
# test_data_path = paste0(base_dir,'/CHORDv2/training/models/2.00_expandedIndels/outer_cv/02/test_data.txt.gz')

#========= Main =========#
source(core_functions_path)

## Make whitelist
training_data <- read.delim(training_data_path, stringsAsFactors=T)

sel_sample_cv_summary <- whitelistSamplesFromCvPreds(sel_sample_cv_dir)
write.table(
   sel_sample_cv_summary,
   paste0(out_dir,'/sel_sample_cv_summary.txt'),
   sep='\t',quote=F,row.names=F
)

whitelist_samples <- sel_sample_cv_summary[sel_sample_cv_summary$whitelist==TRUE,'sample']

## Train final model
training_data_whitelisted <- training_data[rownames(training_data) %in% whitelist_samples,]
write.table(
   training_data_whitelisted,
   gzfile(paste0(out_dir,'/training_data_whitelisted.txt.gz')),
   sep='\t',quote=F
)

rf <- randomForestTrain(training_data_whitelisted, seed=seed)
rf$seed <- seed
saveRDS(rf, paste0(out_dir,'/rf_model.rds'))

#========= Predict on test set if exists =========#
if(!is.na(test_data_path)){
   test_data <- read.delim(test_data_path, stringsAsFactors=T)
   test_pred <- as.data.frame(predict(rf, test_data, 'prob'))
   test_pred$response <- test_data$response
   
   test_pred <- cbind(sample=rownames(test_pred), test_pred)
   rownames(test_pred) <- NULL
   
   write.table(
      test_pred, paste0(out_dir,'/test_preds.txt.gz'),
      sep='\t',quote=F,row.names=F
   )
}