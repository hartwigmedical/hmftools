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
out_path <- args[3]

core_functions_path <- if(is.na(args[4])){ 
   paste0(base_dir,'/CHORDv2/processed/training/scripts/train_rf/trainCoreFunctions.R')
} else {
   args[4] 
}

# training_data_path = paste0(base_dir,'/CHORDv2/training/scripts/train_rf/training_data.txt.gz')
# core_functions_path = paste0(base_dir,'/CHORDv2/training/scripts/train_rf/trainCoreFunctions.R')
# seed = 1

#========= Main =========#
source(core_functions_path)
training_data <- read.delim(training_data_path, stringsAsFactors=T)

rf_cv <- randomForestCv(training_data, train.function=randomForestTrain, seed=seed)
preds <- do.call(rbind,lapply(rf_cv,`[[`,'pred'))
write.table(preds, gzfile(out_path), sep='\t',quote=F)


##
# rds <- readRDS('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORDv2/training/models/2.00_expandedIndels/cv_sel_samples/rf_1_bak.rds')
# m1 <- do.call(rbind,lapply(rds,`[[`,'pred'))
# 
# m2 <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORDv2/training/models/2.00_expandedIndels/cv_sel_samples/pred_1.txt.gz')
