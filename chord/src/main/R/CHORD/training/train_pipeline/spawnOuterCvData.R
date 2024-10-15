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

core_functions_path <- if(is.na(args[4])){ 
   paste0(base_dir,'/CHORDv2/processed/training/scripts/train_rf/trainCoreFunctions.R')
} else {
   args[4] 
}

# training_data_path = paste0(base_dir,'/CHORDv2/training/scripts/train_rf/training_data.txt.gz')
# core_functions_path = paste0(base_dir,'/CHORDv2/training/scripts/train_rf/trainCoreFunctions.R')
# seed = 1
# out_dir = paste0(base_dir,'/CHORDv2/training/models/2.00_expandedIndels/outer_cv/')

#========= Main =========#
source(core_functions_path)

## Make whitelist
training_data <- read.delim(training_data_path, stringsAsFactors=T)
set.seed(seed)
spawnOuterCvData(training_data, out_dir)
