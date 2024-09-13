options(stringsAsFactors=F)

## Path prefixes ================================
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

## Merge cancer types ================================
metadata <- (function(){
   df1 <- read.delim(paste0(base_dir,'/CHORDv2/processed/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))
   df1 <- df1[,c('sample_id','primary_tumor_location','in_training_set','response')]
   colnames(df1) <- c('sample','cancer_type','in_training_set','response')
   df1$cohort <- 'HMF'
   
   df2 <- read.delim(paste0(base_dir,'/CHORDv2/processed/training/scripts/sel_training_samples/PCAWG_2020/metadata_training_samples_ann.txt'))
   df2 <- df2[,c('sample','cancer_type','in_training_set','response')]
   df2$cohort <- 'PCAWG'
   
   list(HMF=df1, PCAWG=df2)
})()

# table(
#    metadata$cohort,
#    metadata$in_training_set
# )

metadata <- list(
   HMF=metadata$HMF[metadata$HMF$in_training_set,],
   PCAWG=metadata$PCAWG[metadata$PCAWG$in_training_set,]
)

holdout_samples$HMF <- (function(){
   metadata_ss <- metadata$HMF
   
   ct_counts <- sort(table(
      metadata_ss[
         metadata_ss$response!='none',
         'cancer_type'
      ]
   ))
   ct_whitelist <- names(ct_counts)[ct_counts>=2]
   
   metadata_ss$cancer_type[ !(metadata_ss$cancer_type %in% ct_whitelist) ] <- 'Other'
   #table(metadata_ss$cancer_type)
   
   out <- split(metadata_ss$sample, metadata_ss$cancer_type)
   names(out) <- gsub(' |/','_',names(out))
   
   return(out)
})()

holdout_samples$PCAWG <- (function(){
   metadata_ss <- metadata$PCAWG
   # ct_counts <- sort(table(
   #    metadata_ss[
   #       metadata_ss$response!='none',
   #       'cancer_type'
   #    ]
   # ))
   # ct_whitelist <- names(ct_counts)[ct_counts>=2]
   
   metadata_ss$cancer_type[ !(metadata_ss$cancer_type %in% names(holdout_samples$HMF)) ] <- 'Other'
   
   out <- split(metadata_ss$sample, metadata_ss$cancer_type)
   names(out) <- gsub(' |/','_',names(out))
   return(out)
})()

## Spawn training ================================
training_data <- read.delim(paste0(base_dir,'/CHORDv2/processed/training/training_data/mTrain_mergedDelMh2-5.txt.gz'))

models_dir <- paste0(base_dir,'/CHORDv2/processed/training/models_ct_holdout/')
for(i in names(holdout_samples)){
   message(i)
   model_subdir <- paste0(models_dir,'/',i,'/')
   dir.create(model_subdir, showWarnings=F)

   ## Feature matrix
   training_data_ss <- training_data[
      !(rownames(training_data) %in% holdout_samples[[i]])
   ,]

   training_data_path <- paste0(model_subdir,'/training_data.txt.gz')
   write.table(training_data_ss, gzfile(training_data_path), sep='\t', quote=F)

   ## Submit script
   submit_sh <- sprintf("#!/bin/bash
sh /hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORDv2/processed/training/scripts/train_rf/trainPipeline.sh \\
%s \\
%s \\
1 0
", sub('/Users/lnguyen','',training_data_path), sub('/Users/lnguyen','',model_subdir))

   fileConn <- file(paste0(model_subdir,'/run_train_pipeline.sh'))
   writeLines(submit_sh, fileConn)
   close(fileConn)

}
rm(training_data)


## Manually submit pipeline sh file














