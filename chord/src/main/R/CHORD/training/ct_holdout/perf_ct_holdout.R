options(stringsAsFactors=F)

library(CHORD)
library(mutSigExtractor)
library(mltoolkit)
library(reshape2)
library(ggplot2)


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

## Performance ================================
## Load features ----------------------------
contexts <- list(
   HMF=readRDS(paste0(base_dir,'/datasets/processed/HMF_DR010_DR047/matrices/contexts_merged.rds')),
   PCAWG=readRDS(paste0(base_dir,'/datasets/processed/PCAWG_2020/matrices/contexts/contexts_merged.rds'))
)

trans_func <- function(x){ 
   ## Merge del mh bimh 2-5
   indels <- x$indel
   indel_types <- c('del.rep','ins.rep','del.mh','ins.mh','del.none','ins.none')
   indels_split <- splitDfRegex(indels, indel_types)
   names(indels_split) <- indel_types
   
   counter <- 0
   indels_custom <- lapply(indels_split, function(i){
      counter <<- counter + 1
      df <- as.data.frame(rowSums(i))
      colnames(df) <- indel_types[counter]
      return(df)
   })
   
   indels_custom$del.mh <- with(indels_split,{
      cbind(
         del.mh['del.mh.bimh.1'],
         'del.mh.bimh.2.5'=rowSums(del.mh[!(colnames(del.mh) %in% 'del.mh.bimh.1')])
      )
   })
   
   ## Add new indels back to list
   l <- x[c('snv','indel','sv')] 
   l$indel <- do.call(cbind, unname(indels_custom))
   m <- transformContexts(l, simplify.types='snv', rel.types='all')
   return(m)
}

contexts <- lapply(contexts, trans_func)

contexts <- within(contexts,{
   HMF <- HMF[rownames(HMF) %in% metadata$HMF$sample,]
   HMF$response <- metadata$HMF[ match(rownames(HMF), metadata$HMF$sample),'response' ]
   
   PCAWG <- PCAWG[rownames(PCAWG) %in% metadata$PCAWG$sample,]
   PCAWG$response <- metadata$PCAWG[ match(rownames(PCAWG), metadata$PCAWG$sample),'response' ]
})

## Load RF models ----------------------------
rf_models <- (function(){
   dirs <- list.dirs(paste0(base_dir,'/CHORDv2/processed/training/models_ct_holdout/'), recursive=F, full.names=T)
   l <- lapply(dirs, function(i){
      readRDS(paste0(i,'/rf_model.rds'))
   })
   names(l) <- basename(dirs)
   return(l)
})()

## Order RF models by BRCA1/2 def freq
cancer_types <- (function(){
   ct_counts <- sort(table(
      metadata$HMF[
         metadata$HMF$response!='none',
         'cancer_type'
      ]
   ))
   names(ct_counts)[ct_counts>=2]
})()
cancer_types <- c('Other',cancer_types)
cancer_types <- gsub(' |/','_',cancer_types)

rf_models <- rf_models[cancer_types]

## Predict, calc perf ----------------------------
calcPerformance <- function(features, holdout.samples){
   #features=contexts$PCAWG
   #holdout.samples=holdout_samples$PCAWG
   
   preds <- lapply(names(rf_models), function(i){
      #i=names(rf_models)[1]
      m <- features[rownames(features) %in% holdout.samples[[i]],]
      
      if(nrow(m)==0){ return(NULL) }
      
      pred <- as.data.frame(predict(
         rf_models[[i]], 
         m[,colnames(m)!='response'], 
         type='prob'
      ))
      
      pred$hrd <- pred$BRCA1 + pred$BRCA2
      pred$response <- m$response
      pred$cancer_type <- i
      
      return(pred)
   })
   names(preds) <- names(rf_models)
   
   preds <- preds[ !sapply(preds, is.null) ]
   
   preds[['(Aggregate)']] <- do.call(rbind, unname(preds))
   #subset(preds$All, response!='none' & hrd<0.5)
   
   
   confusions <- lapply(names(preds), function(i){
      #i='Biliary'
      pred <- preds[[i]]
      response <- ifelse(pred$response!='none',TRUE,FALSE)
      
      confusion <- as.data.frame(rbind(
         confusionMatrix(pred$hrd,response, cutoff=0.5),
         confusionMatrix(pred$BRCA1,response, cutoff=0.5),
         confusionMatrix(pred$BRCA2,response, cutoff=0.5)
      ))
      cbind(
         cancer_type=i,
         pred_class=c('HRD','BRCA1_type','BRCA2_type'),
         confusion
      )
   })
   confusions <- do.call(rbind, confusions)
   #confusions$cutoff <- NULL
   
   perfs <- cbind(
      confusions, 
      calcPerf(confusions, metrics=c('fpr','fnr'), add.start.end.values=F)
   )
   perfs$fpr[is.na(perfs$fpr)] <- 0
   perfs$fnr[is.na(perfs$fnr)] <- 0
   
   perfs <- perfs[,colnames(perfs)!='cutoff']
   
   perfs$fpr_label <- paste0(round(perfs$fpr,2),' (',perfs$fp,' / ',perfs$fp + perfs$tn,')')
   perfs$fnr_label <- paste0(round(perfs$fnr,2),' (',perfs$fn,' / ',perfs$fn + perfs$tp,')')
   
   return(perfs)
}

perf <- (function(){
   df1 <- calcPerformance(contexts$HMF, holdout_samples$HMF)
   df1$cohort <- 'HMF'
   
   df2 <- calcPerformance(contexts$PCAWG, holdout_samples$PCAWG)
   df2$cohort <- 'PCAWG'
   
   rbind(df1,df2)
})()



## Plot ================================
p_perf <- (function(){
   pd <- (function(){
      l <- list(
         fpr=perf[,c('cancer_type','pred_class','fpr','fpr_label','cohort')],
         fnr=perf[,c('cancer_type','pred_class','fnr','fnr_label','cohort')]
      )
      
      l <- lapply(names(l), function(i){
         df <- l[[i]]
         colnames(df)[3:4] <- c('value','label')
         df$metric <- i
         return(df)
      })
      
      do.call(rbind, l)
   })()
   
   pd$cancer_type <- factor(pd$cancer_type, rev(unique(pd$cancer_type)))
   pd$pred_class <- factor(pd$pred_class, unique(pd$pred_class))
   pd$metric <- factor(pd$metric, unique(pd$metric))
   
   ggplot(pd[pd$pred_class=='HRD',], aes(x=cancer_type, y=value, fill=metric, group=metric, label=label)) +
      facet_wrap(~cohort, ncol=1) +
      geom_bar(stat='identity', position=position_dodge(width=0.5), width=0.5) +
      geom_text(aes(y=0), position=position_dodge(width=0.5), vjust=0.5, hjust=0, angle=90) +
      geom_vline(xintercept=1.5, size=0.5, linetype='dotted') +
      labs(y='Metric value', x='Held-out cancer type') +
      theme_bw() +
      theme(
         axis.text.x=element_text(angle=30,hjust=1,vjust=1)
      )
})()

pdf(paste0(base_dir,'/CHORDv2/processed/training/scripts/ct_holdout/perf_ct_holdout.pdf'), 8, 6)
plot(p_perf)
dev.off()











