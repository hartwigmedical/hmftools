library(mutSigExtractor)
library(randomForest)
library(mltoolkit)
library(reshape2)
library(ggplot2)
library(scales)
library(grid)

options(stringsAsFactors=F)

#========= Path prefixes =========#
base_dir <- list(
   hpc='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/',
   mnt='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/',
   local='/Users/lnguyen/Documents/Luan_projects/'
)

for(i in base_dir){
   if(dir.exists(i)){ 
      base_dir <- i 
      break
   }
}

#========= Misc functions =========#
forceDfOrder <- function(df){
   as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
}

#========= Load data =========#
rf <- (function(){
   dir <- paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/')
   rf <- readRDS(paste0(dir,'/rf_model.rds'))
   rf$trans.func <- function(x){
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
   rf$dir <- dir
   return(rf)
})()

# rf <- (function(){
#    #dir <- paste0(base_dir,'/CHORDv2/training/models/2.01b_expandedIndels_wilcox1e-10/')
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.01b_expandedIndels_wilcox1e-10/seed_models/seed02/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types='snv',rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()


metadata <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))
sel_samples <- metadata[metadata$in_training_set,'sample_id'] ## Only use samples where BRCA status is known for certain

contexts <- readRDS(paste0(base_dir,'/datasets/HMF_DR010_DR047/matrices/contexts_merged.rds'))
contexts <- contexts[c('snv','indel','sv')]
contexts <- lapply(contexts, function(i){
   i[rownames(i) %in% sel_samples,]
})

#========= Functions =========#
predictHrd <- function(contexts){
   m <- rf$trans.func(contexts)
   pred <- as.data.frame(predict(rf, m,type='prob'))
   pred <- cbind(sample=rownames(pred),pred); rownames(pred) <- NULL
   pred$hrd <- pred$BRCA1 + pred$BRCA2
   pred$response <- metadata[match(pred$sample, metadata$sample_id),'response']
   
   return(pred)
}

#predictHrd(contexts)

resampleFeatureMatrix <- function(m, scale.factor=1){
   #m=contexts$indel
   #scale.factor=1
   
   if(scale.factor==1){ return(m) }
   
   feature_ids <- as.factor(1:ncol(m)) ## Replace feature names by integers to save memory
   
   resampleFeatureVector <- function(v){
      #v=unlist(m[229,])
      total_counts <- sum(v)
      if(total_counts==0){ ## Prevent divide by 0 errors
         return(
            rep(0, length(feature_ids))
         )
      }
      
      prob <- v/total_counts
      table( sample(feature_ids, total_counts*scale.factor, replace=T, prob=prob) )
   }
   
   out <- t(apply(m,1,resampleFeatureVector))
   colnames(out) <- colnames(m)
   
   return(out)
}

#========= Main =========#
calcDownSamplePerf <- function(seed=1, out_path=NULL){
   
   set.seed(seed)
   
   #--------- Resampling ---------#
   message('## Performing resampling...')
   #contexts_downsampled_path <- paste0(base_dir,'/CHORDv2/analysis/downscale_mut_perf/data/contexts_downsampled.rds')
   
   # if(file.exists(contexts_downsampled_path)){
   #    contexts_downsampled <- readRDS(contexts_downsampled_path)
   # } else {
      
   downsample_factors <- 1.5 ^ c(0:14)
   
   counter <- 0
   contexts_downsampled <- lapply(downsample_factors, function(i){
      #i=1
      counter <<- counter + 1
      message('Performing resampling for ','[',counter,']: ',i)
      
      list(
         snv=resampleFeatureMatrix(contexts$snv, 1/i),
         indel=resampleFeatureMatrix(contexts$indel, 1/i),
         sv=resampleFeatureMatrix(contexts$sv, 1/i)
      )
   })
   names(contexts_downsampled) <- paste0('x',signif(downsample_factors,3))
   #    saveRDS(contexts_downsampled, contexts_downsampled_path)
   # }
   
   #--------- Calculate median mut counts ---------#
   message('## Calculating median mut counts...')
   median_mut_counts <- do.call(rbind,lapply(contexts_downsampled, function(i){
      #i=contexts_downsampled[[1]]
      df <- sapply(i, rowSums)
      apply(df,2,median)
   }))
   rownames(median_mut_counts) <- names(contexts_downsampled)
   
   #--------- Calculate performance ---------#
   message('## Making confusion matrices...')
   preds <- lapply(contexts_downsampled, predictHrd)
   names(preds) <- names(contexts_downsampled)
   
   confusionMatrixMcCustom <- function(df){
      
      confusion <- list()
      
      confusion <- confusionMatrixMC(
         df[c('BRCA1','BRCA2','none')], 
         df$response, neg.response='none',cutoff='all'
      )
      
      confusion$hrd <- confusionMatrix( 
         df$hrd, 
         toBinaryResponse(df$response,c('BRCA1','BRCA2'),1,'none',0) 
      )
      
      confusion <- confusion[c('hrd','BRCA1','BRCA2')]
      names(confusion) <- c('HRD','BRCA1-type HRD','BRCA2-type HRD')
      
      confusion <- confusion[sapply(confusion,function(i){ !is.null(i) })]
      
      return(confusion)
   }
   
   #confusions_path <- paste0(base_dir,'/CHORDv2/analysis/downscale_mut_perf/data/confusions.rds')
   
   # if(file.exists(confusions_path)){
   #    confusions <- readRDS(confusions_path)
   # } else {
   counter <- 0
   confusions <- lapply(preds, function(i){
      counter <<- counter + 1
      message('Making confusion matrix for ','[',counter,']: ',names(preds)[counter])
      #i=preds[[1]]
      confusionMatrixMcCustom(i)
   })
   names(confusions) <- names(preds)
   #    saveRDS(confusions, confusions_path)
   # }
   
   perf <- do.call(rbind,lapply(confusions, function(i){
      #i <- confusions[[1]]
      sapply(i, function(j){
         #i=confusion[[1]]
         m <- calcPerfCompound(j, 'pr', metric.names.as.x.y=T)
         calcAUC(m[,'x'],m[,'y'])
      })
   }))
   
   out <- list(
      #contexts_downsampled=contexts_downsampled,
      #confusions=confusions,
      median_mut_counts=median_mut_counts,
      perf=perf
   )
   
   if(is.null(out_path)){
      return(out)
   } else {
      saveRDS(out, out_path)
   }
}

seeds <- 1:10
out_dir <-  paste0(base_dir,'/CHORDv2/analysis/downscale_mut_perf/')


for(i in seeds){
   message('#### SEED ',i,' ####')
   out_path <- paste0(out_dir,'/data/seed_',i,'.rds')
   if(!file.exists(out_path)){
      calcDownSamplePerf(seed=i, out_path)
   }
}

counter <- 0
rds_files <- lapply(
   paste0(out_dir,'/data/seed_',seeds,'.rds'),
   function(i){ 
      counter <<- counter+1
      message('Reading in: ',counter)
      readRDS(i)[c('perf','median_mut_counts')] 
   }
)

#========= Plot =========#
## Perf
perf_melt <- melt(lapply(rds_files,`[[`,'perf'))
colnames(perf_melt) <- c('downsample_factor','pred_class','metric_value','seed')

perf_pd <- with(perf_melt,{
   aggregate(metric_value, by=list(downsample_factor, pred_class), FUN=function(x){
      c(median=median(x), min=min(x), max=max(x))
   })
})
perf_pd <- cbind(perf_pd[1:2], as.data.frame(perf_pd$x))
colnames(perf_pd)[1:2] <- c('downsample_factor','pred_class')
perf_pd <- forceDfOrder(perf_pd)

## Median mut counts
median_mut_counts_pd <- rds_files[[1]]$median_mut_counts ## All seeds have same median mut count
median_mut_counts_pd <- melt(median_mut_counts_pd)
colnames(median_mut_counts_pd) <- c('downsample_factor','mut_type','median_mut_count')
median_mut_counts_pd <- forceDfOrder(median_mut_counts_pd)


p1 <- ggplot(perf_pd, aes(x=downsample_factor,y=median, group=pred_class, color=pred_class)) + 
   geom_hline(yintercept=0.5, linetype='dotted') +
   geom_point() + geom_line() +
   geom_linerange(aes(ymin=min, ymax=max)) +
   ylim(0,1) +
   scale_color_manual(values=c('HRD'='darkgrey','BRCA1-type HRD'='#f58225', 'BRCA2-type HRD'='#69439d')) +
   labs(y='AUPRC (median, min, max)', x='Downsample factor', color='Prediction class') +
   theme_bw()


p2 <- ggplot(median_mut_counts_pd, aes(x=downsample_factor, y=median_mut_count, group=mut_type)) +
   facet_wrap(~mut_type, ncol=1, scales='free_y') +
   geom_point(aes(color=mut_type)) + 
   geom_line(aes(color=mut_type)) +
   geom_text(aes(label=median_mut_count), vjust=-1) +
   scale_y_continuous(expand=c(0.1 ,0.1 ,0.5 ,0.1)) +
   labs(y='Median number of mutations', x='Downsample factor', color='Mut. type') +
   theme_bw()

p_combined <- plot_grid(p1,p2, align='v', axis='tblr', ncol=1)

pdf(paste0(out_dir,'/downsampling_perf.pdf'),10,7)
grid.draw(p_combined)
dev.off()




