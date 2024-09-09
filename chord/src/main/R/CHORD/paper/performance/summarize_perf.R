library(ggplot2)
library(mutSigExtractor)
library(randomForest)
library(cowplot)
library(grid)
library(gridExtra)
library(mltoolkit)
library(reshape2)
library(openxlsx)


options(stringsAsFactors=F)

#========= Paths =========#
base_dir <- list(
   hpc='/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
   mnt='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
   local='/Users/lnguyen/Documents/projects/P0013_WGS_patterns_Diagn/'
)

for(i in base_dir){
   if(dir.exists(i)){
      base_dir <- i
      break
   }
}

#========= Load data =========#
#--------- Metadata ---------#
#metadata_ann <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))
#sel_samples <- subset(metadata_ann, max_purity_biopsy & !has_msi)$sample_id

metadata <- list()

metadata$HMF <- read.delim(paste0(base_dir,'/CHORDv2/processed/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))
metadata$BRCA_EU <- read.delim(paste0(base_dir,'/CHORD/processed/ICGC/analysis/hrdetect_output/hrdetect.BRCA_EU.txt'))
metadata$PCAWG <- read.delim(paste0(base_dir,'/datasets/processed/PCAWG_2020/metadata/pcawg_metadata_ann.txt.gz'))

metadata$PCAWG_training <- read.delim(paste0(base_dir,'/CHORDv2/processed/training/scripts/sel_training_samples/PCAWG_2020/metadata_training_samples_ann.txt'))
#metadata$PCAWG_training <- subset(metadata$PCAWG_training)

#--------- Features ---------#
contexts <- list()

## HMF
contexts$HMF <- (function(){
   l <- readRDS(paste0(base_dir,'/datasets/processed/HMF_DR010_DR047/matrices/contexts_merged.rds'))
   
   m_sv_mh <- read.delim(paste0(base_dir,'/CHORD/processed/HMF_DR010_DR047/scripts/get_sv_type_and_homology/sv_mh_contexts.txt.gz'))
   m_sv_mh <- m_sv_mh[rownames(l[[1]]),]
   
   l$sv_mh <- m_sv_mh
   
   return(l)
})()

## BRCA-EU
contexts$BRCA_EU <- readRDS(paste0(base_dir,'/CHORD/processed/ICGC/matrices/BRCA-EU/merged_contexts.rds'))
contexts$BRCA_EU.HMF_pipeline <- readRDS(paste0(base_dir,'/CHORD/processed/BRCA_EU/matrices/BRCA_EU_hmfPipeline_2/merged_contexts.rds'))

## PCAWG
contexts$PCAWG <- readRDS(paste0(base_dir,'/datasets/processed/PCAWG_2020/matrices/contexts/contexts_merged.rds'))


####################################################################################################
forceDfOrder <- function(df){
   as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
}

main <- function(rf, hmf.data.only=F, export.preds.only=F){
   
   # rf=(function(){
   #    dir <- paste0(base_dir,'/CHORDv2/processed/training/models/2.02_mergedDelMh2-5/')
   #    #dir <- paste0(base_dir,'/CHORDv2/processed/training/models/2.02_mergedDelMh2-5/seed_models/seed04/')
   #    rf <- readRDS(paste0(dir,'/rf_model.rds'))
   #    rf$trans.func <- function(x){
   #       ## Merge del mh bimh 2-5
   #       indels <- x$indel
   #       indel_types <- c('del.rep','ins.rep','del.mh','ins.mh','del.none','ins.none')
   #       indels_split <- splitDfRegex(indels, indel_types)
   #       names(indels_split) <- indel_types
   # 
   #       counter <- 0
   #       indels_custom <- lapply(indels_split, function(i){
   #          counter <<- counter + 1
   #          df <- as.data.frame(rowSums(i))
   #          colnames(df) <- indel_types[counter]
   #          return(df)
   #       })
   # 
   #       indels_custom$del.mh <- with(indels_split,{
   #          cbind(
   #             del.mh['del.mh.bimh.1'],
   #             'del.mh.bimh.2.5'=rowSums(del.mh[!(colnames(del.mh) %in% 'del.mh.bimh.1')])
   #          )
   #       })
   # 
   #       ## Add new indels back to list
   #       l <- x[c('snv','indel','sv')]
   #       l$indel <- do.call(cbind, unname(indels_custom))
   #       m <- transformContexts(l, simplify.types='snv', rel.types='all')
   #       return(m)
   #    }
   #    rf$dir <- dir
   #    return(rf)
   # })()
   # hmf.data.only=F
   # export.preds.only=F
   
   cv_fold_dirs <- list.dirs(paste0(rf$dir,'/outer_cv/'), recursive=F)
   out_dir <- paste0(rf$dir,'/plots/')
   dir.create(out_dir, recursive=T, showWarnings=F)
   
   #========= Load data from models =========#
   message('Loading data from models...')
   cv_rfs <- lapply(paste0(cv_fold_dirs,'/rf_model.rds'), readRDS)
   
   counter <- 0
   preds <- lapply(contexts, function(i){
      counter <<- counter+1
      #i=contexts[[1]]
      
      if(hmf.data.only & names(contexts)[counter]!='HMF'){ return() }
      
      m <- rf$trans.func(i)
      pred <- as.data.frame(predict(rf, m,type='prob'))
      pred <- cbind(sample=rownames(pred),pred); rownames(pred) <- NULL
      pred$hrd <- pred$BRCA1 + pred$BRCA2
      
      indel_load <- rowSums(i$indel)
      indel_rep_load <- rowSums(i$indel[,grep('rep',colnames(i$indel))])
      sv_load <- rowSums(i$sv)
      
      pred$hr_status <- ifelse(pred$hrd>=0.5,'HR_deficient','HR_proficient')
      pred$hr_status[indel_load<50 | indel_rep_load>=14000] <- 'cannot_be_determined'
      
      pred$hrd_type <- c('BRCA1_type','BRCA2_type')[ max.col(pred[,c('BRCA1','BRCA2')]) ]
      pred$hrd_type[ indel_load<50 | sv_load<30 | indel_rep_load>=14000 ] <- 'cannot_be_determined'
      pred$hrd_type[ pred$hr_status=='HR_proficient' ] <- 'none'
      
      pred$indel_load <- indel_load
      pred$indel_rep_load <- indel_rep_load
      pred$sv_load <- sv_load
      
      pred <- cbind(pred,m)
      
      return(pred)
   })
   if(hmf.data.only){ preds <- preds['HMF'] }
   
   ## Outer CV
   preds$HMF_CV <- (function(){
      
      test_pred <- lapply(paste0(cv_fold_dirs,'/test_preds.txt.gz'),read.delim)
      test_pred <- do.call(rbind,test_pred)
      
      test_pred$hrd <- test_pred$BRCA1 + test_pred$BRCA2
      
      return(test_pred)
   })()
   preds$HMF_CV <- preds$HMF_CV[,c('sample','BRCA1','BRCA2','none','hrd','response')]
   
   ## CV 100x for filtering samples
   sel_sample_cv_summary <- read.delim(paste0(rf$dir,'/sel_sample_cv_summary.txt'))
   
   #========= Add BRCA response =========#
   message('Annotating response to predictions...')
   
   #--------- HMF ---------#
   getAnnotationsHmf <- function(sample_names){
      #sample_names=preds$HMF$sample
      whitelist_samples <- subset(metadata$HMF,!has_msi & max_purity_biopsy)$sample_id
      response <- metadata$HMF[match(sample_names, metadata$HMF$sample_id),'response']
      
      data.frame(
         response,
         whitelist=sample_names %in% whitelist_samples
      )
   }
   
   preds$HMF <- cbind(
      preds$HMF,
      getAnnotationsHmf(preds$HMF$sample)
   )
   
   if(!hmf.data.only){
      
      #--------- BRCA-EU ---------#
      getAnnotationsBrcaEu <- function(sample_names){
         #sample_names=preds$BRCA_EU$sample
         out <- metadata$BRCA_EU[match(sample_names, metadata$BRCA_EU$sample),c('response','used_for_evaluation')]
         colnames(out) <- c('response','whitelist')
         return(out)
      }
      
      preds$BRCA_EU <- cbind(
         preds$BRCA_EU,
         getAnnotationsBrcaEu(preds$BRCA_EU$sample)
      )
      
      preds$BRCA_EU.HMF_pipeline <- cbind(
         preds$BRCA_EU.HMF_pipeline,
         getAnnotationsBrcaEu(preds$BRCA_EU.HMF_pipeline$sample)
      )
      
      #--------- PCAWG ---------#
      getAnnotationsPcawg <- function(sample_names){
         #sample_names=preds$PCAWG$sample
         whitelist_samples <- subset(metadata$PCAWG, donor_wgs_exclusion_white_gray=='Whitelist')$sample
         training_samples <- subset(metadata$PCAWG_training, in_training_set)$sample
         
         response <- metadata$PCAWG[match(sample_names, metadata$PCAWG$sample),'response']
         response_umcu <- metadata$PCAWG_training[match(sample_names, metadata$PCAWG_training$sample),'response']
         
         data.frame(
            response,
            response_umcu,
            whitelist=sample_names %in% whitelist_samples,
            in_training_set=sample_names %in% training_samples
         )
      }
      
      preds$PCAWG <- cbind(
         preds$PCAWG,
         getAnnotationsPcawg(preds$PCAWG$sample)
      )
   }
   
   ## Export raw predictions
   preds <- lapply(preds, function(i){
      rownames(i) <- NULL; return(i)
   })
   saveRDS(preds, paste0(out_dir,'/preds.rds'))
   
   #========= Subset preds =========#
   preds_ss <- list(
      HMF = subset(preds$HMF,hr_status!='cannot_be_determined' & hrd_type!='cannot_be_determined' & whitelist),
      HMF_CV = preds$HMF_CV
   )
   
   if(!hmf.data.only){
      preds_ss$BRCA_EU = subset(
         preds$BRCA_EU, 
         hr_status!='cannot_be_determined' 
         & hrd_type!='cannot_be_determined' 
         & whitelist
      )
      preds_ss$BRCA_EU.HMF_pipeline = subset(
         preds$BRCA_EU.HMF_pipeline,
         hr_status!='cannot_be_determined' 
         & hrd_type!='cannot_be_determined' 
         & whitelist
      )
      
      preds_ss$PCAWG = subset(
         preds$PCAWG, 
         hr_status!='cannot_be_determined' 
         & hrd_type!='cannot_be_determined' 
         & whitelist
      )
      
      preds_ss$PCAWG_umcu = subset(
         preds$PCAWG, 
         hr_status!='cannot_be_determined' 
         & hrd_type!='cannot_be_determined' 
         & whitelist & !is.na(response_umcu)
      )
      preds_ss$PCAWG_umcu = within(preds_ss$PCAWG_umcu,{
         response <- response_umcu
         response_umcu <- NULL
      })
      
      preds_ss$PCAWG_umcu_hi_conf_brca <- (function(){
         df <- subset(
            preds$PCAWG, 
            hr_status!='cannot_be_determined' 
            & hrd_type!='cannot_be_determined' 
            & whitelist & !is.na(response_umcu)
            & in_training_set
         )
         df$response <- df$response_umcu
         return(df)
      })()
   }
   
   
   #========= Export supplementary tables =========#
   ## Export xlsx supplementary table
   prob_colnames <- c('BRCA1','BRCA2','none','hrd')
   pcawg_sel_samples <- preds$PCAWG$sample[ !is.na(preds$PCAWG$response_umcu) ]
   
   getHmfIds <- function(v){
      metadata$HMF$hmf_id[ match(v, metadata$HMF$sample_id) ]
   }

   preds_txt <- (function(){
      #--------- Main ---------#
      l <- preds
      l$HMF$sample <- getHmfIds(l$HMF$sample)
      l$HMF_CV$sample <- getHmfIds(l$HMF_CV$sample)
      l$PCAWG <- l$PCAWG[l$PCAWG$sample %in% pcawg_sel_samples,]

      sel_cols <- c('sample','response',prob_colnames,'hr_status','hrd_type','indel_load','indel_rep_load','sv_load')

      prob_colnames <- c('BRCA1','BRCA2','none','hrd')
      counter <- 0
      out <- do.call(rbind, lapply(l, function(i){
         #i=l$HMF_CV
         counter <<- counter + 1
         df <- i[,colnames(i) %in% sel_cols]

         missing_cols <- sel_cols[ !(sel_cols %in% colnames(df)) ]
         for(i in missing_cols){
            df[,i] <- NA
         }
         
         df <- df[,sel_cols]
         
         colnames(df)[colnames(df) %in% prob_colnames] <-
            paste0('p_',colnames(df)[colnames(df) %in% prob_colnames])
         
         df <- cbind(group=names(l)[counter], df)

         return(df)
      }))
      rownames(out) <- NULL

      #--------- Add metadata ---------#
      ##
      hmf_uniq_samples <- subset(metadata$HMF, max_purity_biopsy, hmf_id, drop=T)
      out$is_max_purity_biopsy <- out$sample %in% hmf_uniq_samples
      out$is_max_purity_biopsy[out$group!='HMF'] <- NA
      
      ## cancer type
      out$group <- factor(out$group, unique(out$group))
      out_split <- split(out, out$group)
      
      out_split$BRCA_EU$cancer_type <- 'Breast'
      out_split$BRCA_EU.HMF_pipeline$cancer_type <- 'Breast'
      
      out_split$HMF$cancer_type <- metadata$HMF$primary_tumor_location[ match(out_split$HMF$sample, metadata$HMF$hmf_id) ]
      out_split$HMF_CV$cancer_type <- metadata$HMF$primary_tumor_location[ match(out_split$HMF_CV$sample, metadata$HMF$hmf_id) ]
      
      out_split$PCAWG$cancer_type <- metadata$PCAWG$cancer_type[ match(out_split$PCAWG$sample, metadata$PCAWG$sample) ]
      
      ## is for perf eval
      out_split$HMF$used_for_perf_eval <- NA
      out_split$HMF_CV$used_for_perf_eval <- TRUE
      
      out_split$BRCA_EU$used_for_perf_eval <- out_split$BRCA_EU$sample %in% preds_ss$BRCA_EU$sample
      out_split$BRCA_EU.HMF_pipeline$used_for_perf_eval <- out_split$BRCA_EU.HMF_pipeline$sample %in% preds_ss$BRCA_EU.HMF_pipeline$sample
      out_split$PCAWG$used_for_perf_eval <- out_split$PCAWG$sample %in% preds_ss$PCAWG_umcu_hi_conf_brca$sample
      
      ## umcu response for PCAWG
      out_split$PCAWG$response <- preds$PCAWG$response_umcu[ match(out_split$PCAWG$sample, preds$PCAWG$sample) ]
      
      out <- do.call(rbind, out_split); rownames(out) <- NULL
      
      ##
      out$in_training_set <- metadata$HMF$in_training_set[ match(out$sample, metadata$HMF$hmf_id) ]
      
      out$used_for_pancancer_analysis <- with(out,{
         v <- ifelse(hr_status=='cannot_be_determined' | hrd_type=='cannot_be_determined',FALSE,TRUE)
         v[!(group %in% c('HMF','PCAWG'))] <- NA
         v[group=='HMF' & !(sample %in% hmf_uniq_samples)] <- FALSE
         return(v)
      })
      
      return(out)
   })()

   write.table(preds_txt, paste0(out_dir,'/preds.txt'), sep='\t', quote=F, row.names=F)

   ## Export features
   features <- (function(){
      raw <- lapply(contexts, function(i){
         do.call(cbind, unname(i[c('snv','indel','sv')]))
      })
      #names(raw) <- paste0(names(raw),'.raw')

      proc <- lapply(contexts, rf$trans.func)
      #rownames(proc$HMF) <- metadata$HMF$hmf_id[ match(rownames(proc$HMF), metadata$HMF$sample_id) ]

      meltCustom <- function(l){
         #l=raw
         rownames(l$HMF) <- metadata$HMF$hmf_id[ match(rownames(l$HMF), metadata$HMF$sample_id) ]
         l$PCAWG <- l$PCAWG[rownames(l$PCAWG) %in% pcawg_sel_samples,]

         do.call(rbind,lapply(names(l), function(i){
            #i=names(l)[[1]]
            df <- as.data.frame(l[[i]])
            df <- cbind(
               group=i,
               sample=rownames(df),
               df
            )
            rownames(df) <- NULL
            return(df)
         }))
      }

      list(
         raw=meltCustom(raw),
         proc=meltCustom(proc)
      )

   })()

   write.xlsx(features, paste0(out_dir,'/features.xlsx'))
   
   ####################################################################################################
   # Sorted prediction probs                                                                          #
   ####################################################################################################
   message('Plotting sorted prediction probs...')
   
   plotPredProbs <- function(
      pred=NULL, probs=NULL, response=NULL, sample.names=NA,
      error=NA, error.as.linerange=F,
      
      class.colors=c(BRCA1='#f58225', BRCA2='#69439d'),
      class.labels=c(BRCA1='BRCA1-type HRD',BRCA2='BRCA2-type HRD'),
      response.colors=c(BRCA1='#f58225', BRCA2='#69439d','none'=NA),
      
      default.class.name='HRD',
      
      cutoff.hline.y=0.5, show.confusion=0, confusion.xpos=NULL,
      rel.heights=c(1, 0.3), ylims=c(0,1), title=NULL, y.lab=NULL, top=NA, bottom=NA, 
      
      hide.info=c(), hide.response=F, do.ordering=T
   ){
      
      # probs=pred_icgc[c('BRCA1','BRCA2')]
      # colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
      # response <- as.factor(pred_icgc$response)
      # sample.names=pred_icgc$sample
      
      # error=smooth.spline(cv100_pred_stats$hrd.sd)$y
      # class.colors=c('BRCA1'='#f9b47c', 'BRCA2.mean'='#a58ec4')
      
      #pred=preds_ss$HMF
      
      #--------- Pre-process ---------#
      if(!is.null(data)){
         response <- pred$response
         probs <- pred[c('BRCA1','BRCA2')]
      }
      
      if(!(is.data.frame(probs) || is.matrix(probs))){
         stop('`probs` must be a data.frame or matrix')
      }
      
      if(is.matrix(probs)){ probs <- as.data.frame(probs) }
      df <- cbind(sum=rowSums(probs), probs, response, error, sample.names)
      
      if(do.ordering){
         df <- df[do.call(order, df),]
      }
      
      #df$sum <- NULL
      df$index <- 1:nrow(df)
      
      df_full <- df
      if(!is.na(top)){ 
         df <- tail(df, top) 
      } else if(!is.na(bottom)){ 
         df <- head(df, bottom) 
      }
      
      df_melt <- melt(subset(df, select=-sum),c('index','response','error','sample.names'))
      colnames(df_melt)[5:6] <- c('class','probs')
      
      #--------- Prediction probs ---------#
      plot_pred <- ggplot() +
         geom_bar(
            data=df_melt, mapping=aes(index, probs, fill=class),
            stat='identity', position='stack', width=1
         ) +
         scale_x_continuous(expand=c(0,0)) +
         scale_y_continuous(limits=ylims) +
         labs(x=paste0(nrow(df_full),' samples'), y='Probability', fill='Prediction class') +
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            axis.title.y=element_text(size=12),
            axis.title.x=element_blank(),
            #axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            legend.justification=c(0,0.5),
            #legend.key=element_rect(size=1.5, color='white'),
            legend.key.size=unit(11,'pt'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=10, face='bold'),
            plot.margin = unit(c(2,0,6,2),'pt')
         )
      
      if(!is.null(class.colors)){
         plot_pred <- plot_pred + scale_fill_manual(values=class.colors, labels=class.labels)
      }
      
      if(!anyNA(error)){
         pd_error <- data.frame(
            index=df$index,
            y=rowSums(df[colnames(df) %in% colnames(probs)])
         )
         
         pd_error$upper <- pd_error$y + df$error
         pd_error$lower <- pd_error$y - df$error
         
         if(error.as.linerange){
            plot_pred <- plot_pred + geom_linerange(data=pd_error, aes(index, ymin=lower, ymax=upper), size=0.5, color='black')
         } else {
            pd_error <- melt(pd_error,c('index','y'))
            plot_pred <- plot_pred + geom_line(data=pd_error, aes(index, value, group=variable), size=0.5, color='black')}
      }
      
      if(!is.null(cutoff.hline.y)){
         plot_pred <- plot_pred + geom_hline(yintercept=cutoff.hline.y, linetype='dashed', size=0.3)
      }
      
      if(!is.null(title)){ plot_pred <- plot_pred + ggtitle(title) }
      if(!is.null(y.lab)){ plot_pred <- plot_pred + ylab(y.lab) }
      
      if(show.confusion>=1){
         
         m_confusion <- confusionMatrix(
            predicted=df_full$sum, 
            #actual=toBinaryResponse(df_full$response,c('BRCA1','BRCA2'),1,'none',0),
            actual=df_full$response %in% c('BRCA1','BRCA2'),
            cutoff=cutoff.hline.y
         )
         m_confusion <- rbind(
            m_confusion[c('tp','fp')],
            m_confusion[c('fn','tn')]
         )
         
         rownames(m_confusion) <- NULL
         colnames(m_confusion) <- NULL
         
         if(show.confusion==1){
            cols <- matrix(c('red','grey50'), nrow(m_confusion), ncol(m_confusion), byrow = T)
         } else if(show.confusion==2){
            m_confusion <- as.matrix(m_confusion[,1] + m_confusion[,2])
            cols <- as.matrix(c('#C00000','#548235'))
         }
         
         tt <- ttheme_minimal(core=list(fg_params = list(col = cols)), base_size=10)
         table_confusion <- tableGrob(m_confusion, theme = tt)
         
         if(is.null(confusion.xpos)){
            table_xpos <- nrow(df_full) * 0.1
         } else {
            table_xpos <- confusion.xpos
         }
         
         #plot + annotation_custom(table_confusion, xmin=table_xpos, xmax=table_xpos, ymin=cutoff, ymax=cutoff)
         plot_pred <- plot_pred + 
            annotation_custom(
               table_confusion, xmin=table_xpos, xmax=table_xpos, 
               ymin=cutoff.hline.y, ymax=cutoff.hline.y
            )
      }
      
      #--------- Response ---------#
      df_response <- df_melt[,c('index','response','sample.names')]
      df_response <- df_response[!duplicated(df_response),]
      
      df_response$sample.names <- factor(df_response$sample.names, unique(df_response$sample.names))
      neg_class <- names(response.colors)[is.na(response.colors)]
      levels(df_response$response)[levels(df_response$response) %in% neg_class] <- NA
      # response.colors2 <- response.colors
      # names(response.colors2)[response.colors2 %in% neg_class] <- ''
      
      if(!is.na(sample.names)){
         plot_response <- ggplot(df_response, aes(sample.names, paste0('n = ',nrow(probs)), fill=response))
      } else {
         plot_response <- ggplot(df_response, aes(index, paste0('n = ',nrow(probs)), fill=response)) +
            scale_x_continuous(expand=c(0,0))
      }
      
      plot_response <- plot_response +
         geom_tile() +
         
         scale_y_discrete(expand=c(0,0)) +
         labs(x=paste0(nrow(df_full),' samples'), fill='Gene deficiency') +
         
         theme_bw() +
         theme(
            axis.line=element_blank(),
            #axis.text.y=element_text(face='bold',size=10),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=12),
            panel.border=element_rect(fill=NA, color='black', size=0.5),
            panel.grid=element_blank(),
            legend.justification=c(0,0.5),
            legend.key=element_rect(size=1.2, color='white'),
            legend.key.size=unit(10,'pt'),
            legend.text=element_text(size=9),
            legend.title=element_text(size=10, face='bold'),
            plot.margin = unit(c(2,0,2,2),'pt')
         )
      
      if(!is.na(sample.names)){
         plot_response <- plot_response + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
         plot_pred <- plot_pred + theme(axis.ticks.x=element_blank())
      }
      
      if(!is.null(response.colors)){
         plot_response <- plot_response + 
            scale_fill_manual(values=response.colors, breaks=names(response.colors)[!is.na(response.colors)])
      }
      
      #--------- Hide info ---------#
      if('left' %in% hide.info){
         theme_hidden_left <- theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
         )
         
         plot_pred <- plot_pred + theme_hidden_left
         plot_response <- plot_response + theme_hidden_left
      }
      
      if('right' %in% hide.info){
         theme_hidden_right <- theme(legend.position='none')
         plot_pred <- plot_pred + theme_hidden_right
         plot_response <- plot_response + theme_hidden_right
      }
      
      if('bottom' %in% hide.info){
         theme_hidden_bottom <- theme(axis.title.x=element_blank())
         plot_pred <- plot_pred + theme_hidden_bottom
         plot_response <- plot_response + theme_hidden_bottom
      }
      
      #--------- Combine ---------#
      if(hide.response){
         plot_pred + theme(axis.title.x=element_text(), axis.text.x=element_text())
      } else {
         plot_grid(plot_pred, plot_response, align='v', axis='tblr', ncol=1, rel_heights=rel.heights)
      }
   }
   
   plotPredProbsSplit <- function(
      pred, top=300, title='CHORD on HMF dataset', 
      show.confusion=0, confusion.xpos=600, cutoff.hline.y=NULL, hide.response=F, axis.title.x=NULL
   ){
      #pred=preds_ss$HMF
      
      #probs <- pred[c('BRCA1','BRCA2')]
      #colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
      
      n_samples <- nrow(pred)
      bottom <- n_samples - top
      
      plots <- list()
      
      plots$left <- plotPredProbs(
         pred, title=paste0('Bottom ', bottom), bottom=bottom, 
         hide.info=c('right','bottom'), rel.heights=c(1, 0.25), 
         show.confusion=show.confusion, confusion.xpos=confusion.xpos, cutoff.hline.y=cutoff.hline.y,
         hide.response=hide.response
      )
      
      plots$right <- plotPredProbs(
         pred, title=paste0('Top ', top), top=top, 
         hide.info=c('left','bottom'), rel.heights=c(1, 0.25), cutoff.hline.y=cutoff.hline.y,
         hide.response=hide.response
      )
      
      if(hide.response){
         plots <- lapply(plots, function(i){ i + theme(axis.title.x=element_blank()) })
      }
      
      if(is.null(axis.title.x)){
         axis.title.x <- sprintf('%s samples            ', n_samples)
      }
      
      out <- arrangeGrob(plots$left, plots$right, nrow=1, widths=c(0.55, 0.45), top=title, bottom=axis.title.x)
      class(out) <- c(class(out),'arrangedGrob')
      return(out)
   }
   
   lp_pred_probs <- list()
   
   lp_pred_probs$HMF = plotPredProbsSplit(
      preds_ss$HMF, title='HMF', 
      cutoff=0.5, show.confusion=2, axis.title.x=paste0(nrow(preds_ss$HMF),' HMF patients'),
      hide.response=F
   )
   lp_pred_probs$HMF_CV = plotPredProbsSplit(preds$HMF_CV, title='HMF CV', cutoff.hline.y=NULL)
   
   if(!hmf.data.only){
      lp_pred_probs$BRCA_EU = plotPredProbs(preds_ss$BRCA_EU, title='BRCA EU', cutoff.hline.y=NULL)
      lp_pred_probs$BRCA_EU.full = plotPredProbs(preds$BRCA_EU, title='BRCA EU, full', show.confusion=2)
      lp_pred_probs$BRCA_EU.HMF_pipeline = plotPredProbs(preds_ss$BRCA_EU.HMF_pipeline, title='BRCA EU, HMF pipeline')
      #lp_pred_probs$PCAWG = plotPredProbsSplit(preds_ss$PCAWG, title='PCAWG')
      
      lp_pred_probs$PCAWG_umcu = plotPredProbsSplit(
         preds_ss$PCAWG_umcu, title='PCAWG (UMCU gene ann)', 
         top=150, cutoff.hline.y=0.5, show.confusion=2, confusion.xpos=150
      )
      
      lp_pred_probs$PCAWG_hi_conf_brca = plotPredProbsSplit(
         preds_ss$PCAWG_umcu_hi_conf_brca, 
         title='PCAWG (UMCU gene ann, high confidence BRCA status)', top=100,
         cutoff.hline.y=NULL
      )
      
      
      lp_pred_probs$HMF_PCAWG_umcu = (function(){
         common_cols <- intersect(colnames(preds_ss$HMF), colnames(preds_ss$PCAWG_umcu))
         df <- rbind(preds_ss$HMF[,common_cols], preds_ss$PCAWG_umcu[,common_cols])
         plotPredProbsSplit(
            df, title='HMF + PCAWG (UMCU gene ann)', 
            top=400, cutoff.hline.y=0.5, show.confusion=2, confusion.xpos=500, hide.response=T
         )
      })()
   }

   pdf(paste0(out_dir,'/sorted_probs.pdf'), 8,4)
   for(i in lp_pred_probs){
      if('arrangedGrob' %in% class(i)){ grid.newpage() }
      grid.draw(i)
   }
   dev.off()
   
   
   ####################################################################################################
   # VCF vs BAM predictions                                                                           #
   ####################################################################################################
   message('Plotting BRCA EU bam vs vcf predictions...')
   
   plotPredsVcfsVsBams <- function(pred_bam, pred_vcf){
      #pred_vcf=preds_ss$BRCA_EU
      #pred_bam=preds_ss$BRCA_EU.HMF_pipeline
      
      l <- list(pred_bam=pred_bam, pred_vcf=pred_vcf)
      common_samples <- Reduce(intersect,lapply(l,`[[`,'sample'))
      common_samples <- l$pred_bam$sample[order(l[[1]]$hrd)] ## Order by bam predictions
      
      l <- lapply(l, function(i){ i[match(common_samples, i$sample),] })
      
      brca_colors <- c(BRCA1='#f58225',BRCA2='#69439d',none=NA)
      
      #========= Sorted probs =========#
      plots <- lapply(l, function(i){
         plotPredProbs(i, cutoff.hline.y=NULL, hide.response=T, do.ordering=F) + 
            theme(
               axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank()
            )
      })
      
      plots$pred_bam <- plots$pred_bam + 
         geom_label(aes(x=1, y=1, hjust=0, vjust=1, label='Variants called from BAM\nfiles with HMF pipeline'), fill='white')
      
      plots$pred_vcf <- plots$pred_vcf +
         geom_label(aes(x=1, y=1, hjust=0, vjust=1, label='Variants downloaded from ICGC'), fill='white')
      
      #========= HRD score difference =========#
      df_diff <- forceDfOrder(data.frame(
         sample=l$pred_bam$sample,
         diff=l$pred_bam$hrd - l$pred_vcf$hrd
      ))
      
      plots$diff <- ggplot(df_diff, aes(sample,diff)) + 
         geom_bar(stat='identity') +
         geom_hline(yintercept=0, size=0.25, color='red') +
         ylab('Difference') +
         theme_bw() +
         theme(
            panel.grid.minor.y=element_blank(),
            #panel.grid.major.y=element_blank(),
            #axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
         )
      
      #========= Response =========#
      df_response <- forceDfOrder(l[[1]][c('sample','response')]) 
      plots$response <- ggplot(df_response, aes(x=sample,y=1, fill=response)) +
         geom_tile() +
         scale_fill_manual(name='Gene deficiency', values=brca_colors, breaks=c('BRCA1','BRCA2')) +
         scale_x_discrete(expand=c(0,0), name=paste0('\n', nrow(df_response),' ICGC (BRCA-EU) samples')) +
         scale_y_discrete(expand=c(0,0)) +
         theme_bw() +
         theme(
            panel.border=element_rect(fill=NA),
            panel.grid.major=element_blank(),
            axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
            axis.title.y=element_blank()
         )
      
      plots_rearranged <- plots[c('pred_bam','diff','pred_vcf','response')]
      
      plot_grid(plotlist=plots_rearranged, align='v',axis='lr',ncol=1, rel_heights=c(1,0.5,1,0.8))
   }
   
   if(!hmf.data.only){
      pdf(paste0(out_dir,'/brca_eu_bam_vs_vcf_preds.pdf'), 10, 6)
      grid.draw(plotPredsVcfsVsBams(preds_ss$BRCA_EU.HMF_pipeline, preds_ss$BRCA_EU))
      dev.off()
   }
   
   ####################################################################################################
   # Performance                                                                                      #
   ####################################################################################################
   message('Plotting performance...')
   
   #========= Confusion =========#
   confusionMatrixMcCustom <- function(df){
      #df=preds$BRCA_EU
      
      confusion <- list()
      
      confusion <- confusionMatrix(
         df[c('BRCA1','BRCA2','none')],
         oneHotEncode(as.factor(df$response)),
         cutoff.interval=0.002
      )

      confusion$hrd <- confusionMatrix( 
         df$hrd, 
         ifelse(df$response=='none',FALSE,TRUE),
         cutoff.interval=0.002
      )
      
      confusion <- confusion[c('hrd','BRCA1','BRCA2')]
      names(confusion) <- c('HRD','BRCA1-type HRD','BRCA2-type HRD')
      
      confusion <- confusion[sapply(confusion,function(i){ !is.null(i) })]
      
      return(confusion)
   }
   
   #========= Compound performance =========#
   plotPerfCompoundCustom <- function(
      confusion, compound.metric='pr', plot.title=NULL,
      #line.colors=c('black','#F8766D','#619CFF'),
      line.colors=c('#7e7e7e','#f58225','#69439d'),
      line.size=c(1, 0.5, 0.5),
      legend.position=c(0.5, 0.1)
   ){
      
      l <- lapply(confusion, function(i){
         #i=confusion$HRD
         m <- calcPerfCompound(i, compound.metric, metric.names.as.x.y=T)
         m <- as.data.frame(m)
         m$auc <- calcAUC(m$x,m$y)
         return(m)
      })
      df <- as.data.frame(do.call(rbind,l))
      rownames(df) <- NULL
      
      ##
      df$class <- unlist(lapply(names(l), function(i){ rep(i, nrow(l[[i]])) }))
      df$class <- sprintf('%s (%s)', df$class, format(round(df$auc,3),nsmall=3))
      df$class <- factor(df$class,unique(df$class))
      
      ##
      axis_titles <- switch(
         compound.metric,
         roc=c('False positive rate','True positive rate'),
         pr=c('Recall (TPR)','Precision (PPV)'),
         npv_tnr=c('True negative rate','Negative predictive value')
      )
      
      plot <- ggplot(df, aes(x,y, color=class, size=class)) + 
         geom_path() +
         labs(x=axis_titles[1],y=axis_titles[2], color='Prediction class (AUC)') +
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=10),
            legend.title=element_text(size=10),
            legend.text=element_text(size=9),
            legend.key.height=unit(10,'pt'),
            legend.key.width=unit(12,'pt'),
            legend.position=legend.position,
            legend.justification=c(0.5, 0),
            legend.background=element_rect(color='black',size=0.2,fill=alpha('white',0.7))
         )
      
      if(!is.null(line.colors)){ plot <- plot + scale_color_manual(values=line.colors) }
      if(!is.null(line.size)){ plot <- plot + scale_size_manual(values=line.size, guide=F) }
      if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
      
      return(plot)
   }
   
   #========= Simple performance =========#
   plotPerfCustom <- function(
      confusion, metrics=c('tpr','fpr'), plot.title=NULL,
      line.colors=NULL, line.types=NULL,
      legend.position=c(0.5, 0.35), na.rough.fix=F, show.cutoff=F
   ){
      # confusion=confusions$HMF_CV
      # metrics=c('tpr','fpr')
      # metrics='prec'
      
      l <- lapply(names(confusion), function(i){
         #i='BRCA1-like HRD'
         df <- calcPerf(confusion[[i]], metrics=metrics, melt=T, add.start.end.values=F)
         if(na.rough.fix){ df[is.na(df)] <- 0 }
         df <- unique(df)
         df$class <- i
         return(df)
      })
      df <- as.data.frame(do.call(rbind,l))
      rownames(df) <- NULL
      
      ## Use human readable metric names 
      metric_names <- c(
         tpr='True positive rate',
         tnr='True negative rate',
         fpr='False positive rate',
         fnr='False negative',
         prec='Precision',
         f1='F1-score',
         mcc=' MCC'
      )
      
      for(i in names(metric_names)){
         df$metric <- gsub(i, metric_names[i], df$metric, fixed=T)
      }
      
      ##
      for(i in 1:ncol(df)){
         if(is.character(df[[i]])){ df[[i]] <- factor(df[[i]], unique(df[[i]])) }
      }
      
      ##
      p <- ggplot(df, aes(cutoff, value, color=metric, linetype=class)) + 
         geom_path() +
         ylim(0,1) +
         labs(x='Classif. cutoff',y='Metric value',color='Metric',linetype='Prediction class') +
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=10),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8),
            legend.position=legend.position,
            legend.direction='vertical',
            legend.key.height=unit(10,'pt'),
            legend.key.width=unit(12,'pt'),
            legend.background=element_rect(fill=NA),
            legend.spacing.y=unit(1,'pt'),
            legend.box.background=element_rect(color='black',size=0.2,fill=alpha('white',0.7))
         )
      
      if(!is.null(line.colors)){ p <- p + scale_color_manual(values=line.colors) }
      if(!is.null(line.types)){ p <- p + scale_linetype_manual(values=line.types) }
      if(!is.null(plot.title)){ p <- p + ggtitle(plot.title) }
      
      
      ## Calculate optimal cutoff
      if(show.cutoff){
         if(all(c('tpr','fpr') %in% metrics)){
            m <- calcPerf(confusion$HRD, metrics)
            # m <- cbind(
            #    m, diff=m[,'tpr'] - m[,'fpr']
            # )
            max_index <- which.max(m[,'tpr'] - m[,'fpr'])
            cutoff <- m[max_index,'cutoff']
         } else if('f1' %in% metrics){
            m <- calcPerf(confusion$HRD, 'f1')
            cutoff <- m[which.max(m[,'f1']),'cutoff']
         } else {
            warning('Cutoff calculation not available for the chosen metrics')
         }
         
         p <- p + geom_vline(xintercept=cutoff, linetype='dotted')
      }
      
      return(p)
   }
   
   message('> Making confusion matrices...')
   confusions <- lapply(
      preds_ss[names(preds_ss)!='HMF'], ## Exclude whole HMF
      confusionMatrixMcCustom
   )
   saveRDS(confusions, paste0(out_dir,'/confusions.rds'))
   
   message('> Making plots...')
   lp_perf <- lapply(names(confusions), function(i){
      #i=confusions$HMF_CV
      m <- confusions[[i]]
      
      grobs <- list(
         plotPerfCompoundCustom(m, compound.metric='pr'),
         plotPerfCompoundCustom(m, compound.metric='roc'),
         ggplot() + theme_void(),
         
         plotPerfCustom(m, c('tpr','fpr')),
         plotPerfCustom(m, 'prec'),
         plotPerfCustom(m, 'f1')
      )
      grobs <- lapply(grobs, function(i){
         i + theme(plot.margin=ggplot2::margin(8,8,8,8))
      })
      
      arrangeGrob(
         grobs=grobs,
         ncol=3, top=i, padding = unit(0.2,'line')
      )
   })
   
   pdf(paste0(out_dir,'/performance.pdf'), 10, 6.5)
   for(i in lp_perf){
      grid.newpage()
      grid.draw(i)
   }
   dev.off()
   
   ####################################################################################################
   # CHORD vs HRDetect                                                                                #
   ####################################################################################################
   
   mkPredPairs <- function(pred, metadata){
      df1 <- pred[,c('sample','hrd','hr_status')]
      colnames(df1)[2:3] <- c('p_chord','passed_chord_qc')
      df1$passed_chord_qc <- ifelse(df1$passed_chord_qc!='cannot_be_determined',TRUE,FALSE)
      
      df2 <- metadata[,c('sample','p_hrdetect','response')]
      
      out <- merge(df1,df2,'sample',all=F)
      out <- subset(out, !is.na(p_chord) & !is.na(p_hrdetect) & passed_chord_qc)
      return(out)
   }
   
   pred_pairs <- list()
   
   pred_pairs$BRCA_EU <- mkPredPairs(preds$BRCA_EU, metadata$BRCA_EU)
   pred_pairs$PCAWG <- mkPredPairs(preds$PCAWG, metadata$PCAWG)
   
   pred_pairs$PCAWG_umcu <- mkPredPairs(preds_ss$PCAWG_umcu, metadata$PCAWG_training)
   pred_pairs$PCAWG_umcu_hi_conf_brca <- mkPredPairs(preds_ss$PCAWG_umcu_hi_conf_brca, metadata$PCAWG_training)
   
   ## Subset for samples with UMCU annotation for fair comparison
   pred_pairs$PCAWG <- pred_pairs$PCAWG[pred_pairs$PCAWG$sample %in% pred_pairs$PCAWG_umcu$sample,]
   
   
   ## Which samples are in both BRCA EU and PCAWG   
   # metadata_pcawg_ss <- metadata$PCAWG[,c('sample','submitter_sample_id')]
   # metadata_pcawg_ss <- metadata_pcawg_ss[grep('^PD',metadata_pcawg_ss$submitter_sample_id),]
   # metadata_pcawg_ss$in_brca_eu_dataset <- metadata_pcawg_ss$submitter_sample_id %in% paste0(metadata$BRCA_EU$sample,'a')
   # pcawg_sample_blacklist <- metadata_pcawg_ss$sample[metadata_pcawg_ss$in_brca_eu_dataset]
   # 
   # pred_pairs$PCAWG$in_brca_eu_dataset <- pred_pairs$PCAWG$sample %in% pcawg_sample_blacklist
   # pred_pairs$PCAWG$brca_eu_id <- metadata_pcawg_ss[match(pred_pairs$PCAWG$sample, metadata_pcawg_ss$sample),'submitter_sample_id']
   # pred_pairs$PCAWG$brca_eu_id[!pred_pairs$PCAWG$in_brca_eu_dataset] <- NA
   
   #dim(subset(pred_pairs$PCAWG, passed_chord_qc & !in_brca_eu_dataset))
   #dim(subset(pred_pairs$PCAWG, passed_chord_qc))
   
   plotChordVsHrdetect <- function(
      df, title=NULL,
      class.colors=c(BRCA1='#f58225', BRCA2='#69439d',none='grey'),
      chord.cutoff=0.5, hrdetect.cutoff=0.7
   ){
      # df=pred_pairs$PCAWG
      # title=NULL
      # class.colors=c(BRCA1='#f58225', BRCA2='#69439d',none='grey')
      # chord.cutoff=0.5
      # hrdetect.cutoff=0.7
      
      quadrant_counts <- with(df,{
         l <- list(
            nw = p_hrdetect < hrdetect.cutoff & p_chord >= chord.cutoff,
            ne = p_hrdetect >= hrdetect.cutoff & p_chord >= chord.cutoff,
            
            sw = p_hrdetect < hrdetect.cutoff & p_chord < chord.cutoff,
            se = p_hrdetect >= hrdetect.cutoff & p_chord < chord.cutoff
         )
         v <- unlist(lapply(l, sum))
         matrix(v,nrow=2,ncol=2,byrow=T)
      })
      
      p <- ggplot(df, aes(y=p_chord, x=p_hrdetect, color=response)) + 
         
         geom_hline(yintercept=chord.cutoff, linetype='dotted') +
         geom_vline(xintercept=hrdetect.cutoff, linetype='dotted') +
         
         geom_point(data=subset(df, response=='none')) +
         geom_point(data=subset(df, response!='none')) +
         scale_color_manual(values=class.colors, breaks=c('BRCA1','BRCA2')) +
         
         xlim(0,1) + ylim(0,1) +
         labs(color='Gene deficiency',y='P(HRD) CHORD', x='P(HRD) HRDetect') +
         
         draw_grob(
            tableGrob(quadrant_counts, theme=ttheme_minimal(base_colour='darkred')),
            x=hrdetect.cutoff, y=chord.cutoff, hjust=0.5, vjust=0.5,
         ) +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5)
         )
      
      if(!is.null(title)){ p <- p + ggtitle(title) }
      
      return(p)
   }
   
   # lp_pred_pairs <- list(
   #    plotChordVsHrdetect(pred_pairs$PCAWG, title='PCAWG'),
   #    plotChordVsHrdetect(pred_pairs$PCAWG_hi_conf_brca, title='PCAWG (high confidence BRCA status)'),
   #    plotChordVsHrdetect(pred_pairs$BRCA_EU, title='BRCA-EU')
   # )
   
   lp_pred_pairs <- lapply(names(pred_pairs), function(i){
      plotChordVsHrdetect(pred_pairs[[i]], title=i)
   })
   
   pdf(paste0(out_dir,'/chord_vs_hrdetect.pdf'),6,5)
   for(i in lp_pred_pairs){ plot(i) }
   dev.off()
   
   write.xlsx(pred_pairs, paste0(out_dir,'/preds_chord_vs_hrdetect.xlsx'))
   
   
   ####################################################################################################
   # Feature importance                                                                               #
   ####################################################################################################
   message('Plotting feature importance...')
   
   #========= Functions =========#
   getFeatImp <- function(model, class=NULL, f_feature_rename=NULL){
      df <- as.data.frame(randomForest::importance(rf, class=class, type=1))
      colnames(df) <- 'mda'
      
      feature_names <- rownames(df)
      if(!is.null(f_feature_rename)){
         feature_names <- f_feature_rename(feature_names)
      }
      
      df <- cbind(feature=feature_names,df); rownames(df) <- NULL
      df[order(df$mda, decreasing=F),]
   }
   
   getFeatImpCv <- function(cv.models, class=NULL, melt.df=T, f_feature_rename=NULL){
      #cv.models=cv_rfs
      
      l <- lapply(cv.models, function(i){
         randomForest::importance(i, class=class, type=1)[,1]
      })
      
      ## Fill in NA in folds with missing features vs all features across folds
      all_features <- unique(unlist(lapply(l, names), use.names=F))
      l <- lapply(l, function(i){
         #i=l[[1]]
         
         missing_features <- all_features[!(all_features %in% names(i))]
         i[missing_features] <- NA
         
         ## Force consistent order
         return(i[all_features])
      })
      
      df <- as.data.frame(do.call(cbind,l))
      colnames(df) <- 1:length(l)
      df <- df[order(rowMeans(df, na.rm=T), decreasing=F),]
      
      feature_names <- rownames(df)
      if(!is.null(f_feature_rename)){
         feature_names <- f_feature_rename(feature_names)
      }
      
      df <- cbind(feature=feature_names, df); rownames(df) <- NULL
      
      if(!melt.df){
         return(df)
      }
      
      df_melt <- suppressMessages({ melt(df, by='feature') })
      colnames(df_melt)[2:3] <- c('fold','mda')
      return(df_melt)
   }
   
   plotFeatImp <- function(imp, imp.cv, plot.title=NULL, hide.axis.title.x=F, seed=1){
      #imp=l_imp[[1]]$pred
      #imp.cv=l_imp[[1]]$cv
      
      set.seed(seed)
      
      imp.cv$feature <- factor(imp.cv$feature, unique(imp.cv$feature))
      imp$feature <- factor(imp$feature, unique(imp.cv$feature))
      
      ## Note: some feature importances will be NA for certain CV folds (i.e. do to feature selection)
      ## ggplot will remove these geom_points automatically and return a warning
      plot <- ggplot(imp.cv, aes(feature, mda)) + 
         geom_boxplot(width=0.8, color='grey', size=0.3, outlier.color=NA) +
         geom_jitter(width=0.1, color='black', size=1, shape=1) +
         geom_point(data=imp, aes(feature, mda), shape=108, size=1, stroke=7, color='red') +
         scale_y_continuous(limits=c(0,NA)) +
         labs(y='Feature importance\n(Mean decrease in accuracy)') +
         
         coord_flip() +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            panel.grid.major.y=element_line(linetype='dashed'),
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=9),
            axis.text.y=element_text(size=9)
         )
      
      if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
      if(hide.axis.title.x){ plot <- plot + theme(axis.title.x=element_blank()) }
      
      return(plot)
   }
   
   #========= Main =========#
   rf_classes <- c('HRD','BRCA2','BRCA1')
   l_imp <- lapply(rf_classes, function(i){
      #print(i)
      if(i=='HRD'){
         rf_class <- NULL
      } else {
         rf_class <- i
      }
      
      imp <- list()
      imp$pred <- getFeatImp(rf, class=rf_class)
      imp$cv <- getFeatImpCv(cv_rfs, class=rf_class)
      
      return(imp)
   })
   names(l_imp) <- rf_classes
   
   ## 3x1
   p_imp <- plot_grid(
      plotFeatImp(l_imp[[1]]$pred, l_imp[[1]]$cv, plot.title='HRD', hide.axis.title.x=T),
      plotFeatImp(l_imp[[2]]$pred, l_imp[[2]]$cv, plot.title='BRCA2-type HRD'),
      plotFeatImp(l_imp[[3]]$pred, l_imp[[3]]$cv, plot.title='BRCA1-type HRD', hide.axis.title.x=T),
      
      nrow=1, align='hv', axis='tblr'
   )
   
   ## Calculate pdf height
   imp_merged <- do.call(rbind, lapply(l_imp, function(i){
      rbind(i$pred, i$cv[c('feature','mda')])
   }))
   n_features <- length(unique(imp_merged$feature))
   
   pdf(paste0(out_dir,'/feature_importance.pdf'), 10, n_features*0.2+0.4) ## Formula takes into account plot/axis title heights
   grid.draw(p_imp)
   dev.off()
   
   
   ####################################################################################################
   # CV feature selection                                                                             #
   ####################################################################################################
   message('Plotting repeated CV feature selection summary...')
   
   plotTimesHrd <- function(
      sel_sample_cv_summary,
      response.colors=c('BRCA1'='#f58225', 'BRCA2'='#69439d','none'='grey'),
      
      hline.upper.yintercept=60, hline.upper.label="Below: 'BRCA1' and 'BRCA2' samples filtered",
      hline.lower.yintercept=40, hline.lower.label="Above: 'none' samples filtered",
      
      rel.heights=c(1, 0.2),
      show.top=200
   ){
      ## Prep data
      df <- sel_sample_cv_summary
      df$response <- factor(df$response, c('none','BRCA1','BRCA2'))
      df <- df[do.call(order, df[c('n_times_pos','response')]),]
      df$index <- 1:nrow(df)
      
      n_samples <- nrow(df)
      
      ## Make summary table
      counts_table <- aggregate(df$sample, by=list(df$response), FUN=length)
      colnames(counts_table) <- c('Class','Total')
      counts_table$Blacklist <- aggregate(!df$whitelist, by=list(df$response), FUN=sum)$x
      counts_table$Whitelist <- aggregate(df$whitelist, by=list(df$response), FUN=sum)$x
      
      rownames(counts_table) <- counts_table$Class
      counts_table$Class <- NULL
      
      counts_table <- rbind(
         counts_table,
         Sum=colSums(counts_table)
      )
      counts_table <- cbind(Class=rownames(counts_table), counts_table); rownames(counts_table) <- NULL
      
      if(!is.na(show.top)){ df <- tail(df, show.top) }
      
      ## Main plot
      p <- ggplot(df, aes(index, n_times_pos, color=response)) + 
         geom_point() +
         
         #scale_x_continuous(minor_breaks=df$index) +
         scale_y_continuous(breaks=seq(0, max(df$n_times_pos), 10)) +
         scale_color_manual(values=response.colors, breaks=names(response.colors)) +
         
         geom_hline(yintercept=hline.upper.yintercept, linetype='dotted') +
         annotate(
            'text', x=min(df$index), y=hline.upper.yintercept, hjust=0, vjust=2, size=3,
            label=hline.upper.label
         ) +
         
         geom_hline(yintercept=hline.lower.yintercept, linetype='dotted') +
         annotate(
            'text', x=min(df$index), y=hline.lower.yintercept, hjust=0, vjust=-1, size=3,
            label=hline.lower.label
         ) +
         
         labs(
            #y='No. times predicted\nHRD (prob. \u2265 0.5)',
            y=expression("No. times P(HRD) ">=" 0.5"),
            x='Sample', color='Gene deficiency', 
            title='No. times predicted HRD after 100x 10-fold CV on training set (top 200)'
         ) +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank()
         )
      
      ## Add table grob
      p <- p + draw_grob(
         tableGrob(counts_table, rows=NULL, theme=ttheme_default(base_size=9)), 
         x=n_samples, y=50, hjust=50, vjust=1
      )
      
      return(p)
   }
   
   pdf(paste0(out_dir,'/sel_sample_cv_summary.pdf'), 8, 4)
   grid.draw( plotTimesHrd(sel_sample_cv_summary) )
   dev.off()
}

####################################################################################################
#========= Models =========#
chord_v2_models <- list()

#--------- CHORD v1 clone ---------#
# chord_v2_models[['2.00a']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.00a_CHORDv1_oldTrainSet/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types=c('snv','indel'),rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.00a']])

# chord_v2_models[['2.00b']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.00b_CHORDv1/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types=c('snv','indel'),rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.00b']])

#--------- Expanded indels ---------#
# chord_v2_models[['2.01a']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.01a_expandedIndels/')
#    #dir <- paste0(base_dir,'/CHORDv2/training/models/2.01a_expandedIndels/seed_models/seed04/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types='snv',rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.01a']])

# chord_v2_models[['2.01b']] <- (function(){
#    #dir <- paste0(base_dir,'/CHORDv2/training/models/2.01b_expandedIndels_wilcox1e-10/')
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.01b_expandedIndels_wilcox1e-10/seed_models/seed05/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types='snv',rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.01b']])

#--------- Merged del mh bimh 2-5 ---------#
chord_v2_models[['2.02']] <- (function(){
   dir <- paste0(base_dir,'/CHORDv2/processed/training/models/2.02_mergedDelMh2-5/')
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
main(chord_v2_models[['2.02']])

# models_ct_holdout_dir <- paste0(base_dir,'/CHORDv2/processed/training/models_ct_holdout/')
# for(dir in list.dirs(models_ct_holdout_dir, full.names=T, recursive=F)){
#    message('##',dir)
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- chord_v2_models[['2.02']]$trans.func
#    rf$dir <- dir
#    main(rf)
# }


#--------- Sigs ---------#
chord_v2_models[['2.04b']] <- (function(){
   dir <- paste0(base_dir,'/CHORDv2/processed/training/models/2.04b_mergedDelMh2-5_snvSigs_svSigs/')
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
      l$snv <- fitToSignatures(SBS_SIGNATURE_PROFILES_V2,l$snv, verbose=F)
      l$sv <- fitToSignatures(SV_SIGNATURE_PROFILES,l$sv, verbose=F)
      m <- transformContexts(l, rel.types='all')
      return(m)
   }
   rf$dir <- dir
   return(rf)
})()
main(chord_v2_models[['2.04b']])

#--------- Sv mh ---------#
# chord_v2_models[['2.03a']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.03a_svMh/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){
#       l <- x[c('snv','indel','sv_mh')]
#       names(l) <- c('snv','indel','sv')
#       m <- transformContexts(l, simplify.types=c('snv','indel'), rel.types='all')
#       return(m)
#    }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.03a']], hmf.data.only=T)

# chord_v2_models[['2.03b']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.03b_svMh_expandedIndels/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){
#       l <- x[c('snv','indel','sv_mh')]
#       names(l) <- c('snv','indel','sv')
#       m <- transformContexts(l, simplify.types='snv', rel.types='all')
#       return(m)
#    }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.03b']], hmf.data.only=T)

# chord_v2_models[['2.03c']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.03c_svMh_mergedDelMh2-5/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){
#       #x=contexts[[1]]
#       ## Merge del mh bimh 2-5
#       indels <- x$indel
#       indel_types <- c('del.rep','ins.rep','del.mh','ins.mh','del.none','ins.none')
#       indels_split <- splitDfRegex(indels, indel_types)
#       names(indels_split) <- indel_types
#       
#       counter <- 0
#       indels_custom <- lapply(indels_split, function(i){
#          counter <<- counter + 1
#          df <- as.data.frame(rowSums(i))
#          colnames(df) <- indel_types[counter]
#          return(df)
#       })
#       
#       indels_custom$del.mh <- with(indels_split,{
#          cbind(
#             del.mh['del.mh.bimh.1'],
#             'del.mh.bimh.2.5'=rowSums(del.mh[!(colnames(del.mh) %in% 'del.mh.bimh.1')])
#          )
#       })
#       
#       ## Add new indels back to list
#       l <- x[c('snv','indel','sv_mh')] 
#       names(l) <- c('snv','indel','sv')
#       l$indel <- do.call(cbind, unname(indels_custom))
#       m <- transformContexts(l, simplify.types='snv', rel.types='all')
#       return(m)
#    }
#    rf$dir <- dir
#    return(rf)
# })()
#main(chord_v2_models[['2.03c']], hmf.data.only=T)




