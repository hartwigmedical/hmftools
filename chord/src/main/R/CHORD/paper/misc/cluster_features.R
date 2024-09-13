library(openxlsx)
library(mutSigExtractor)
library(Rtsne)
library(RColorBrewer)
library(ggplot2)

options(stringsAsFactors=F)

#========= Path prefixes =========#
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

#========= Misc functions =========#
forceDfOrder <- function(df){
   as.data.frame(lapply(df, function(i){
      if(is.numeric(i)){ i }
      else { factor(i, unique(i)) }
   }))
}

#========= Load data =========#
pred <- read.xlsx(
   paste0(base_dir,'/CHORDv2/processed/analysis/supp_tables_final/S1_preds.xlsx'),
   sheet='CHORD'
)
pred <- subset(pred, group %in% c('HMF','PCAWG'))

#---------- Features ----------#
contexts <- read.xlsx(
   paste0(base_dir,'/CHORDv2/processed/training/models/2.02_mergedDelMh2-5/plots/features.xlsx'),
   sheet='raw'
)
contexts <- subset(contexts, group %in% c('HMF','PCAWG'))
rownames(contexts) <- contexts$sample
contexts <- contexts[,-c(1:2)] 

contexts <- splitDfRegex(contexts, patterns=list(snv='\\[', indel='[.]', sv='[A-Z]{3}'))

trans.func <- function(features_split){
   ## Merge del mh bimh 2-5
   indels <- features_split$indel
   indel_types <- c('del.rep','ins.rep','del.mh','ins.mh','del.none','ins.none')
   indels_split <- mutSigExtractor::splitDfRegex(indels, indel_types)
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
   l <- features_split[c('snv','indel','sv')]
   l$indel <- do.call(cbind, unname(indels_custom))
   
   ## Simplify SNVs to the 6 types of base substitions (C>A, C>G, C>T, T>A, T>C, T>G) by ignoring
   ## the immediate flanking nucleotides
   ## Convert absolute counts to relative counts. Calculated per variant type.
   m <- mutSigExtractor::transformContexts(l, simplify.types='snv', rel.types='all')
   return(m)
}

contexts_trans <- trans.func(contexts)

#========= Rtsne =========#
pancancer_analysis_samples <- subset(pred,used_for_pancancer_analysis)$sample
m <- contexts_trans[rownames(contexts_trans) %in% pancancer_analysis_samples,]

#m <- m[!grepl('^HMF',rownames(m)),]

set.seed(1)
tnse_out <- Rtsne(m, perplexity=9, check_duplicates=F, verbose=T)

df_tsne_contexts <- (function(){
   pd <- as.data.frame(tnse_out$Y)
   colnames(pd) <- c('x','y')
   
   pd <- cbind(sample=rownames(m), pd)
   
   ## Add metadata
   pd$p_BRCA1 <- pred$p_BRCA1[ match(pd$sample, pred$sample) ]
   pd$p_BRCA2 <- pred$p_BRCA2[ match(pd$sample, pred$sample) ]
   pd$p_hrd <- pred$p_hrd[ match(pd$sample, pred$sample) ]
   
   pd$hr_status <- pred$hr_status[ match(pd$sample, pred$sample) ]
   pd$hr_status <- factor(pd$hr_status, c('cannot_be_determined','HR_proficient','HR_deficient'))
   
   pd$in_training_set <- pred$in_training_set[ match(pd$sample, pred$sample) ]
   pd$in_training_set[is.na(pd$in_training_set)] <- FALSE
   
   pd$cancer_type <- pred$cancer_type[ match(pd$sample, pred$sample) ]
   
   ## Assign cancer type colors
   tab <- table(pred$cancer_type)
   low_freq_cancer_types <- names(tab)[tab<50]
   
   pd$cancer_type_simple <- pd$cancer_type
   pd$cancer_type_simple[pd$cancer_type_simple %in% low_freq_cancer_types] <- 'Other'
   
   cancer_types <- table(pd$cancer_type_simple)
   cancer_types <- sort(cancer_types, decreasing=T)
   cancer_types <- names(cancer_types)
   color_pal <- c(
      brewer.pal(9,'Set1'),
      brewer.pal(12,'Set3'),
      brewer.pal(9,'Pastel1')
   )
   cancer_type_colors <<- structure(
      color_pal[1:length(cancer_types)],
      names=cancer_types
   )
   
   pd <- do.call(rbind,split(pd, pd$hr_status))
   rownames(pd) <- NULL
   pd$sample <- factor(pd$sample, pd$sample)
   
   ##
   pd <- cbind(
      pd, 
      m[match(pd$sample,rownames(m)),1:13]
   )
   
   return(pd)
})()

#--------- tSNE plot ---------#
# HLINE <- -20
# VLINE <- -13


p_tsne_contexts <- ggplot(df_tsne_contexts, aes(x, y)) +
   #geom_point(aes(color=hr_status, fill=cancer_type_simple, size=in_training_set), shape=21, stroke=0.7) +
   geom_point(aes(color=hr_status, fill=cancer_type_simple), shape=21, stroke=0.7) +
   
   # geom_vline(xintercept=HLINE) +
   # geom_hline(yintercept=VLINE) +
   
   #geom_segment(aes(x=-Inf, xend=HLINE, y=VLINE, yend=VLINE), color='black', linetype='dotted') + ## hline
   #geom_segment(aes(x=HLINE, xend=HLINE, y=-Inf, yend=VLINE), color='black', linetype='dotted') + ## vline
   
   scale_fill_manual(name='[Fill] Cancer type', values=cancer_type_colors) +
   scale_color_manual(name='[Outline] CHORD HR status', values=c('HR_deficient'='black','HR_proficient'='lightgrey','cannot_be_determined'='white')) +
   #scale_size_manual(name='[Size] In training set', values=c('TRUE'=2,'FALSE'=1)) +
   
   scale_y_continuous(breaks=seq(-100,100,by=10)) +
   scale_x_continuous(breaks=seq(-100,100,by=10)) +
   
   ylab('tSNE dimension 2') + xlab('tSNE dimension 1') +
   
   theme_bw() +
   theme(
      legend.position='left'
   )

pdf(paste0(base_dir,'/CHORDv2/processed/analysis/rebuttal_natgen/plots/tsne_features.pdf'),10,7)
plot(p_tsne_contexts)
dev.off()

#--------- Get HRD group ---------#
sample_subsets <- list()

sample_subsets$tsne_hrd <- subset(
   df_tsne_contexts,
   (y >= 23 & x >= 28)
   | (y >= 30 & x >= 25)
   | (y >= 35 & x >= 20)
   | (y >= 45 & x >= 4),
   sample, drop=T
)

sample_subsets$chord_hrd <- subset(df_tsne_contexts, hr_status=='HR_deficient', sample, drop=T)

sample_subsets$chord_hrd[ !(sample_subsets$tsne_hrd %in% sample_subsets$chord_hrd) ]

p_tsne_hrd_cancer_types <- (function(){
   pd <- do.call(rbind, lapply(names(sample_subsets), function(i){
      #i=sample_subsets$tsne_hrd
      
      cancer_type_totals <- table(df_tsne_contexts$cancer_type)
      
      hrd_counts <- structure(rep(0, length(cancer_type_totals)),names=names(cancer_type_totals))
      tab <- table(subset(df_tsne_contexts, sample %in% sample_subsets[[i]], cancer_type, drop=T))
      hrd_counts[names(tab)] <- tab
      
      df <- data.frame(
         cancer_type=names(cancer_type_totals),
         abs=as.numeric(hrd_counts),
         rel=as.numeric(hrd_counts)/as.numeric(cancer_type_totals),
         total=as.numeric(cancer_type_totals),
         group=i
      )
      
      df <- df[order(df$rel, decreasing=T),]
      return(df)
   }))
   
   pd <- forceDfOrder(pd)
   pd$label <- with(pd,{
      paste0(
         abs,' / ',total,' ' #,
         #'(',signif(100*rel,2),'%)'
      )
   })
   
   ggplot(pd, aes(x=cancer_type, y=rel)) +
      facet_grid(group~., switch='y') +
      geom_bar(stat='identity', color='black', fill='#639E60') +
      geom_text(aes(y=0.01, label=label), angle=90, hjust=0, vjust=0.5) +
      scale_y_continuous(labels=function(x){ paste0(x*100,'%') }, name='Frequency') +
      theme_bw() +
      theme(
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         axis.title.x=element_blank(),
         strip.placement='outside'
      )
   
})()

pdf(paste0(base_dir,'/CHORDv2/processed/analysis/rebuttal_natgen/plots/tsne_hrd_cancer_types.pdf'),9,4)
plot(p_tsne_hrd_cancer_types)
dev.off()
