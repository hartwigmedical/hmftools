library(CHORD)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(randomForest)
library(cowplot)
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
      if(is.numeric(i)){ i }
      else { factor(i, unique(i)) }
   }))
}

selectCommonSamplesInList <- function(l, by='rownames'){
   if(by=='rownames'){
      common_samples <- Reduce(intersect, lapply(l,rownames))
      lapply(l, function(i){ i[common_samples,] })
   } else {
      common_samples <- Reduce(intersect, lapply(l,`[[`, by))
      lapply(l, function(i){ i[match(common_samples, i[[by]]),] })
   }
}

#========= Load features =========#
## SV mutations don't have clonal/subclonal split
contexts_full <- readRDS(paste0(base_dir,'/datasets/HMF_DR010_DR047/matrices/contexts_merged.rds'))

## Read clonality contexts from Arne
contexts_clonality_raw <- (function(){
   dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/analysis/matrices/matrices_clonal_subclonal/matrices/'
   #dir2 <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_DR010_DR047/analysis/rebuttal_natgen/subclonal_analysis/data/'
   
   paths <- list(
      all=list(
         snv=paste0(dir,'HMF_mutSigExtractor_SNV_ALL.rds'),
         indel=paste0(dir,'HMF_mutSigExtractor_INDEL_ALL.rds')#,
         #mnv=paste0(dir2,'HMF_mutSigExtractor_MNV2_ALL.rds')
      ),
      
      subclonal=list(
         snv=paste0(dir,'HMF_mutSigExtractor_SNV_SUBCLONAL.rds'),
         indel=paste0(dir,'HMF_mutSigExtractor_INDEL_SUBCLONAL.rds')#,
         #mnv=paste0(dir2,'HMF_mutSigExtractor_MNV2_SUBCLONAL.rds')
      ),
      
      clonal=list(
         snv=paste0(dir,'HMF_mutSigExtractor_SNV_CLONAL.rds'),
         indel=paste0(dir,'HMF_mutSigExtractor_INDEL_CLONAL.rds')#,
         #mnv=paste0(dir2,'HMF_mutSigExtractor_MNV2_CLONAL.rds')
      )
   )
   
   lapply(paths, function(i){
      #i=paths$all
      i <- lapply(i, function(j){
         #j=i[[3]]
         m <- readRDS(j)
         m <- t(m)
         m <- m[rownames(m)!='empty',]
         as.data.frame(m)
      })
      
      ## Append SV contexts
      i$sv <- as.matrix(contexts_full$sv) 
      selectCommonSamplesInList(i)
   })
})()


#========= Sample groups =========#
metadata <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_pred_ann.txt.gz'))

#UNIQ_SAMPLES <- metadata[metadata$max_purity_biopsy,'sample']

WHITELIST_SAMPLES <- subset(metadata, in_training_set, sample, drop=T)

GROUPS <- list(
   had_radio=metadata[metadata$has_radiotherapy_pre_treatment=='Yes','sample_id'],
   brca_def=metadata[metadata$response %in% c('BRCA1','BRCA2'),'sample']
)

#========= Feature importance =========#
contexts_clonality_trans <- lapply(contexts_clonality_raw, function(i){
   #i=contexts_clonality_raw[['subclonal']]
   l <- i[c('snv','indel','sv')]
   
   snv_load <- rowSums(l$snv)
   indel_load <- rowSums(l$indel)
   sel_samples <- intersect(
      names(snv_load)[snv_load>=100],
      names(indel_load)[indel_load>=50]
   )
   sel_samples <- intersect(WHITELIST_SAMPLES, sel_samples)
   
   l <- transformContexts(l, simplify.types='snv',rel.types='all', export.list=T)
   
   l <- lapply(l, function(j){
      #j=l[[2]]
      m <- j[rownames(j) %in% sel_samples,]
      #m <- m/rowSums(m)
      #m[is.na(m)] <- 0
      #rownames(m) <- rownames(j)
      return(m)
   })
   
   as.matrix(do.call(cbind, unname(l)))
   
})

#--------- Train random forests ---------#
detFeatImp <- function(x, y, seed=1, balance.ratio=NULL, rf.name=NULL){
   
   # x=contexts_clonality_trans$all
   # y=rownames(x) %in% GROUPS$had_radio
   # top.n=30
   # seed=1
   set.seed(seed)
   if(!is.factor(y)){ y <- factor(y, c('TRUE','FALSE')) }
   
   ## Remove features that:
   ## - decrease in pos class compared to neg class
   ## - have zero change
   x_split <- split(as.data.frame(x), y)
   median_diffs <- apply(x_split$`TRUE`,2,median) - apply(x_split$`FALSE`,2,median)
   sel_features <- names(median_diffs)[ sign(median_diffs)==1 ]
   
   df <- cbind(as.data.frame(x), response=y)
   
   ## Balance classes
   
   if(!is.null(balance.ratio)){
      df_split <- split(df, df$response) #balance.ratio
      names(df_split) <- c('pos','neg')
      
      class_counts <- lapply(df_split, nrow)
      
      if(class_counts$pos > class_counts$neg){
         target_count <- class_counts$neg*balance.ratio
         if(target_count >= class_counts$neg){ target_count <- class_counts$neg }
         
         df_split$pos <- df_split$pos[
            sample(1:class_counts$pos, target_count)
            ,]
      } else if(class_counts$neg > class_counts$pos){
         target_count <- class_counts$pos*balance.ratio
         if(target_count >= class_counts$pos){ target_count <- class_counts$pos }
         
         df_split$neg <- df_split$neg[
            sample(1:class_counts$neg, target_count)
            ,]
      }
      df <- rbind(df_split$pos, df_split$neg)
   }
   
   ## Run RF and get importance
   rf <- randomForest(df[,sel_features], df$response, importance=T, ntree=500)
   
   if(!is.null(rf.name)){
      rf$name <- rf.name
   }
   
   return(rf)
}

l_rf <- lapply(names(contexts_clonality_trans), function(i){
   
   message('\n## ',i,' variants')
   #which_variants <- i
   
   rf_list <- list()
   x <- contexts_clonality_trans[[i]]
   #y <- rownames(x) %in% GROUPS$brca_def
   #y <- rownames(x) %in% GROUPS$had_radio
   
   message('Positive control') #---------
   rf_list[[1]] <- detFeatImp(
      x, rownames(x) %in% GROUPS$brca_def,
      rf.name=c(paste0(i,' variants'),'all samples; BRCA def?')
   )
   
   message('Radio vs no radio') #---------
   rf_list[[2]] <- detFeatImp(
      x, rownames(x) %in% GROUPS$had_radio,
      rf.name=c(paste0(i,' variants'),'all samples; had radio?')
   )
   
   ## BRCA prof
   rf_list[[3]] <- (function(){
      x2 <- x[!(rownames(x) %in% GROUPS$brca_def),]
      detFeatImp(
         x2, rownames(x2) %in% GROUPS$had_radio,
         rf.name=c(paste0(i,' variants'),'BRCA prof samples; had radio?')
      )
   })()
   
   return(rf_list)
})
names(l_rf) <- names(contexts_clonality_trans)

#saveRDS(l_rf, paste0(base_dir,'/CHORDv2/analysis/radiotherapy_features/data/l_rf.rds'))

#--------- Plot importance ---------#
plotFeatImp <- function(rf, top.n=15, plot.title=NULL){
   #rf=l_rf[[1]][[1]]
   
   imp <- as.data.frame(rf$importance)
   imp <- imp[order(imp$MeanDecreaseGini, decreasing=T),]
   
   ## Return
   imp <- data.frame(
      feature=rownames(imp),
      imp=imp$MeanDecreaseGini
   )
   
   if(top.n>nrow(imp)){ top.n <- nrow(imp) }
   imp_ss <- forceDfOrder(imp[1:top.n,])
   imp_ss$label <- paste0(' ',imp_ss$feature,' [',round(imp_ss$imp,1),']')
   
   oob_err_rate <- paste0(
      round(rf$err.rate[nrow(rf$err.rate),'OOB']*100, 2),'% / ',
      round(rf$err.rate[nrow(rf$err.rate),'TRUE']*100, 2),'% / ',
      round(rf$err.rate[nrow(rf$err.rate),'FALSE']*100, 2),'%'
   )
   
   p <- ggplot(imp_ss, aes(x=feature,y=imp)) +
      geom_bar(stat='identity', fill='#EE9E9B', color='black') +
      geom_text(aes(y=0, label=label), hjust=0, vjust=0.5, angle=90, size=3) +
      #scale_x_discrete(position='top') +
      #xlab(paste0('Top ',top.n, ' / ', nrow(imp),' features || OOB err all/pos/neg: ',oob_err_rate)) +
      ylab('Feature importance\n(Mean decrease in Gini)') +
      theme_bw() +
      theme(
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank()
      )
   
   if(is.null(plot.title)){
      p <- p + ggtitle(rf$name[1], rf$name[2])
   }
   
   return(p)
}

counter <- 0
for(i in l_rf){
   counter <- counter + 1
   file_suffix <- names(l_rf)[counter]

   pdf(paste0(base_dir,'/CHORDv2/analysis/radiotherapy_features/plots/rf_imp_',file_suffix,'_variants.pdf'), 5, 2.5)
   for(j in i){ plot(plotFeatImp(j)) }
   dev.off()

   message(file_suffix,' done')
}

# #========= Manual comparison =========#
# features <- colnames(contexts_clonality_trans$all)
# 
# 
# df_melt <- melt(contexts_clonality_trans)
# colnames(df_melt) <- c('sample','feature','contrib','which_variants')
# 
# sel_features <- paste0('del.mh.bimh.',1:5)
# df_melt <- subset(df_melt, feature %in% sel_features)
# 
# df_melt$brca_def <- ifelse(df_melt$sample %in% GROUPS$brca_def, 'brca_def','brca_prof')
# df_melt$had_radio <- ifelse(df_melt$sample %in% GROUPS$had_radio, 'had_radio','no_radio')
# 
# df_melt$group <- paste0(df_melt$brca_def,', ',df_melt$had_radio)
# 
# pd <- subset(df_melt, which_variants=='clonal')
# pd <- forceDfOrder(pd)
# 
# p <- ggplot(pd, aes(group, contrib)) +
#    facet_wrap(~feature, scales='free_y',ncol=5) +
#    geom_jitter(width=0.25, shape=21, color='darkgrey') +
#    geom_boxplot(aes(fill=group),outlier.shape=NA,alpha=0.9) +
#    
#    stat_compare_means(
#       comparisons=combn(levels(pd$group),2,simplify=F),
#       method.args=list(alternative='greater'),
#       label='p.format', size=3
#    ) +
#    #ggtitle(paste0(i,' variants')) +
#    
#    theme_bw() +
#    theme(
#       axis.title.x=element_blank(),
#       axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#       legend.position='top',
#       legend.justification=c(0,0.5)
#    )
# 
# 
# 
# l_p_feature_comparisons <- lapply(names(contexts_clonality_trans), function(i){
#    ## Prepare data for plotting
#    #m <- contexts_clonality_trans$subclonal[WHITELIST_SAMPLES,sel_features]
#    m <- contexts_clonality_trans[[i]][rownames(contexts_clonality_trans[[i]]) %in% WHITELIST_SAMPLES,sel_features]
#    df <- melt(m)
#    colnames(df) <- c('sample','feature','contrib')
# 
#    df <- cbind(
#       df,
#       as.data.frame(lapply(GROUPS, function(i){ df$sample %in% i }))
#    )
# 
#    ## Assign confusion groups
#    df$confuse_groups <- (function(){
#       sel_groups <- c('clonal_hrd','subclonal_fp','subclonal_tn')
#       v <- sel_groups[ apply(df[,sel_groups],1,which.max) ]
#       factor(v, rev(sel_groups))
#    })()
# 
#    p_had_radio <- ggplot(df, aes(had_radio, contrib)) +
#       facet_wrap(~feature, scales='free_y',ncol=5) +
#       geom_jitter(width=0.25, shape=21, color='darkgrey') +
#       geom_boxplot(aes(fill=had_radio),outlier.shape=NA,alpha=0.9) +
# 
#       stat_compare_means(
#          method.args=list(alternative='greater'),
#          label='p.format', size=3
#       ) +
#       ggtitle(paste0(i,' variants')) +
# 
#       theme_bw() +
#       theme(
#          axis.title.x=element_blank(),
#          legend.position='top',
#          legend.justification=c(0,0.5)
#       )
# 
#    p_confuse_groups <- ggplot(df, aes(confuse_groups, contrib)) +
#       facet_wrap(~feature, scales='free_y',ncol=5) +
#       geom_jitter(width=0.25, shape=21, color='darkgrey') +
#       geom_boxplot(aes(fill=confuse_groups),outlier.shape=NA,alpha=0.9) +
# 
#       scale_y_continuous(expand=c(0.01, 0.01, 0.01, 0.08)) +
#       stat_compare_means(
#          comparisons=combn(levels(df$confuse_groups),2,simplify=F),
#          method.args=list(alternative='less'), ## For some reason, with all vs all, 'less'/'greater' is inverted
#          label='p.format', size=3
#       ) +
#       ggtitle(paste0(i,' variants')) +
# 
#       theme_bw() +
#       theme(
#          axis.title.x=element_blank(),
#          axis.text.x=element_text(angle=15, hjust=1),
#          legend.position='top',
#          legend.justification=c(0,0.5)
#       )
# 
#    plot_grid(p_had_radio, p_confuse_groups, ncol=1)
# })
# 
# # pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/rebuttal_natgen/subclonal_analysis/plots/wilcox_v2.pdf'), 10, 12)
# # for(i in l_p_feature_comparisons){ grid.draw(i) }
# # dev.off()

