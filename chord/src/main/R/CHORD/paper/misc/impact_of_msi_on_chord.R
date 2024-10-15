library(mutSigExtractor)
library(ggplot2)
library(scales)
library(grid)
#library(ggpubr)

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
metadata <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_pred_ann.txt.gz'))

chord_dir <- paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/')
pred <- readRDS(paste0(chord_dir,'/plots/preds.rds'))
pred <- pred$HMF

#pred <- pred[pred$sample %in% metadata[metadata$used_for_pancancer_analysis,'sample'],]

selected_samples <- pred$sample

#--------- Contexts ---------#
contexts <- list()

contexts$abs <- readRDS(paste0(base_dir,'/datasets/HMF_DR010_DR047/matrices/contexts_merged.rds'))
contexts$abs <- transformContexts(contexts$abs, simplify.types=c('snv','indel'), export.list=T)

contexts$rel <- transformContexts(contexts$abs, rel.types='all', export.list=T)

#--------- Annotation ---------#
#pred_hmf <- readRDS(paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/plots/preds.rds'))$HMF

diplotypes <- (function(){
   path <- '/Users/lnguyen/Documents/R_cache/hmf_gene_diplotypes_max.txt.gz'
   
   if(!file.exists(path)){
      file.copy(
         paste0(base_dir,'/datasets/HMF_DR010_DR047/scripts/annotate_genes/gene_diplotypes_with_amps_HMF_DR010_DR047.txt.gz'),
         path
      )
   }
   
   read.delim(path)
})()

#sel_genes <- c('BRCA2','BRCA1','RAD51C','PALB2')
sel_genes <- c('BRCA2','BRCA1')
diplotypes_ss <- diplotypes[
   diplotypes$hgnc_symbol %in% sel_genes 
   #& diplotypes$sample %in% selected_samples
,]

diplotypes_split <- split(diplotypes_ss, diplotypes_ss$sample)

gene_defs <- unlist(lapply(diplotypes_split, function(i){
   gene_def_init <- 'none'
   gene_def <- i$hgnc_symbol[i$hit_score==10 & !(i$diplotype_origin %in% c('germ_som','som_som'))]
   if(length(gene_def)!=0){
      gene_def_init <- gene_def[1]
   }
   return(gene_def_init)
}))

annotation <- data.frame(
   sample=names(gene_defs),
   response=gene_defs,
   row.names=NULL
)

#annotation$custom_blacklist_checked <- annotation$sample %in% unlist(custom_blacklist_checked)


#========= Plot  =========#
#--------- Prep data ---------#
df <- pred[c('sample','hrd')]
df$p_rank <- order(df$hrd)
df <- merge(df,annotation,sort=F)

df$response2 <- df$response
#df$response2[df$custom_blacklist_checked] <- 'none'

df$del.mh <- contexts$abs$indel[match(df$sample, rownames(contexts$abs$indel)),'del.mh'] 
df$del.mh_rel <- contexts$rel$indel[match(df$sample, rownames(contexts$rel$indel)),'del.mh'] 

df <- cbind(
   df,
   (function(){
      m <- as.data.frame(contexts$abs$indel[,c('del.rep','ins.rep')] )
      m$indel.rep <- apply(m,1,sum)
      m$indel.all <- rowSums(contexts$abs$indel)
      m$indel.rep_rel <- m$indel.rep / m$indel.all
      
      m[match(df$sample, rownames(m)),]
   })()
)

df <- df[order(df$response, decreasing=T),] ## Force plotting of BRCA1/2 deficient samples last
df <- as.data.frame(lapply(df, function(i){
   if(!is.numeric(i)){ i <- factor(i, unique(i)) }
   return(i)
}))

## Plot points with gene deficiencies last
df <- do.call(rbind,split(df,df$response2!='none'))

df$hr_status <- ifelse(df$hrd >= 0.5, 'hrd', 'hrp')
df$has_msi <- df$indel.rep>=14000
df$is_fn_hrp <- with(df,{ hr_status=='hrp' & response!='none' })

df_ss <- subset(df, is_fn_hrp)

## Ordering
df <- forceDfOrder(df)

gene_colors <- c(BRCA1='#f58225', BRCA2='#69439d', RAD51C='green', PALB2='magenta', none='grey')
df$response2 <- factor(df$response2, names(gene_colors))

#subset(df, indel.rep>14000 & hrd<0.5 & response!='none')$sample


#--------- indel.rep vs. HRD score ---------#
plot_p_chord_vs_msi <- ggplot(df, aes(indel.rep, hrd, color=response2)) + 
   geom_point() +
   geom_vline(xintercept=14000, linetype='dotted') +
   annotate('text',x=14000*1.1, y=0.5, label='MSI positive:\nindel.rep > 14000', hjust=0) +
   
   scale_color_manual(
      values=gene_colors,
      breaks=names(gene_colors)
   ) +
   scale_x_log10(labels=comma) +
   labs(y='Probability of HRD', x='indel.rep (abs. contrib.)', color='Gene deficiency') +
   #guides(guide_legend(reverse=F)) +
   theme_bw()
   


pdf(paste0(base_dir,'/CHORDv2/analysis/misc/plots/p_chord_vs_indel_rep.pdf'), 8, 5)
grid.draw(plot_p_chord_vs_msi)
dev.off()

# #--------- del.mh ---------#
# df_melt <- (function(){
#    id_vars <- c('del.mh','del.mh_rel')
#    #df_ss <- df[df$in_training_set==T,]
#    do.call(rbind,lapply(id_vars, function(i){
#       out <- df[,c(i,'response2','has_msi')]
#       colnames(out)[1] <- 'value'
#       out$variable <- i
#       return(out)
#    }))
# })()
# 
# ggplot(df_melt,aes(response2, value)) +
#    facet_wrap(~variable, ncol=1, scales='free_y', strip.position='left') +
#    geom_boxplot(outlier.color=NA) +
#    geom_jitter(aes(color=has_msi), width=0.3) +
#    scale_color_manual(
#       values=c('TRUE'='red', 'FALSE'='darkgrey'),
#       breaks=c('TRUE','FALSE')
#    ) +
#    #scale_y_continuous(expand=c(0,0,0.5,0)) +
#    stat_compare_means(ref.group='none', label='p.format') +
#    labs(y='del.mh contribution', x='Gene deficiency', color='Has MSI') +
#    theme_bw() +
#    theme(
#       
#    )
# 
# 
# 
# ggplot(df,aes(response2, del.mh)) +
#    geom_boxplot(outlier.color=NA) +
#    geom_jitter(aes(color=has_msi), width=0.3) +
#    scale_color_manual(
#       values=c('TRUE'='red', 'FALSE'='darkgrey'),
#       breaks=c('TRUE','FALSE')
#    ) +
#    stat_compare_means(ref.group='none') +
#    labs(y='del.mh (rel. contrib.)', x='Gene deficiency', color='Has MSI') +
#    theme_bw()
# 
# 
# 
# ggplot(df,aes(del.mh, p_chord, color=response2)) +
#    geom_point() +
#    scale_color_manual(
#       values=gene_colors,
#       breaks=names(gene_colors)
#    ) +
#    #scale_x_continuous(labels=comma) +
#    labs(y='Probability of HRD', x='del.mh (abs. contrib.)', color='Gene deficiency') +
#    #guides(guide_legend(reverse=F)) +
#    theme_bw()








#========= PCA =========#
ann_msi_samples <- annotation[ annotation$has_msi, ]
subset(ann_msi_samples, response!='none')

m <- do.call(cbind, unname(contexts$abs))
m <- m[ann_msi_samples$sample,c('del.rep','ins.rep','del.mh','ins.mh','del.none','ins.none')]

pca <- prcomp(m, center=T, scale=T)

df <- as.data.frame(pca$x)
df <- cbind(sample=rownames(df), df); rownames(df) <- NULL

df$response <- ann_msi_samples$response[ match(df$sample, ann_msi_samples$sample) ]

ggplot(df, aes(PC1,PC2, color=response)) +
   geom_point()

# plot_p_chord_vs_del.mh_vs_msi <- ggplot(df, aes(del.mh_rel, p_chord, fill=response2, color=has_msi)) + 
#    geom_point(shape=21, size=2.5) +
#    #geom_vline(xintercept=msi_cutoff, linetype='dotted') +
#    #annotate('text',x=msi_cutoff*1.1, y=0.5, label='MSI positive:\nindel.rep > 14000', hjust=0) +
#    
#    scale_fill_manual(values=c(BRCA1='#f58225', BRCA2='#69439d', none='grey')) +
#    scale_color_manual(values=c('TRUE'='red', 'FALSE'='black')) +
#    #scale_x_log10(labels=comma) +
#    labs(y='Probability of HRD', x='del.mh (rel. contrib.)', fill='Gene deficiency', color='Has MSI') +
#    theme_bw()
# 
# #pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/p_chord_vs_del.mh_vs_msi.pdf'), 8, 5)
# pdf('/Users/lnguyen/Desktop/p_chord_vs_del.mh_vs_msi.pdf', 8, 5)
# grid.draw(plot_p_chord_vs_del.mh_vs_msi)
# dev.off()








