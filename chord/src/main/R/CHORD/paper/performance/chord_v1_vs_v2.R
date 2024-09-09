library(CHORD)
library(reshape2)
library(ggplot2)
#library(ggpubr)
library(randomForest)
library(cowplot)
library(grid)
library(gridExtra)
#library(Boruta)
#library(glmnet)

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
## HMF
contexts_clonality_hmf <- (function(){
   dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/analysis/matrices/matrices_clonal_subclonal/matrices/'

   ## SV mutations don't have clonal/subclonal split
   contexts_full <- readRDS(paste0(base_dir,'/CHORD/HMF_DR010_DR047/matrices/merged_contexts.rds'))
   contexts_sv_mh <- read.delim(paste0(base_dir,'/CHORD/HMF_DR010_DR047/scripts/get_sv_type_and_homology/sv_mh_contexts.txt.gz'))
   
   paths <- list(
      
      clonal=list(
         snv=paste0(dir,'HMF_mutSigExtractor_SNV_CLONAL.rds'),
         indel=paste0(dir,'HMF_mutSigExtractor_INDEL_CLONAL.rds')
      ),
      
      subclonal=list(
         snv=paste0(dir,'HMF_mutSigExtractor_SNV_SUBCLONAL.rds'),
         indel=paste0(dir,'HMF_mutSigExtractor_INDEL_SUBCLONAL.rds')
      ),
      
      all=list(
         snv=paste0(dir,'HMF_mutSigExtractor_SNV_ALL.rds'),
         indel=paste0(dir,'HMF_mutSigExtractor_INDEL_ALL.rds')
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
      i$sv_mh <- as.matrix(contexts_sv_mh)
      
      selectCommonSamplesInList(i)
   })
})()

contexts_clonality_pcawg <- (function(){
   dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/PCAWG_2020/matrices/'

   ## SV mutations don't have clonal/subclonal split
   contexts_full <- readRDS(paste0(base_dir,'/datasets/PCAWG_2020/matrices/contexts_merged.rds'))
   contexts_sv_mh <- read.delim(paste0(base_dir,'/datasets/PCAWG_2020/scripts/get_sv_type_and_homology/sv_mh_contexts.txt.gz'))
   
   paths <- list(
      
      clonal=list(
         snv=paste0(dir,'PCAWG_mutSigExtractor_SNV_CLONAL.rds'),
         indel=paste0(dir,'PCAWG_mutSigExtractor_INDEL_CLONAL.rds')
      ),
      
      subclonal=list(
         snv=paste0(dir,'PCAWG_mutSigExtractor_SNV_SUBCLONAL.rds'),
         indel=paste0(dir,'PCAWG_mutSigExtractor_INDEL_SUBCLONAL.rds')
      ),
      
      all=list(
         snv=paste0(dir,'PCAWG_mutSigExtractor_SNV_ALL.rds'),
         indel=paste0(dir,'PCAWG_mutSigExtractor_INDEL_ALL.rds')
      )
   )
   
   lapply(paths, function(i){
      #i=paths$all
      i <- lapply(i, function(j){
         #j=i[[3]]
         m <- readRDS(j)
         m <- t(m)
         m <- m[rownames(m)!='empty',]
         
         ##Fix sample names
         rownames(m) <- gsub('^X','',rownames(m))
         rownames(m) <- gsub('[.]','-',rownames(m))
         
         as.data.frame(m)
      })
      
      ## Append SV contexts
      i$sv <- as.matrix(contexts_full$sv) 
      i$sv_mh <- as.matrix(contexts_sv_mh)
      
      selectCommonSamplesInList(i)
   })
})()

#========= Load other data =========#
metadata_hmf <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))
metadata_hmf$biallelic_status[grep('^[0-3];',metadata_hmf$biallelic_status)] <- '0;none'

metadata_pcawg <- read.delim(paste0(base_dir,'/datasets/PCAWG_2020/metadata/pcawg_metadata_ann.txt.gz'))

metadata_pcawg$biallelic_status <- factor(
   metadata_pcawg$biallelic_status,
   c('none','GMutation/SMutation','GMutation/SDeletion','SDeletion/SMutation','SDeletion/SSV','SDeletion/SDeletion')
)

getMetadata <- function(sample_names,col){
   metadata[match(sample_names,metadata$sample_id),col]
}

#========= Main =========#
plotV1vsV2Comparison <- function(rf, contexts_clonality, metadata, out_dir=NULL){
   # rf = (function(){
   #    dir <- paste0(base_dir,'/CHORDv2/training/models/2.00a_CHORDv1_oldTrainSet/')
   #    rf <- readRDS(paste0(dir,'/rf_model.rds'))
   #    rf$trans.func <- function(x){ transformContexts(x, simplify.types=c('snv','indel'),rel.types='all') }
   #    rf$dir <- dir
   #    return(rf)
   # })()
   # rf=chord_v2_models[['2.02']]
   # contexts_clonality=contexts_clonality_hmf
   # metadata=metadata_hmf
   
   # contexts_clonality=contexts_clonality_pcawg
   # metadata=metadata_pcawg
   # out_dir=paste0(base_dir,'/CHORDv2/training/models/2.00b_CHORDv1/plots/')
   
   #========= CHORD predictions =========#
   preds <- list()
   
   preds$v1 <- lapply(contexts_clonality, function(i){
      #i=contexts_clonality$all
      chord_v1 <- readRDS(paste0(base_dir,'/CHORDv2/training/models/2.00b_CHORDv1//rf_model.rds'))
      m <- do.call(cbind, unname(i[c('snv','indel','sv')]))
      pred <- chordPredict(
         rf.model=chord_v1,
         m, trans.func=function(x){ transformContexts(x, simplify.types=c('snv','indel'),rel.types='all') },
         verbose=F
      )
      return(pred)
   })
   
   preds$v2 <- lapply(contexts_clonality, function(i){
      #i=contexts_clonality$all
      #m=rf$trans.func(i)
      pred <- chordPredict(
         i, 
         rf.model=rf,
         trans.func=rf$trans.func,
         verbose=F
      )
      return(pred)
   })
   
   #========= Plot clonality prediction comparisons =========#
   FILL_COLORS <- c(
      ## HMF
      '0;none'='white',
      '4;2x_pathogenic'='green',
      '5;loh+vus'='cyan',
      '6;loh+likely_pathogenic'='orange',
      '7;loh+pathogenic'='red',
      '8;cn_loss'='magenta',
      
      ## PCAWG
      'none'='white',
      'GMutation/SMutation'='green',
      'GMutation/SDeletion'='orange',
      'SDeletion/SMutation'='red',
      'SDeletion/SSV'='pink',
      'SDeletion/SDeletion'='magenta'
   )
   
   mkPredPairs <- function(l, rm.failed.qc=T){
      #l=preds$v1
      
      pairs <- combn(names(l),2,simplify=F)
      
      l <- lapply(pairs, function(i){
         #i=pairs[[1]]
         
         l_pred <- list(
            x = l[[ i[1] ]],
            y = l[[ i[2] ]]
         )
         
         l_pred <- selectCommonSamplesInList(l_pred, 'sample')
         
         if(rm.failed.qc){
            sel_rows <- l_pred$x$hr_status != 'cannot_be_determined' & l_pred$y$hr_status != 'cannot_be_determined'
            l_pred <- lapply(l_pred, function(j){ j[sel_rows,] })
         }
         
         out <- data.frame(
            sample = l_pred$x$sample,
            p_hrd.x = l_pred$x$p_hrd,
            p_hrd.y = l_pred$y$p_hrd
         )
         out$brca_def <- metadata[match(out$sample, metadata$sample),'response']
         out$brca_biall_status <- metadata[match(out$sample, metadata$sample),'biallelic_status']
         out$in_training_set <- metadata[match(out$sample, metadata$sample),'in_training_set']
         
         if('has_radiotherapy_pre_treatment' %in% colnames(metadata)){
            out$has_radiotherapy_pre_treatment <- metadata[match(out$sample, metadata$sample),'has_radiotherapy_pre_treatment']
            out$has_radiotherapy_pre_treatment[nchar(out$has_radiotherapy_pre_treatment)==0] <- NA
            
            out$has_systemic_pre_treatment <- metadata[match(out$sample, metadata$sample),'has_systemic_pre_treatment']
            out$has_systemic_pre_treatment[nchar(out$has_systemic_pre_treatment)==0] <- NA
            
            out$has_any_pre_treatment <- ifelse(
               out$has_radiotherapy_pre_treatment=='Yes' | out$has_systemic_pre_treatment=='Yes',
               'Yes','No'
            )
         }
         
         return(out)
      })
      names(l) <- unlist(lapply(pairs, paste, collapse='_'))
      
      if(length(l)==1){
         return(l[[1]])
      } else {
         return(l)
      }
   }
   
   calcQuadrantCounts <- function(df, x.cutoff=0.5, y.cutoff=0.5){
      with(df,{
         l <- list(
            nw = p_hrd.x < x.cutoff & p_hrd.y >= y.cutoff,
            ne = p_hrd.x >= x.cutoff & p_hrd.y >= y.cutoff,
            
            sw = p_hrd.x < x.cutoff & p_hrd.y < y.cutoff,
            se = p_hrd.x >= x.cutoff & p_hrd.y < y.cutoff
         )
         v <- unlist(lapply(l, sum))
         matrix(v,nrow=2,ncol=2,byrow=T)
      })
   }
   
   plotPredPair <- function(df, axis.labs=NULL, subtitle=NULL, x.cutoff=0.5, y.cutoff=0.5, swap.x.y=F){
      # pair.name='clonal_subclonal'
      # df=pred_pairs$v1$clonal_subclonal
      # x.cutoff=0.5
      # y.cutoff=0.5
      
      # df <- do.call(rbind,split(df, grepl('none', df$brca_biall_status)))
      # df$sample <- factor(df, unique(df))
      
      if(swap.x.y){
         x_col <- colnames(df)=='p_hrd.x'
         y_col <- colnames(df)=='p_hrd.y'
         
         colnames(df)[x_col] <- 'p_hrd.y'
         colnames(df)[y_col] <- 'p_hrd.x'
      }
      
      quadrant_counts <- calcQuadrantCounts(df, x.cutoff, y.cutoff)
      
      p <- ggplot(df, aes(p_hrd.x, p_hrd.y)) +
         geom_hline(yintercept=y.cutoff, linetype='dotted') +
         geom_vline(xintercept=x.cutoff, linetype='dotted') +
         
         # geom_point(aes(fill=brca_def, color=in_training_set),shape=21) +
         # scale_fill_manual(values=c(BRCA1='red',BRCA2='blue',none='white')) +
      
         geom_point(
            data=df[grepl('none',df$brca_biall_status),],
            aes(fill=brca_biall_status, color=if(!is.null(df$in_training_set)){ in_training_set } else { TRUE }),
            shape=21
         ) +
         
         geom_point(
            data=df[!grepl('none',df$brca_biall_status),],
            aes(fill=brca_biall_status, color=if(!is.null(df$in_training_set)){ in_training_set } else { TRUE }),
            shape=21
         ) +
         
         scale_color_manual(values=c('TRUE'='black','FALSE'='grey'), name='in_training_set') +
         scale_fill_manual(values=FILL_COLORS) +
         
         draw_grob(
            tableGrob(quadrant_counts, theme=ttheme_minimal(base_colour='darkred')),
            x=x.cutoff, y=y.cutoff, hjust=0.5, vjust=0.5,
         ) +
         
         ggtitle(
            paste0('n = ',nrow(df)),
            subtitle=if(is.null(subtitle)){ waiver() } else { subtitle }
         ) +
         
         theme_bw() +
         theme(
            legend.position='bottom',
            legend.direction='vertical'
         )
      
      if(is.null(df$in_training_set)){ p <- p + guides(color=F) }
      
      if(!is.null(axis.labs)){
         p <- p + xlab(axis.labs[1]) + ylab(axis.labs[2])
      }
      
      return(p)
   }
   
   #--------- CHORD v1 vs v2 ---------#
   pred_pairs <- lapply(preds, mkPredPairs)
   
   ## Reduction in radiotherapy samples
   calcTreatmentEnrichment <- function(df, by.col='has_radiotherapy_pre_treatment'){
      #df=pred_pairs$v2$clonal_subclonal
      
      df <- df[!is.na(df[[by.col]]),]
      
      nw <- table(subset(df, p_hrd.x<0.5 & p_hrd.y>=0.5)[[by.col]])
      sw <- table(subset(df, p_hrd.x<0.5 & p_hrd.y<0.5)[[by.col]])
      
      m <- rbind(
         c(nw['Yes'],sw['Yes']),
         c(sum(nw),sum(sw))
      )
      
      ratios <- m[1,]/m[2,]
      names(ratios) <- c('yes_nw','yes_sw')
      
      fisher <- unname(fisher.test(m,alternative='greater')$p.value)
      
      list(m=m, ratios=ratios, p.value=fisher)
   }
   
   calcTreatmentEnrichment(pred_pairs$v2$clonal_subclonal)
   calcTreatmentEnrichment(pred_pairs$v1$clonal_subclonal)
   
   calcTreatmentEnrichment(pred_pairs$v2$clonal_subclonal, 'has_systemic_pre_treatment')
   calcTreatmentEnrichment(pred_pairs$v1$clonal_subclonal, 'has_systemic_pre_treatment')
   
   calcTreatmentEnrichment(pred_pairs$v2$clonal_subclonal, 'has_any_pre_treatment')
   calcTreatmentEnrichment(pred_pairs$v1$clonal_subclonal, 'has_any_pre_treatment')
   
   p_v1_vs_v2 <- (function(){
      df <- mkPredPairs(
         list(v1=preds$v1$all, v2=preds$v2$all)
      )
      plotPredPair(df, axis.labs=c('p_hrd v1','p_hrd v2'), subtitle='HRD score v1 vs v2')
   })()
   
   #--------- Clonality preds ---------#
   # plotPredClonalityCombn <- function(l){
   #    #l=pred_pairs$v1
   #    lapply(names(l),function(i){
   #       #i=names(l)[[1]]
   #       axis_titles <- strsplit(i,'_')[[1]]
   #       plotPredPair(l[[i]], axis_titles)
   #    })
   #    #plot_grid(plotlist=plotlist, nrow=1, align='h',axis='tblr')
   # }
   
   p_clonality_combn <- list(
      plotPredPair(pred_pairs$v1$clonal_subclonal, axis.labs=c('clonal','subclonal'), subtitle='Subclonal vs clonal (v1)'),
      plotPredPair(pred_pairs$v2$clonal_subclonal, axis.labs=c('clonal','subclonal'), subtitle='Subclonal vs clonal (v2)'),
      
      plotPredPair(pred_pairs$v1$subclonal_all, swap.x.y=T,  axis.labs=c('all','subclonal'), subtitle='Subclonal vs all (v1)'),
      plotPredPair(pred_pairs$v2$subclonal_all, swap.x.y=T,  axis.labs=c('all','subclonal'), subtitle='Subclonal vs all (v2)')
   )
   
   #========= Output =========#
   out <- plot_grid(
      p_clonality_combn[[1]], p_clonality_combn[[2]], p_v1_vs_v2,
      p_clonality_combn[[3]], p_clonality_combn[[4]],
      
      nrow=2, align='h',axis='tblr'
   )
   
   if(!is.null(out_dir)){
      out_dir <- paste0(rf$dir,'/plots/')
      dir.create(out_dir, recursive=T, showWarnings=F)
      
      pdf(paste0(out_dir,'/chord_v1_vs_v2_preds_clonality.pdf') ,13,10)
      grid.draw(out)
      dev.off()
   } else {
      return(out)
   }

}

plotV1vsV2ComparisonWrapper <- function(rf){
   #rf <- chord_v2_models[['2.02']]
   out_dir <- paste0(rf$dir,'/plots/')
   dir.create(out_dir, recursive=T, showWarnings=F)
   
   p1 <- plotV1vsV2Comparison(rf, contexts_clonality_hmf, metadata_hmf)
   p2 <- plotV1vsV2Comparison(rf, contexts_clonality_pcawg, metadata_pcawg)
   
   pdf(paste0(out_dir,'/chord_v1_vs_v2_preds_clonality.pdf') ,13,12)
   grid.draw(p1)
   grid.draw(p2)
   dev.off()
}


#========= Models =========#
chord_v2_models <- list()

# #--------- CHORD v1 clone ---------#
# chord_v2_models[['2.00a']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.00a_CHORDv1_oldTrainSet/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types=c('snv','indel'),rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()

chord_v2_models[['2.00b']] <- (function(){
   dir <- paste0(base_dir,'/CHORDv2/training/models/2.00b_CHORDv1/')
   rf <- readRDS(paste0(dir,'/rf_model.rds'))
   rf$trans.func <- function(x){ transformContexts(x, simplify.types=c('snv','indel'),rel.types='all') }
   rf$dir <- dir
   return(rf)
})()

# #--------- Expanded indels ---------#
# chord_v2_models[['2.01a']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.01a_expandedIndels/')
#    #dir <- paste0(base_dir,'/CHORDv2/training/models/2.01a_expandedIndels/seed_models/seed04/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types='snv',rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()
# #plotV1vsV2ComparisonWrapper(chord_v2_models[['2.01a']])
# 
# chord_v2_models[['2.01b']] <- (function(){
#    dir <- paste0(base_dir,'/CHORDv2/training/models/2.01b_expandedIndels_wilcox1e-10/')
#    rf <- readRDS(paste0(dir,'/rf_model.rds'))
#    rf$trans.func <- function(x){ transformContexts(x, simplify.types='snv',rel.types='all') }
#    rf$dir <- dir
#    return(rf)
# })()
# #plotV1vsV2ComparisonWrapper(chord_v2_models[['2.01b']])

#--------- Merged del mh bimh 2-5 ---------#
chord_v2_models[['2.02']] <- (function(){
   dir <- paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/')
   #dir <- paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/seed_models/seed04/')
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
plotV1vsV2ComparisonWrapper(chord_v2_models[['2.02']])

# #--------- Sv mh ---------#
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
# 
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
# 
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
# 
# for(i in 1:length(chord_v2_models)){
#    #i=1
#    message('Processing: ',names(chord_v2_models)[i])
#    rf <- chord_v2_models[[i]]
#    out_dir <- paste0(rf$dir,'/plots/')
#    dir.create(out_dir, recursive=T, showWarnings=F)
#    
#    p1 <- plotV1vsV2Comparison(rf, contexts_clonality_hmf, metadata_hmf)
#    p2 <- plotV1vsV2Comparison(rf, contexts_clonality_pcawg, metadata_pcawg)
#    
#    pdf(paste0(out_dir,'/chord_v1_vs_v2_preds_clonality.pdf') ,13,6)
#    grid.draw(p1)
#    grid.draw(p2)
#    dev.off()
#    #if(i==1){ break }
# }





