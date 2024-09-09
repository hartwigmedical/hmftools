library(RColorBrewer)
library(reshape2)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(openxlsx)

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

#========= Load data =========#
pred <- read.delim(paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/pred_analysis_samples.txt'))
rank_order_clust <- readRDS(paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/rank_order_clust.rds'))
l_m_diplotypes <- readRDS(paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/l_m_diplotypes.rds'))
metadata <- read.delim(paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/metadata_for_analysis.txt'))


####################################################################################################
# Make plot data                                                                                   #
####################################################################################################
#========= Master dataframe =========#
summary_hrd <- (function(){
   
   selected_samples <- pred$sample
   
   #--------- Tag samples with genetic cause ---------#
   df <- data.frame(
      sample=names(rank_order_clust$clusters),
      cluster=rank_order_clust$clusters,
      is_hrd=TRUE,
      row.names=NULL
   )
   df$genotype <- names(rank_order_clust$cluster_names)[ match(df$cluster,rank_order_clust$cluster_names) ]
   
   df$genotype_monoall <- (function(){
      m <- rank_order_clust$df_ranked
      gene_cols <- grep('^p_',colnames(m),invert=T)
      apply(m[,gene_cols], 1, function(i){
         colnames(m[,gene_cols])[ which.max(i) ]
      })
   })()
   
   #--------- Add biallelic effect ---------#
   flattenDiplotypeMatrix <- function(m){
      #m=l_m_diplotypes$eff
      
      genes <- unique(sapply(strsplit(colnames(m),'_'),`[`,1))
      gene_cols <- lapply(genes, function(i){ 
         grep(i, colnames(m))
      }); names(gene_cols) <- genes
      
      ## Reorder eff matrix to match with df
      if(!all(df$sample %in% rownames(m))){
         warning("df$sample and rownames(m) are not completely intersecting")
      }
      
      m <- m[match(df$sample,rownames(m)),]
      
      ## Select diplotype (value pairs) for correct gene
      m$genotype <- df$genotype_monoall
      # m$genotype <- gsub('unknown_', '', m$genotype)
      # m$genotype <- gsub('-type', '', m$genotype)
      
      m_flat <- t(apply(as.matrix(m),1,function(i){
         genotype <- i['genotype']
         if(genotype %in% genes){
            i[ gene_cols[[genotype]] ]
         } else {
            c('none','none')
         }
      })); colnames(m_flat) <- c('a1','a2')
      
      return(m_flat)
   }
   
   m_diplotypes_flat <- (function(){
      m_eff <- flattenDiplotypeMatrix(l_m_diplotypes$eff)
      
      m_a_origin <- flattenDiplotypeMatrix(l_m_diplotypes$a_origin)
      colnames(m_a_origin) <- paste0(colnames(m_a_origin),'.origin')
      
      m_a_score <- flattenDiplotypeMatrix(l_m_diplotypes$score)
      colnames(m_a_score) <- paste0(colnames(m_a_score),'.score')
      
      data.frame(m_eff, m_a_origin, m_a_score)
   })()
   
   df$diplotype_origin <- apply(m_diplotypes_flat,1,function(i){
      a1 <- i[1]
      a2 <- i[2]
      a1.origin <- i[3]
      a2.origin <- i[4]
      
      
      if(a1=='deep_deletion'){
         a1
      } else if(a1=='loh'){
         if(a2!='none'){
            paste0('loh+',a2.origin)
         } else {
            'loh+unknown'
         }
      } else if(a1!='none' & a2!='none'){
         paste0(a1.origin,'+',a2.origin)
      } else {
         'unknown'
      }
   })
   
   ## Set diplotype origin to loh+unknown for samples in cluster 4 and 6
   df$diplotype_origin <- unlist(Map(function(genotype,diplotype_origin){
      if(grepl('^unknown',genotype) & diplotype_origin %in% c('loh+germ','loh+som') ){
         return('loh+unknown')
      } else {
         return(diplotype_origin)
      }
   },df$genotype,df$diplotype_origin))
   
   ## Format diplotype origin
   df <- within(df,{
      diplotype_origin <- gsub('[+]', ' + ', diplotype_origin)
      diplotype_origin <- gsub('deep_deletion', 'Deep deletion', diplotype_origin)
      diplotype_origin <- gsub('loh', 'LOH', diplotype_origin)
   })
   
   ## Format genotypes (pt.2)
   df$genotype <- gsub('unknown_BRCA1-type','unkn. (BRCA1-type HRD)',df$genotype)
   df$genotype <- gsub('unknown_BRCA2-type','unkn. (BRCA2-type HRD)',df$genotype)
   
   #--------- Add HRP samples to df ---------#
   df_prof <- df[1,]
   counter <- 0
   for(i in df){
      counter <- counter+1
      if(is.numeric(i)){ df_prof[counter] <- 0 }
      else if(is.logical(i)){ df_prof[counter] <- FALSE }
      else { df_prof[counter] <- 'NA' }
   }
   
   prof_samples <- selected_samples[ !(selected_samples %in% df$sample) ]
   df_prof <- df_prof[rep(1,length(prof_samples)),]
   df_prof$sample <- prof_samples
   
   df <- rbind(df, df_prof)
   
   if(!all(df$sample %in% selected_samples)){
      warning("df$sample and selected_samples are not completely intersecting")
   }
   
   #--------- Add metadata ---------#
   df$cancer_type <- metadata$cancer_type[ match(df$sample,metadata$sample) ]
   df$cohort <- metadata$cohort[ match(df$sample,metadata$sample) ]
   
   sample_names <- df$sample
   sample_names_converted <- metadata$hmf_id[ match(sample_names, metadata$sample) ]
   sample_names_converted[is.na(sample_names_converted)] <- sample_names[is.na(sample_names_converted)]
   
   df$sample <- sample_names_converted
   
   return(df)
})()
write.table(
   summary_hrd, paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/summary_hrd.txt'),
   sep='\t',quote=F,row.names=F
)


#========= Make plot data =========#
countFeatures <- function(df){
   
   #df=subset(summary_hrd, cohort=='HMF')

   #========= Main counting function =========#
   countFeature1ByFeature2 <- function(
      df, feature.split, feature.count, rm.rownames=T, 
      calc.rel.counts=NULL, 
      calc.unsplit.counts=F, unsplit.name='Pan-cancer'
   ){
      # df=summary_hrd
      # feature.split='cancer_type'
      # feature.count='is_hrd'
      # row.order=hrd_by_cancer_type$feature_split
      
      df_split <- split(df, df[,feature.split])
      df_split <- df_split[unique(df[,feature.split])] ## preserve original order
      
      ## Main
      values_uniq <- unique(df[,feature.count])
      values_uniq <- structure(rep(0,length(values_uniq)),names=values_uniq)
      out <- do.call(rbind, lapply(df_split, function(i){
         #i=df_split[[1]]
         tab <- table(i[,feature.count])
         v <- values_uniq
         v[names(tab)] <- tab
         return(v)
      }))
      
      out <- as.data.frame(out)
      
      ## Calculate pancancer counts
      if(calc.unsplit.counts){
         tab <- table(df[,feature.count])
         v <- values_uniq
         v[names(tab)] <- tab
         out <- rbind(v, out)
         rownames(out)[1] <- 
            if(is.null(unsplit.name)){ 'All' }
         else { unsplit.name }
      }
      
      if(!is.null(calc.rel.counts)){
         if(calc.rel.counts=='macro'){ out <- out/nrow(df) }
         if(calc.rel.counts=='micro'){ out <- out/rowSums(out) }
      }
      
      ## feature.split as column
      if(rm.rownames){
         out <- cbind(feature_split=rownames(out),out)
         rownames(out) <- NULL
      }
      
      return(out)
   }
   
   
   #========= HRD by cancer type =========#
   message('Calculating HRD freq per cancer type...')
   calcHrdByCancerType <- function(summary.hrd, min.hrd.freq=5){
      #summary.hrd=summary_hrd
      #summary.hrd=subset(summary_hrd, cluster %in% c(1,2,3,5))
      
      df <- as.data.frame(
         countFeature1ByFeature2(summary.hrd, 'cancer_type', 'is_hrd', rm.rownames=F, calc.unsplit.counts=T)
      )
      df$total <- df[['FALSE']] + df[['TRUE']]
      colnames(df)[1:2] <- c('abs','rel')
      df$rel <- df$abs/df$total
      
      df <- data.frame(feature_split=rownames(df), df)
      
      ## 
      df$split_var <- ifelse(df$abs <= min.hrd.freq, 'low', 'high')
      df$split_var[df$feature_split=='Pan-cancer'] <- 'except'
      
      ## Deal split sort cancer types with low/high absolute number of hrd samples
      df_split <- split(df, df$split_var)
      out <- do.call(rbind, lapply(df_split, function(i){
         i[order(i$rel, decreasing=T),]
      }))
      out$feature_split <- factor(out$feature_split,out$feature_split)
      out$split_var <- NULL
      
      ## Store metadata
      class(out) <- c(class(out), min.hrd.freq=min.hrd.freq)
      rownames(out) <- NULL
      
      return(out)
   }
   
   hrd_by_cancer_type <- calcHrdByCancerType(df)
   
   #========= Count features by cancer type =========#
   message('Counting HRD causes per cancer type...')
   #--------- Counts ---------#
   sel_features <- c('genotype','diplotype_origin')
   feature_counts <- lapply(sel_features, function(i){
      out <- countFeature1ByFeature2(
         df,
         feature.split='cancer_type',
         feature.count=i,
         calc.unsplit.counts=T
      )
      out <- out[ match(hrd_by_cancer_type$feature_split, out$feature_split), ]
      out$feature_split <- factor(out$feature_split,out$feature_split)
      out <- out[colnames(out)!='NA']
   })
   names(feature_counts) <- sel_features
   
   #--------- Simplify LOH + mut groups  ---------#
   feature_counts$biall_hit_type <- (function(){
      df <- feature_counts$diplotype_origin
      mut_mut_names <- c('germ + som','som + som')
      
      df[,'2x small mut.'] <- rowSums(df[,colnames(df) %in% mut_mut_names,drop=F])
      df <- df[,!(colnames(df) %in% mut_mut_names)]
      
      colnames(df)[3:4] <- c('LOH + somatic mut.','LOH + germline mut.')
      
      return(df)
   })()
   
   #--------- Biallelic hit origin  ---------#
   feature_counts$mut_origin <- (function(){
      df <- feature_counts$diplotype_origin
      
      germ_and_som <- rowSums(df[,colnames(df) %in% c('LOH + germ','germ + som'),drop=F])
      somatic_2x <- rowSums(df[,colnames(df) %in% c('Deep deletion','LOH + som','som + som'),drop=F])
      unknown <- rowSums(df[,colnames(df) %in% c('LOH + unknown','unknown'),drop=F])
      
      data.frame(
         feature_split=df$feature_split, 
         'Germline + somatic'=germ_and_som,
         '2x somatic'=somatic_2x, 
         'Unknown'=unknown, 
         check.names=F
      )
   })()

   #========= Export counts =========#
   #--------- Plot data ---------#
   message('Returning plot data...')
   
   feature_counts$diplotype_origin <- NULL
   
   out <- 
      c(
         hrd_freq=list(hrd_by_cancer_type),
         feature_counts
      )
      
   # saveRDS(
   #    c(
   #       hrd_freq=list(hrd_by_cancer_type),
   #       feature_counts
   #       #lapply(feature_counts, melt,'feature_split')
   #    ),
   #    paste0(data_out_dir,'/feature_counts.rds')
   # )
   
   # #--------- Supplementary data ---------#
   # message('Exporting supplementary table...')
   # calcRelCounts <- function(df){
   #    #df=feature_counts$biall_hit_type
   #    rel_counts <- df[,-1]/rowSums(df[,-1])
   #    cbind(feature_split=df[,1], rel_counts)
   # }
   # 
   # l_xlsx <- list()
   # 
   # df_export <- df
   # 
   # df_export$genotype_monoall <- NULL
   # 
   # colnames(df_export)[colnames(df_export)=='diplotype_origin'] <- 'hit_type'
   # 
   # l_xlsx[['Raw data']] <- df_export
   # 
   # l_xlsx[['Frequency of HRD']] <- hrd_by_cancer_type
   # 
   # l_xlsx[['Gene deficiency']] <- feature_counts$genotype
   # l_xlsx[['Gene deficiency (rel.)']] <- calcRelCounts(feature_counts$genotype)
   # 
   # l_xlsx[['Biall. hit type']] <- feature_counts$biall_hit_type
   # l_xlsx[['Biall. hit type (rel.)']] <- calcRelCounts(feature_counts$biall_hit_type)
   # 
   # l_xlsx[['Origin of biall. hit']] <- feature_counts$mut_origin
   # l_xlsx[['Origin of biall. hit (rel.)']] <- calcRelCounts(feature_counts$mut_origin)
   # 
   # formatXlsxTable <- function(df){
   #    # isNotFloatingPoint <- function(x){
   #    #    all(!grepl('^\\d+[.]',x))
   #    # }
   #    
   #    out <- as.data.frame(lapply(df,function(i){
   #       if(is.character(i)){ i }
   #       else if(is.numeric(i)){
   #          i[is.na(i)] <- 0
   #          i
   #       }
   #       # else if(is.numeric(i) & isNotFloatingPoint(i)){ as.character(i) }
   #       # else if(is.numeric(i)){ as.character(formatC(i, format='f', flag='0', digits=7)) }
   #       else { i }
   #    }), check.names=F)
   #    
   #    colnames(out)[colnames(out)=='feature_split'] <- 'Cancer type'
   #    return(out)
   # }
   # 
   # l_xlsx <- lapply(l_xlsx, formatXlsxTable)
   # 
   # write.xlsx(
   #    l_xlsx, paste0(data_out_dir,'/hrd_by_cancer_type.xlsx')
   # )
}

feature_counts <- list()

feature_counts$hmf <- countFeatures(subset(summary_hrd, cohort=='HMF'))
feature_counts$pcawg <- countFeatures(subset(summary_hrd, cohort=='PCAWG'))

saveRDS(
   feature_counts,
   paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/feature_counts.rds')
)

####################################################################################################
# Plotting                                                                                         #
####################################################################################################
#========= Prep data =========#
forceDfOrder <- function(df){
   as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
}

## For cancers with 0 counts for HMF or PCAWG, move up cancer types with counts >0 in either HMF
## or PCAWG
cancer_type_order <- (function(){
   df1 <- feature_counts$pcawg$hrd_freq[,c('feature_split','abs','rel')]
   df2 <- feature_counts$hmf$hrd_freq[,c('feature_split','abs','rel')]
   
   df <- merge(df1,df2,'feature_split',all=T)
   df[is.na(df)] <- 0
   
   min_hrd_freq <- 5
   df <- do.call(
      rbind,
      lapply(split(df, df$abs.x < min_hrd_freq & df$abs.y < min_hrd_freq), function(i){
         i[order(i$rel.x,i$rel.y, decreasing=T),]
      })
   )
   
   v <- as.character(df$feature_split)
   v <- v[v!='Pan-cancer']
   v <- c('Pan-cancer',v)
   return(v)
})()


#========= Make supplementary table =========#
feature_counts_xlsx <- (function(){
   feature_names <- names(feature_counts[[1]])
   
   l_abs <- lapply(feature_names, function(i){
      #i=1
      
      l_df <- lapply(c('hmf','pcawg'), function(j){
         #j='pcawg'
         df <- feature_counts[[j]][[i]]
         colnames(df)[1] <- 'cancer_type'
         
         ## Fill cancer types with no data
         df$cancer_type <- as.character(df$cancer_type)
         rownames(df) <- df$cancer_type
         
         missing_cancer_types <- !(cancer_type_order %in% df$cancer_type)
         df <- df[cancer_type_order,]
         df$cancer_type[missing_cancer_types] <- cancer_type_order[missing_cancer_types]
         rownames(df) <- NULL
         
         df$cancer_type <- factor(df$cancer_type, unique(df$cancer_type))
         
         ## Add cohort name
         df <- cbind(cancer_type=df$cancer_type, cohort=toupper(j), df[,-1])
         
         return(df)
      })
      
      out <- do.call(rbind, l_df)
      out <- out[order(out$cancer_type),]
      
      return(out)
   })
   names(l_abs) <- feature_names
   
   l_rel <- lapply(l_abs, function(i){
      #i=l_abs[[1]]
      df <- as.data.frame(i[,-c(1:2)])
      
      missing_cancer_types <- is.na(df[,1])
      
      df <- df/rowSums(df)
      df[is.na(df)] <- 0
      
      df <- as.data.frame(lapply(df,function(j){
         j[missing_cancer_types] <- NA
         return(j)
      }), check.names=F)
      
      cbind(i[,1:2], df)
   })
   
   l_merged <- list()
   l_merged$raw <- summary_hrd
   
   for(i in names(l_abs)){
      #i='genotype'
      l_merged[[paste0(i,'.abs')]] <- l_abs[[i]]
      l_merged[[paste0(i,'.rel')]] <- l_rel[[i]]
   }
   
   l_merged$hrd_freq.rel <- NULL
   names(l_merged)[names(l_merged)=='hrd_freq.abs'] <- 'hrd_freq'
   
   return(l_merged)
})()

write.xlsx(
   feature_counts_xlsx, paste0(base_dir,'/CHORDv2/processed/analysis/pancancer_overview/data/hrd_by_cancer_type.xlsx')
)

## Restructure list
plot_data <- (function(){
   listnames <- names(feature_counts[[1]])
   
   l <- lapply(listnames, function(i){
      #i='hrd_freq'
      #print(i)
      
      df1 <- feature_counts$hmf[[i]]
      df2 <- feature_counts$pcawg[[i]]
      
      if(i!='hrd_freq'){
         df1 <- melt(df1, 'feature_split')
         df2 <- melt(df2, 'feature_split')
      }
      
      df1$cohort <- 'HMF'
      df2$cohort <- 'PCAWG'
      
      out <- rbind(df1,df2)
      
      colnames(out)[1] <- 'cancer_type'
      
      out <- forceDfOrder(out)
      out$cancer_type <- factor(out$cancer_type, levels=cancer_type_order)
      out$cohort <- factor(out$cohort, levels=c('PCAWG','HMF'))
      
      return(out)
   })
   
   names(l) <- listnames
   return(l)
})()

#========= Plotting constants =========#
THEME_DEFAULT <- theme(
   panel.spacing=unit(0,'pt'),
   
   panel.grid.major.x=element_blank(),
   panel.grid.minor.x=element_blank(),
   panel.grid.major.y=element_blank(),
   panel.grid.minor.y=element_blank(),
   axis.title.x=element_blank(),
   axis.text.x=element_blank(),
   axis.ticks.x=element_blank(),
   strip.text.x=element_text(angle=90, hjust=1, vjust=0.5),
   
   axis.title.y=element_text(size=11),
   
   panel.border=element_rect(color='darkgrey'),
   strip.background=element_rect(color='grey',fill='#e9e9e9'),
   
   legend.title=element_text(size=0),
   legend.spacing.x=unit(0.1,'cm'),
   legend.key.size=unit(0.6,'line'),
   legend.text=element_text(size=9)
)

## Order of elements will be enforced in the plots
GROUP_FILLS <- list()

GROUP_FILLS$genotype <- list(
   BRCA1_type=c(
      'BRCA1'='#2C867C', ## Teal
      'unkn. (BRCA1-type HRD)'='#BCE6DF' ## Light teal
   ),
   
   BRCA2_type=c(
      'unkn. (BRCA2-type HRD)'='#f7eed3', ## pale brown
      'PALB2'='#D7B66A', ## light brown
      'RAD51C'='#B06E23', #'#b04b23', ## brown
      'BRCA2'='#6c3809' #'#783f0b' ## Dark brown
   )
)

GROUP_FILLS$biall_hit_type <- c(
   'unknown'='white',
   '2x small mut.'='#FEF7C2', ## light yellow
   
   'LOH + unknown'='#d1c0df',
   'LOH + germline mut.'='#8459ab', ## light purple
   'LOH + somatic mut.'='#512678', ## purple
   
   'Deep deletion'='#ff42a8' ##pink
)

GROUP_FILLS$mut_origin <- c(
   Unknown='white',
   `2x somatic`='#8DB7D3',
   `Germline + somatic`='#3163A0'
)

#========= HRD freq =========#
plotFreqHrdByCancerType <- function(
   df, show.labels=T, strip.position='bottom', hide.x.strips=F
){
   # df=plot_data$hrd_freq
   # strip.position='bottom'
   # hide.x.strips=F
   # show.labels='perc'
   
   if(show.labels != FALSE){
      df$label <- with(df,{
         perc <- paste0( signif(100*rel,2), '%' )
         
         if(show.labels==TRUE | show.labels=='full'){
            label <- paste0(perc, ' ', abs,' / ',total)   
         } else if(show.labels=='perc'){
            label <- perc
         } else if(show.labels=='ratio'){
            label <- paste0(abs,' / ',total)
         }
         #label[rel==0] <- ''
         return(label)
      })
   }
   
   p <- ggplot(df, aes(x=cohort, y=rel)) +
      facet_wrap(.~cancer_type, nrow=1, strip.position=strip.position) +
      
      geom_bar(aes(fill=cohort), stat='identity') +
      scale_fill_manual(name='Cohort', values=c('#96c697','#51a152'), labels=c('PCAWG (primary)','HMF (metastatic)')) +
      #scale_fill_manual(name='Cohort', values=c('#96c697','#abb1b9'), labels=c('PCAWG (primary)','HMF (metastatic)')) +
      
      scale_y_continuous(
         name='Freq. HRD',
         labels=function(x){ paste0(x*100,'%') }
      ) +
      
      theme_bw() +
      THEME_DEFAULT
   
   if(show.labels != FALSE){
      label_baseline <- max(df$rel)*0.02
      p <- p + geom_text(aes(y=label_baseline,label=label), size=2.5, angle=90, hjust=0, vjust=0.5)
   }
   
   if(strip.position=='top'){
      p <- p + theme(strip.text.x=element_text(angle=90, hjust=0, vjust=0.5))
   }
   
   if(hide.x.strips){
      p <- p + theme(
         strip.background=element_blank(),
         strip.text.x=element_blank()
      )
   }
   
   return(p)
}

#========= Genetic causes =========#
plotFeatureByCancerType <- function(
   df, fill.colors=NULL, 
   y.axis.title=NULL, strip.position='bottom', hide.x.strips=F
){
   # df=plot_data$genotype
   # fill.colors=GROUP_FILLS$genotype
   # y.axis.title='Gene deficiency'
   # strip.position='top'
   # hide.x.strips=F
   
   if(!is.null(fill.colors)){
      fill_colors <- if(is.list(fill.colors)){ unlist(unname(fill.colors)) } else { fill.colors }
      df$variable <- factor(df$variable, levels=names(fill_colors))
   }
   
   p <- ggplot(df, aes(x=cohort, y=value)) +
      facet_wrap(.~cancer_type, nrow=1, strip.position=strip.position) +
      
      geom_bar(
         aes(fill=variable), stat='identity', position='fill', 
         width=1, color='black', size=0.2
      ) +
      scale_fill_manual(values=fill_colors) +
      scale_x_discrete(expand=c(0.6, 0.6)) +
      
      scale_y_continuous(
         name=if(!is.null(y.axis.title)){ y.axis.title } else { waiver() },
         labels=function(x){ paste0(x*100,'%') }
      ) +
      
      theme_bw() +
      THEME_DEFAULT
   
   if(strip.position=='top'){
      p <- p + theme(strip.text.x=element_text(angle=90, hjust=0, vjust=0.5))
   }
   
   if(hide.x.strips){
      p <- p + theme(
         strip.background=element_blank(),
         strip.text.x=element_blank()
      )
   }
   
   return(p)
}

p_overview <- (function(){
   
   l <- list()
   
   l$hrd_freq <- 
      plotFreqHrdByCancerType(
         plot_data$hrd_freq, show.labels=T, strip.position='top'
      )
   
   l$genotype <- (function(){
      p <- plotFeatureByCancerType(
         plot_data$genotype,
         fill.colors=GROUP_FILLS$genotype,
         y.axis.title='Gene deficiency',
         hide.x.strips=T
      )
      
      ## Hacky solution to split up BRCA2/BRCA1-type HRD keys and add title
      fill_pal <- with(GROUP_FILLS$genotype, { 
         c(rev(BRCA2_type), rev(BRCA1_type))
      })
      
      legend_breaks <- with(GROUP_FILLS$genotype,{
         names(c(rev(BRCA2_type),BRCA1_type))
      })
      
      legend_labels <- legend_breaks
      legend_labels[grep('^unkn',legend_labels)] <- 'Unknown             '
      
      p <- p + 
         scale_fill_manual(
            breaks=legend_breaks,
            labels=legend_labels,
            values=fill_pal
         ) +
         guides(fill=guide_legend(title='BRCA2-type HRD   BRCA1-type HRD',ncol=2, nrow=4)) +
         theme(legend.title=element_text(hjust=0, size=9, face='bold'))
      
      return(p)
   })()
   
   l$biall_hit_type <- 
      plotFeatureByCancerType(
         plot_data$biall_hit_type,
         fill.colors=GROUP_FILLS$biall_hit_type,
         y.axis.title='Biallelic hit type',
         hide.x.strips=T
      )
   
   l$mut_origin <- 
      plotFeatureByCancerType(
         plot_data$mut_origin,
         fill.colors=GROUP_FILLS$mut_origin,
         y.axis.title='Origin of biallelic hit',
         hide.x.strips=F
      )
   
   ## Adjust legend position inside plot
   l <- lapply(l, function(p){
      p + theme(
         legend.position=c(0.68, 0.5), 
         legend.justification=c(0, 0.5),
         legend.background=element_rect(fill='#f1f1f1', color='lightgrey')
      )
   })
   l$hrd_freq <- l$hrd_freq + theme(legend.position=c(0.68, 0.7))
   
   plot_grid(
      plotlist=l, align='v', axis='tblr', ncol=1,
      rel_heights=c(1.7, 1, 1, 1.7)
   )
})()

pdf(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/plots/hrd_cause_by_cancer_type.pdf'),9,9)
grid.draw(p_overview)
dev.off()


#========= Custom fisher tests HMF vs PCAWG=========#
(function(){
   #--------- HRD freq ---------#
   ## Ovary
   fisher.test(
      rbind(
         c(42,34),
         c(142,65)
      ),
      alternative='less'
   )
   
   ## Breast
   fisher.test(
      rbind(
         c(75,28),
         c(636,118)
      ),
      alternative='less'
   )
   
   ## Prostate
   fisher.test(
      rbind(
         c(41,13),
         c(314,232)
      ),
      alternative='greater'
   )
   
   ## Pancreas
   fisher.test(
      rbind(
         c(11,21),
         c(83,287)
      ),
      alternative='greater'
   )
   
   
   #--------- Gene def ---------#
   ## BRCA2-type HRD, prostate
   fisher.test(
      rbind(
         c(39,5),
         c(41,13)
      ),
      alternative='greater'
   )
   
   ## BRCA1-type HRD, ovarian
   fisher.test(
      rbind(
         c(20,26),
         c(42,34)
      ),
      alternative='less'
   )
   
   ## BRCA1-type HRD, breast
   fisher.test(
      rbind(
         c(30,15),
         c(75,28)
      ),
      alternative='greater'
   )
   
   #--------- Biallelic hit type ---------#
   ## Pancancer deep dels
   fisher.test(
      rbind(
         c(26,6),
         c(206,104)
      ),
      alternative='greater'
   )
   
   #--------- Hit origin ---------#
   ## Prostate, pure somatic
   fisher.test(
      rbind(
         c(26,3),
         c(41,13)
      ),
      alternative='greater'
   )
})()


####################################################################################################
# Cancer types, biallelic loss clusters vs. non-biallelic loss clusters                            #
####################################################################################################
cancer_type_order <- levels(plot_data$hrd_freq$cancer_type)
cancer_type_order <- cancer_type_order[cancer_type_order!='Pan-cancer']

p_cancer_type_counts_split <- (function(){
   l <- list(
      `Biallelic loss\n(clusters: 1, 2, 3, 5)`=table(subset(summary_hrd, cluster %in% c(1,2,3,5))$cancer_type),
      `Non-biallelic loss\n(clusters: 4, 6)`=table(subset(summary_hrd, cluster %in% c(4,6))$cancer_type)
   )

   l <- lapply(l, function(i){
      #i <- l[[1]]
      missing_cancer_types <- unique(summary_hrd$cancer_type)[ !(unique(summary_hrd$cancer_type) %in% names(i)) ]
      i[missing_cancer_types] <- 0
      as.table(i)
   })

   pd <- do.call(rbind, lapply(names(l), function(i){
      #i <- names(l)[[1]]
      df <- data.frame(
         cancer_type=names(l[[i]]),
         abs=as.vector(l[[i]])
      )

      df <- df[match(cancer_type_order, df$cancer_type),]

      df$cancer_type <- as.character(df$cancer_type)
      df$cancer_type[!(df$cancer_type %in% c('Ovary','Pancreas','Prostate','Breast','Biliary','Urinary tract'))] <- 'Other'

      df_split <- split(df, df$cancer_type=='Other')

      df_split[['TRUE']] <- data.frame(
         cancer_type='Other',
         abs=sum(df_split[['TRUE']]$abs)
      )

      df <- do.call(rbind, df_split)

      df$rel <- df$abs/sum(df$abs)
      df$label <- paste0(
         df$abs,'\n',
         '(',signif(df$rel*100,2),'%)'
      )

      df$group <- i


      return(df)
   }))

   pd <- as.data.frame(lapply(pd, function(i){
      if(is.numeric(i)){ i }
      else { factor(i, unique(i)) }
   }))
   
   ggplot(pd, aes(cancer_type, abs, group=group, label=label)) +
      geom_bar(stat='identity') +
      geom_text(vjust=-0.5, size=2.5) +
      facet_grid(group~.) +
      ylab('Number of patients') +
      ylim(0, max(pd$abs)*1.3) +
      theme_bw() +
      theme(
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         axis.title.x=element_blank()
      )
   
})()


pdf(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/plots/cancer_type_counts_yes_vs_no_biall_loss.pdf'),7,5)
grid.draw(p_cancer_type_counts_split)
dev.off()
