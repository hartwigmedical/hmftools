library(reshape2)
library(cowplot)#; theme_set(theme_grey())
library(ggplot2)
library(grid)
library(RColorBrewer)
library(openxlsx)

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

script_dir <- paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/')

#========= Common data =========#
#--------- Other data ---------#
eff_metadata <- read.delim(paste0(script_dir,'/eff_metadata.txt'))
hr_gene_likelihood <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/hr_gene_likelihood.txt'))

metadata <- read.delim(paste0(script_dir,'/data/metadata_for_analysis.txt'))

#--------- CHORD ---------#
chord_dir <- paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/')

chord_features <- (function(){
  chord <- readRDS(paste0(chord_dir,'/rf_model.rds'))
  names(chord$forest$xlevels)
})()

pred <- read.delim(paste0(script_dir,'/data/pred_analysis_samples.txt'))
pred_hrd <- pred[pred$hr_status=='HR_deficient',]

#--------- Genes ---------#
genes_bed <- read.delim(paste0(base_dir,'/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/genes_chord_paper.bed'))
colnames(genes_bed)[1] <- 'chrom'

hr_genes_extended <- (function(){
  v <- genes_bed[genes_bed$is_hr_gene,'ensembl_gene_id']
  
  #v1 <- c('BRCA2','BRCA1','RAD51C','PALB2')
  v1 <- subset(genes_bed, hgnc_symbol %in% c('BRCA2','BRCA1','RAD51C','PALB2'))$ensembl_gene_id
  v2 <- sort(v[!(v %in% v1)])
  
  #c(v1,v2)
  v2
})() 

#--------- Diplotypes ---------#
diplotypes_hrd_hrGenes <- read.delim(paste0(script_dir,'/data/diplotypes_hrd_hrGenes.txt.gz'))

diplotypes_filt_ss <- (function(){
  path <- '/Users/lnguyen/Documents/R_cache/diplotypes_filt_ss.txt.gz'

  if(!file.exists(path)){
    df <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/diplotypes_filt.txt.gz'))
    df <- df[df$ensembl_gene_id %in% hr_genes_extended & df$sample %in% diplotypes_hrd_hrGenes$sample,]
    
    write.table(df, gzfile(path), sep='\t', quote=F, row.names=F)
  } else {
    df <- read.delim(path)
  }

  return(df)
})()

#========= Cohort specific data =========#
# #--------- HMF ---------#
# hmf_data <- list()
# 
# hmf_data$diplotypes_hrd_hrdGenes <- read.delim(paste0(script_dir,'/data/diplotypes_hrd_hrGenes.txt.gz'))
# hmf_data$metadata <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_pred_ann.txt.gz'))
# 
# hmf_data$pred_hrd <- subset(preds$HMF, hr_status=='HR_deficient' & sample %in% hmf_data$diplotypes_hrd_hrdGenes$sample)
# 
# hmf_data$diplotypes_filt_ss <- (function(){
#   path <- '/Users/lnguyen/Documents/R_cache/diplotypes_filt.txt.gz'
#   
#   if(!file.exists(path)){
#     file.copy(
#       paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/diplotypes_filt.txt.gz'),
#       path
#     )
#   }
#   
#   df <- read.delim(path)
#   df <- df[df$ensembl_gene_id %in% hr_genes_extended  & df$sample %in% hmf_data$diplotypes_hrd_hrdGenes$sample,]
#   
#   return(df)
# })()
# 
# #--------- PCAWG ---------#
# pcawg_data <- list()
# 
# pcawg_data$diplotypes_hrd_hrdGenes <- read.delim(paste0(script_dir,'/data_pcawg/diplotypes_hrd_hrGenes.txt.gz'))
# pcawg_data$metadata <- read.delim(paste0(base_dir,'/datasets/PCAWG_2020/metadata/pcawg_metadata_ann.txt.gz'))
# colnames(pcawg_data$metadata)[colnames(pcawg_data$metadata)=='cancer_type'] <- 'primary_tumor_location'
# 
# pcawg_data$pred_hrd <- subset(
#   preds$PCAWG, 
#   hr_status=='HR_deficient' 
#   & hrd_type!='cannot_be_determined'
#   & sample %in% pcawg_data$diplotypes_hrd_hrdGenes$sample
# )
# 
# pcawg_data$diplotypes_filt_ss <- (function(){
#   path <- '/Users/lnguyen/Documents/R_cache/diplotypes_filt_pcawg.txt.gz'
#   
#   if(!file.exists(path)){
#     file.copy(
#       paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data_pcawg/diplotypes_filt.txt.gz'),
#       path
#     )
#   }
#   
#   df <- read.delim(path)
#   df <- df[df$ensembl_gene_id %in% hr_genes_extended  & df$sample %in% pcawg_data$diplotypes_hrd_hrdGenes$sample,]
#   
#   return(df)
# })()




#========= Main =========#
main(
  diplotypes_hrd_hrGenes, pred_hrd, metadata, diplotypes_filt_ss, pred,
  data_out_dir=paste0(script_dir,'/data/'), 
  plot_out_dir=paste0(script_dir,'/plots/'), 
  convert_cpct_ids=T, 
  pcawg_apply_custom_hrd_type=T
){
  
  ## HMF 
  # diplotypes_hrd_hrGenes=hmf_data$diplotypes_hrd_hrdGenes
  # pred_hrd=hmf_data$pred_hrd
  # metadata=hmf_data$metadata
  # diplotypes_filt_ss=hmf_data$diplotypes_filt_ss
  # pred=hmf_data$pred_hrd
  # convert_cpct_ids=T
  # pcawg_apply_custom_hrd_type=F
  # 
  # data_out_dir=paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/')
  # plot_out_dir=paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/plots/')
  
  ## PCAWG
  # diplotypes_hrd_hrGenes=pcawg_data$diplotypes_hrd_hrdGenes
  # pred_hrd=pcawg_data$pred_hrd
  # metadata=pcawg_data$metadata
  # diplotypes_filt_ss=pcawg_data$diplotypes_filt_ss
  # pred=pcawg_data$pred_hrd
  # convert_cpct_ids=F
  # pcawg_apply_custom_hrd_type=T
  # 
  # data_out_dir=paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data_pcawg/')
  # plot_out_dir=paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/plots_pcawg/')

  
  ####################################################################################################
  # Make diplotype matrix and clustering                                                             #
  ####################################################################################################
  
  #========= Put data into wide format =========#
  message('Making heat map layers; Converting diplotypes from long to wide format...')
  mkDiplotypeMatrix <- function(
    
    ## Main
    diplotypes, output='eff', gene.list=NULL, ## Make sure to order genes by descending impact
    
    ## Options
    floor.values.where.score.eq0=T,
    boost.high.impact.eff.scores=T,
    simplify.max.score.origin=T,
    simplify.eff=T, snpeff.scoring=eff_metadata
  ){
    #diplotypes=diplotypes_filt_ss
    #gene.list=c('BRCA2','BRCA1','PALB2','RAD51C')
    #output='chrom'
    
    #--------- Init ---------#
    COMMON_COLS <- c('sample','hgnc_symbol')
    
    if(output=='eff' | is.na(output) | is.null(output)){
      #NA_FILL <- 'none'
      
      diplotypes_ss <- diplotypes[c(COMMON_COLS,'a1.eff','a2.eff')]
      
      if(simplify.eff){
        counter <- 1
        for(i in snpeff.scoring$ann){
          diplotypes_ss[diplotypes_ss==i] <- snpeff.scoring$ann_s1[counter]
          counter <- counter + 1
        }
      }
      
    } else if(output=='score'){
      #NA_FILL <- 0
      
      diplotypes_ss <- diplotypes #[c(COMMON_COLS,'a1.max_score','a2.max_score','a1','a2')]
      
      if(boost.high.impact.eff.scores){
        diplotypes_ss[diplotypes_ss$a1.eff=='deep_deletion',c('a1.max_score','a2.max_score')] <- 5.3
        diplotypes_ss[diplotypes_ss$a1.eff=='loh','a1.max_score'] <- 5.2
        diplotypes_ss[diplotypes_ss$a2.eff=='frameshift_variant','a2.max_score'] <- 5.1
      }
      
      diplotypes_ss <- diplotypes_ss[c(COMMON_COLS,'a1.max_score','a2.max_score')]
      
    } else if(output=='a_origin'){
      #NA_FILL <- 'none'
      
      diplotypes_ss <- diplotypes[c(COMMON_COLS,'diplotype_origin')]
      
      ## Split diplotype_origin column into 2
      a_origin <- do.call(rbind,strsplit(diplotypes_ss$diplotype_origin,'_'))
      colnames(a_origin) <- c('a1.origin','a2.origin')
      diplotypes_ss <- cbind(diplotypes_ss, a_origin)
      diplotypes_ss$diplotype_origin <- NULL
      
    } else if(output=='max_score_origin'){
      #NA_FILL <- 'none'
      
      diplotypes_ss <- diplotypes
      
      if(simplify.max.score.origin){
        diplotypes_ss <- within(diplotypes_ss,{
          a1.max_score_origin[grepl('clinvar|enigma',a1.max_score_origin)] <- 'known'
          a2.max_score_origin[grepl('clinvar|enigma',a2.max_score_origin)] <- 'known'
          
          a1.max_score_origin[a1.max_score < 5] <- 'none'
          a2.max_score_origin[a2.max_score < 5] <- 'none'
          
          a1.max_score_origin[grepl('^snpeff$',a1.max_score_origin)] <- 'predicted'
          a2.max_score_origin[grepl('^snpeff$',a2.max_score_origin)] <- 'predicted'
        })
      }
      
      diplotypes_ss <- diplotypes_ss[c(COMMON_COLS,'a1.max_score_origin','a2.max_score_origin')]
      
    } else {
      #NA_FILL <- 'none'
      
      sel_cols <- c(paste0('a1.',output), paste0('a2.',output))
      
      if(any(!(sel_cols %in% colnames(diplotypes)))){
        stop("Please specify valid output type: 'eff', 'score', 'a_origin', 'max_score_origin', or existing allele column suffixes")
      }
      
      diplotypes_ss <- diplotypes[c(COMMON_COLS, sel_cols)]
    }
    
    NA_FILL <- 
      if(is.numeric(diplotypes_ss[,3])){
        0
      } else {
        'none'
      }
    
    #--------- Floor values where a1/a2.max_score == 0 (due to PON filtering) ---------#
    if(floor.values.where.score.eq0){
      diplotypes_ss[diplotypes$a1.max_score==0, 3] <- NA_FILL
      diplotypes_ss[diplotypes$a2.max_score==0, 4] <- NA_FILL
    }
    
    #--------- Gene subset ---------#
    genes <- if(is.null(gene.list)){
      #warning("gene.list was not provided; gene impact order will not be considered. Defaulting to unique(hgnc_symbol). ")
      sort(unique(diplotypes$hgnc_symbol))
    } else {
      gene.list
    }
    
    #--------- Convert from long to wide table ---------#
    #gene_id_type <- if(grepl('^ENSG',genes[1])){ 'ensembl_gene_id' } else { 'hgnc_symbol' }
    
    l <- lapply(genes, function(i){
      #i='PTEN'
      df <- diplotypes_ss[diplotypes_ss[['hgnc_symbol']]==i,]
      colnames(df)[3:4] <- paste0(i,'_',1:2)
      
      df$hgnc_symbol <- NULL
      
      return(df)
    })
    
    out <- Reduce(function(x,y){ merge(x,y,by='sample',all.x=T, sort=F) }, l)
    
    #--------- Move sample col to rownames ---------#
    ## Move sample col to rownames
    rownames(out) <- out$sample
    out$sample <- NULL
    
    ## Deal with NA values created during merging
    out[is.na(out)] <- NA_FILL 
    
    return(out)
  }
  
  mkDiplotypeMatrixLayers <- function(diplotypes, ...){
    output_types <- c('eff','score','a_origin','max_score_origin','chrom','pos','hgvs_c')
    l <- lapply(output_types, function(i){
      m <- mkDiplotypeMatrix(diplotypes, output=i,...)
      #m <- m[match(rownames(pred),rownames(m)),]
      return(m)
    })
    names(l) <- output_types
    return(l)
  }
  
  l_m_diplotypes <- mkDiplotypeMatrixLayers(
    diplotypes_hrd_hrGenes, 
    #gene.list=c('BRCA2','BRCA1','RAD51C','PALB2')
    gene.list=with(hr_gene_likelihood,{ hgnc_symbol[keep==T] })
  )
  
  #========= Custom rank order clustering like method =========#
  message('Performing rank order clustering...')
  #--------- Functions ---------#
  mkGeneRankingMatrix <- function(diplotypes.score, min.biall.score=c(5,3), min.summed.score=6){
    #diplotypes.score=l_m_diplotypes$score
    col_names <- colnames(diplotypes.score)
    genes <- unique(gsub('_1|_2','',col_names))
    
    ## Add allele 1 score to allele 2 score
    out <- do.call(cbind,lapply(genes, function(i){
      df_ss <- diplotypes.score[,grep(i,col_names)]
      
      ## Filter out low hit score to get a clean ordering
      df_ss[,1][ df_ss[,1] < min.biall.score[1] ] <- 0
      df_ss[,2][ df_ss[,2] < min.biall.score[2] ] <- 0
      
      ## Calculate summed score
      df_ss[,1] + df_ss[,2]
    }))
    
    ## Filter out low hit score to get a clean ordering
    out[out <= min.summed.score] <- 0
    colnames(out) <- genes
    rownames(out) <- rownames(diplotypes.score)
    out <- as.data.frame(out)
    
    return(out)
  }
  
  mkGeneRankingMatrixChordPred <- function(df){
    #df=pred
    df_ss <- df[c('sample','BRCA2','BRCA1','hrd')] ## Put BRCA2 has more impact than BRCA1 base on fisher test
    rownames(df_ss) <- df_ss$sample
    
    df_ss[c('sample','hrd')] <- NULL
    
    out <- t(apply(df_ss,1, function(i){
      max_value <- max(i)
      as.integer(i==max_value)
    }))
    colnames(out) <- paste0('p_', colnames(df_ss))
    out <- as.data.frame(out)
    return(out)
  }
  
  rle2Clusters <- function(rle.out, names=NULL){
    counter <- 0
    clusters <- unlist(lapply(rle.out$lengths, function(i){
      counter <<- counter + 1
      rep(counter, i)
    }))
    names(clusters) <- names(names)
    return(clusters)
  }
  
  rankOrderClust <- function(df, force.1.value.per.row=F){
    
    #df=mkGeneRankingMatrix(l_m_diplotypes$score)
    #df=m_ss
    #df=m_score
    #df=m_groups
    
    ## Force single value per row
    if(force.1.value.per.row){
      ncols <- ncol(df)
      df <- as.data.frame(t(apply(df,1,function(i){
        if(all(i==0) || sum(i==0) == ncols-1){ ## Is all zero, or has single value
          return(i)
        } else {
          max_index <- which.max(i)
          i[-max_index] <- 0
          return(i)
        }
      })))
    }
    
    ## Get order
    #sample_order <- do.call(order, cbind(col_max,df))
    sample_order <- do.call(order, df)
    sample_order <- rev(sample_order)
    
    ## Apply ordering
    df_ranked <- df[sample_order,]
    #df_ranked$index <- 1:nrow(df_ranked)
    
    ## Determine cluster groups
    col_max <- apply(df_ranked,1,which.max)
    cluster_breaks <- rle(col_max)
    
    clusters <- rle2Clusters(cluster_breaks, names=col_max)
    
    list(
      df_ranked=df_ranked,
      clusters=clusters,
      order=rownames(df_ranked)
    )
  }
  
  mergeRankOrderClustOutput <- function(l.roclust){
    #l.roclust <- l_roclust
    
    df_ranked <- do.call(
      rbind, 
      unname(lapply(l.roclust,`[[`,'df_ranked')) 
    )
    
    ## Ensure if 2 list items have '1' as a cluster, that these won't be merged
    counter <- 0
    l_clusters <- lapply(l.roclust, function(i){ 
      counter<<-counter+1
      paste0(counter,i$cluster)
    })
    
    clusters <- rle2Clusters(
      rle(unlist(l_clusters, use.names=F))
    )
    names(clusters) <- rownames(df_ranked)
    
    list(df_ranked=df_ranked, clusters=clusters, order=rownames(df_ranked))
  }
  
  rank_order_clust <- (function(){
    
    ## Rank order by CHORD prediction, then hit score data
    m_pred <- mkGeneRankingMatrixChordPred(pred_hrd)
    
    if(pcawg_apply_custom_hrd_type){
      ## Force these samples to be BRCA2 deficient. Based on genetics, these samples were HRD type
      ## misclassified
      sel_samples <- c(
        '81bc7f0c-865d-4801-a935-2ab04170df53',
        '8282283d-247a-431d-9421-0fcc52f0a897',
        '53534b3c-cd15-4d68-a9b1-6902bb234c45'
      )
      m_pred[rownames(m_pred) %in% sel_samples,'p_BRCA1'] <- 0
      m_pred[rownames(m_pred) %in% sel_samples,'p_BRCA2'] <- 1
      #subset(pred_hrd,sample %in% sel_samples)
    }
    
    
    m_score_low <- mkGeneRankingMatrix(l_m_diplotypes$score, min.biall.score=c(5,0), min.summed.score=0)
    m_score_filt <- mkGeneRankingMatrix(l_m_diplotypes$score, min.biall.score=c(5,3), min.summed.score=6)
    
    #--------- 1st round of rank ordering ---------#
    ## Split into BRCA2/BRCA1 groups and further into high/low biall score
    m_groups <- merge(
      m_pred,
      data.frame(is_not_all_0 = apply(m_score_filt, 1, function(i){ !all(i==0) })),
      by='row.names'
    )
    m_groups$is_not_all_0 <- as.integer(m_groups$is_not_all_0)
    rownames(m_groups) <- m_groups[,1]; m_groups[,1] <- NULL
    m_groups <- rankOrderClust(m_groups)$df_ranked
    
    group_strings <- apply(m_groups,1,paste, collapse='')
    
    group_names <- unique(group_strings)
    names(group_names) <- apply(unique(m_groups),1,function(i){
      paste( names(i)[as.logical(i)], collapse='.' )
    })
    
    hrd_groups <- list(
      'BRCA2'=c('BRCA2','RAD51C','PALB2'),
      'BRCA1'=c('BRCA1')
    )
    
    #--------- 2nd round of rank ordering on diplotype scores ---------#
    counter <- 0
    l_roclust <- lapply(group_names, function(i){
      counter <<- counter+1
      #print(counter)
      # i=group_names[2]
      # counter=2
      group_samples <- names(group_strings)[group_strings==i]
      group_name <- names(group_names[counter])
      
      is_definite_biall <- substr(i,3,3)=='1'
      
      if(is_definite_biall){
        m_score <- m_score_filt
      } else {
        m_score <- m_score_low
        for(j in names(hrd_groups)){
          if(grepl(j,group_name)){
            hrd_group <- hrd_groups[[j]]
            break
          }
        }
        m_score[!(colnames(m_score) %in% hrd_group)] <- 0
      }
      
      m_score <- m_score[rownames(m_score) %in% group_samples,]
      out <- rankOrderClust(m_score, force.1.value.per.row=T)#$df_ranked
      
      if(!is_definite_biall){
        out$clusters <- rep(1,length(out$clusters))
      }
      
      #return(out$df_ranked)
      return(out)
    })
    
    ## Merge the output of the seperate clustering
    roclust <- mergeRankOrderClustOutput(l_roclust)
    
    #--------- Post-processing ---------#
    ## Add CHORD info
    roclust$df_ranked <- merge(roclust$df_ranked, m_pred, by='row.names', sort=F)
    rownames(roclust$df_ranked) <- roclust$df_ranked$Row.names
    roclust$df_ranked <- roclust$df_ranked[c(colnames(m_pred),colnames(m_score_filt))] ## Select/reorder columns
    
    ## Fix genotype ~ HRD type mismatches for BRCA1/2-like HRD clusters (set mismatched score to 0)
    chord_classes <- colnames(m_pred)
    gene_names <- colnames(m_score_low)
    
    hrd_type <- apply(roclust$df_ranked[,chord_classes],1, function(i){
      chord_classes[as.logical(i)]
    })

    #View(roclust$df_ranked)
    
    hrd_type <- gsub('p_','',hrd_type)
    modif_samples <- c()
    for(sample in names(hrd_type)){
 
      hrd_gene_group <- hrd_groups[[ hrd_type[[sample]] ]]
      
      sel_cols <- !( colnames(roclust$df_ranked) %in% c(hrd_gene_group, chord_classes) )
      row <- unlist(roclust$df_ranked[sample, sel_cols])
      
      if(sum(row)!=0){
        roclust$df_ranked[sample, sel_cols] <- 0
        modif_samples <- c(modif_samples, sample)
      }
    }
    
    ## Modify m_groups after fixing mismatches
    m_groups[modif_samples,'is_not_all_0'] <- 0
    
    ## Reassign cluster number for BRCA1/2-like HRD clusters
    #Get samples with unknown genotype
    unk_geno_samples <- list(
      rownames(m_groups)[ m_groups$p_BRCA2==1 & m_groups$is_not_all_0==0 ],
      rownames(m_groups)[ m_groups$p_BRCA1==1 & m_groups$is_not_all_0==0 ]
    )
    names(unk_geno_samples) <- names(hrd_groups)
    
    #Get corresponding cluster numbers
    unk_clust_nums <- lapply(l_roclust[c('p_BRCA2','p_BRCA1')], function(i){ 
      samples <- i$order 
      unique( roclust$clusters[samples] )
    })
    names(unk_clust_nums) <- names(hrd_groups)
    
    for(i in names(hrd_groups)){
      #i='BRCA2'
      roclust$clusters[ names(roclust$clusters) %in% unk_geno_samples[[i]] ] <- unk_clust_nums[[i]]
    }
    
    ## Reorder
    updated_order <- order(roclust$clusters)
    roclust$df_ranked <- roclust$df_ranked[updated_order,]
    roclust$clusters <- roclust$clusters[updated_order]
    roclust$order <- roclust$order[updated_order]
    
    ## Make sequential cluster numbers
    roclust$clusters <- structure(
      rle2Clusters(rle(roclust$clusters)),
      names=names(roclust$clusters)
    )
    
    return(roclust)
  })()
  
  #--------- Assign genes to clusters ---------#
  assignGeneToCluster <- function(rank.order.clust, min.col.medians=9){
    df <- cbind(rank_order_clust$df_ranked, cluster=rank_order_clust$clusters)
    cluster_names <- unique(df$cluster)
    
    ## Assign genes to clusters
    l <- lapply(cluster_names, function(i){
      #i=1
      df_ss <- df[df$cluster==i,]
      df_ss$cluster <- NULL
      medians <- apply(df_ss,2,median)
      if(max(medians) >= min.col.medians){
        colnames(df_ss)[which.max(medians)]
      } else {
        'unknown'
      }
    })
    
    out <- cluster_names
    names(out) <- unlist(l)
    
    ## Assign HRD type to unknown
    pred_cols <- grep('^p_\\w+',colnames(df),value=T)
    out_ss <- out[names(out)=='unknown']
    df_ss <- unique(df[ df$cluster %in% out_ss, c(pred_cols, 'cluster') ]) ## Use unique to optimize for speed
    
    df_ss$hrd_type <- colnames(df_ss)[ max.col(subset(df_ss,select=-cluster), 'first') ]
    df_ss$hrd_type <- paste0('unknown_',gsub('p_','',df_ss$hrd_type),'-type')
    
    names(out_ss) <- df_ss$hrd_type[ match(out_ss, df_ss$cluster) ]
    names(out)[match(out_ss,out)] <- names(out_ss)
    
    return(out)
  }
  
  rank_order_clust$cluster_names <- assignGeneToCluster(rank_order_clust)
  
  #--------- Sort samples by rank order clustering ---------#
  l_m_diplotypes <- lapply(l_m_diplotypes, function(i){
    i[match(rank_order_clust$order,rownames(i)),]
  })
  
  #========= Select only high impact biallelic hits to get a clean heat map =========#
  message('Performing rank order clustering...')
  mkDiplotypeMatrixFilter <- function(
    df.ranked,
    
    ## Which genes correspond to which HRD type
    hrd_groups=list(
      p_BRCA2=c('BRCA2','RAD51C','PALB2'),
      p_BRCA1=c('BRCA1')
    ),
    
    min.biall.score=NULL, reverse.bool=F
  ){
    #df.ranked=rank_order_clust$df_ranked
    
    df <- df.ranked
    
    ## Convert presence of pathogenic diplotype to 1
    df[df!=0] <- 1
    
    ## Separate chord predictions and diplotype scores
    df_split <- list(
      pred=df[,names(hrd_groups)],
      diplotypes=df[,names(df)!=names(hrd_groups)]
    )
    
    l_grouping <- hrd_groups[apply(df_split$pred,1,which.max)]
    
    counter <- 1
    mask_per_gene <- t(apply(df_split$diplotypes,1,function(v){
      params_sel <- l_grouping[[counter]]
      if(!is.na(params_sel) && all(v==0)){
        v[ names(v) %in% params_sel ] <- 1
        #v[ names(v)!=params_sel ] <- 0
      }
      counter <<- counter+1
      return(v)
    }))
    
    mask_per_diplotype <- mask_per_gene[,rep(1:ncol(mask_per_gene),each=2)]
    colnames(mask_per_diplotype) <- paste0(
      colnames(mask_per_diplotype),
      rep(c('_1','_2'),ncol(mask_per_gene))
    )
    
    if(!is.null(min.biall.score)){
      low_score_samples <- apply(df.ranked[unlist(hrd_groups, use.names=F)], 1, max) < min.biall.score
      mask_per_diplotype[low_score_samples,] <- 0
    }
    
    if(reverse.bool){ mask_per_diplotype <- 1-mask_per_diplotype }
    
    return(mask_per_diplotype)
  }
  
  filter_mask <- mkDiplotypeMatrixFilter(rank_order_clust$df_ranked, reverse.bool=T)
  
  #========= Export diplotype matrices and clustering =========#
  alpha_mask <- mkDiplotypeMatrixFilter(rank_order_clust$df_ranked, min.biall.score=8) ## Convert to integer matrix and invert
  alpha_mask[alpha_mask==0] <- 0.35
  l_m_diplotypes$alpha <- alpha_mask
  
  ## Export 
  saveRDS(l_m_diplotypes, paste0(data_out_dir,'/l_m_diplotypes.rds'))
  saveRDS(rank_order_clust, paste0(data_out_dir,'/rank_order_clust.rds'))
  
  formatLMDiplotypes <- function(l_m_diplotypes){
    #l=l_m_diplotypes_extended
    l <- l_m_diplotypes

    ## Combine chrom and pos matrices
    m_chrom <- as.matrix(l$chrom)
    m_chrom <- apply(m_chrom,2,function(i){ gsub(' ','',i) }) ## Remove whitespace
    m_chrom[m_chrom=='none'] <- '0'


    m_pos <- as.matrix(l$pos)
    m_pos[m_pos==0] <- NA
    m_pos[l$eff=='deep_deletion'] <- NA

    m_chrom_pos <-  matrix(
      paste(m_chrom,m_pos,sep=':'),
      nrow=nrow(m_chrom), ncol=ncol(m_chrom), dimnames=dimnames(m_chrom)
    )

    m_chrom_pos[m_chrom_pos=='0:NA'] <- 'NA'
    l$chrom_pos <- m_chrom_pos
    l$chrom <- NULL
    l$pos <- NULL

    ## Misc post-processing
    l$alpha <- NULL
    l$score <- round(l$score)

    ## Rename some matrices
    names(l)[names(l)=='eff'] <- 'variant'
    names(l)[names(l)=='score'] <- 'p_score'
    names(l)[names(l)=='a_origin'] <- 'variant_origin'
    names(l)[names(l)=='max_score_origin'] <- 'evidence_type'
    
    l <- lapply(l, function(i){
      #i=l[[1]]
      i <- as.data.frame(i)
      sample_names <- rownames(i)
      
      if(convert_cpct_ids){
        sample_names_converted <- metadata$hmf_id[ match(sample_names, metadata$sample) ]
        sample_names_converted[is.na(sample_names_converted)] <- sample_names[is.na(sample_names_converted)]
      }
      
      m <- cbind(index=1:nrow(i), sample=sample_names_converted,cluster=rank_order_clust$clusters,i)
      rownames(m) <- NULL
      return(m)
    })
    
    return(l)
  }
  
  l_m_diplotypes_export <- formatLMDiplotypes(l_m_diplotypes)
  write.xlsx(l_m_diplotypes_export, paste0(data_out_dir,'/l_m_diplotypes.xlsx'))
  
  
  #========= Including more HR genes =========#
  message('Making heat map layers for extended HR genes...')
  l_m_diplotypes_extended <- (function(min.biall.score=9, transparent.alpha=0.35){
    l <- mkDiplotypeMatrixLayers(diplotypes_filt_ss)
    
    unknown_cluster_samples <- names(rank_order_clust$clusters)[ rank_order_clust$clusters %in% c(4,6) ]
    
    m <- l$score
    gene_allele_indexes <- split(
      1:length(colnames(l$score)), 
      rep(1:(length(colnames(m))/2), each = 2)
    )
    
    ## Make biallele hits with high score opaque
    l$alpha <- do.call(cbind,lapply(gene_allele_indexes, function(i){
      m_gene <- m[,i]
      
      m_new_alpha <- do.call(rbind,lapply(rowSums(m_gene), function(j){
        if(j>=min.biall.score){ c(1,1) }
        else { c(transparent.alpha,transparent.alpha) }
      }))
      
      colnames(m_new_alpha) <- colnames(m_gene)
      rownames(m_new_alpha) <- rownames(m_gene)
      
      return(m_new_alpha)
    }))
    
    ## For samples with biallelic hit in known genes, set to transparent
    l$alpha[!(rownames(l$alpha) %in% unknown_cluster_samples),] <- transparent.alpha 
    
    ## Select genes with high biallele score in at least one sample
    sel_gene_allele_indexes <- unlist(lapply(gene_allele_indexes, function(i){
      m_gene <- m[unknown_cluster_samples,i]
      if(any(rowSums(m_gene)>=min.biall.score)){
        i
      } else {
        NULL
      }
    }), use.names=F)
    
    l <- lapply(l,function(i){ i[,sel_gene_allele_indexes] })
    
    # lapply(l, rownames)
    # names(l)
    # names(l_m_diplotypes)
    
    ## Merge with diplotypes matrix with 4 HR genes
    counter <- 0
    lapply(l, function(i){
      counter <<- counter+1
      m_name <- names(l)[counter]
      out <- i[rownames(l_m_diplotypes$eff),]
      out <- cbind(l_m_diplotypes[[m_name]], out)
      return(out)
    })  
  })()
  
  l_m_diplotypes_extended_export <- formatLMDiplotypes(l_m_diplotypes_extended)
  write.xlsx(l_m_diplotypes_extended_export, paste0(data_out_dir,'/l_m_diplotypes_extended.xlsx'))
  
  
  ####################################################################################################
  # Plotting (main)                                                                                  #
  ####################################################################################################
  
  #========= Plot diplotypes =========#
  plotDiplotypeHeatmap <- function(
    l.m.diplotypes=l_m_diplotypes,
    clusters=rank_order_clust$clusters,
    eff.metadata=eff_metadata,
    auto.color='#c2c2c2',
    title=NULL
  ){
    
    # l.m.diplotypes <- l_m_diplotypes2
    # clusters <- sample_clusters
    
    #--------- Prep plot data ---------#
    df_heatmap <- (function(){
      common_colnames <- c('sample','allele')
      
      counter<-0
      l <- lapply(l.m.diplotypes, function(i){ 
        #print(counter)
        #i <- l.m.diplotypes[[4]]
        counter <<- counter+1
        df <- melt(as.matrix(i))
        colnames(df) <- c(common_colnames, names(l.m.diplotypes)[counter])
        return(df)
      })
      
      df <- Reduce(function(x,y){ merge(x, y, by=common_colnames, all.x=T, sort=F) }, l)
      df$sample <- factor(df$sample, unique(df$sample))
      return(df)
    })()
    
    #--------- Add cluster info ---------#
    #if(!is.null(clusters)){ df_heatmap$cluster <- clusters[match(df_heatmap$sample, names(clusters))] }
    
    #--------- Format values ---------#
    ## Order eff from least to most pathogenic
    effs <- unique(eff.metadata$ann_s1[eff.metadata$ann_s1 %in% df_heatmap$eff])
    eff_order <- unique(effs)
    eff_order <- eff_order[eff_order %in% unique(df_heatmap$eff)]
    
    ## Convert eff to human readable eff
    df_heatmap$eff <- eff.metadata$ann_h_s[ match(df_heatmap$eff, eff.metadata$ann_s1) ]
    eff_order <- unique( 
      eff.metadata$ann_h_s[ match(eff_order, eff.metadata$ann_s1) ] 
    )
    
    df_heatmap$eff <- factor(df_heatmap$eff, levels=eff_order)
    levels(df_heatmap$eff)[ levels(df_heatmap$eff) == 'none' ] <- ' '
    
    ## Format booleans
    df_heatmap <- within(df_heatmap,{
      a_origin[a_origin!='som'] <- ' '
      a_origin[a_origin=='som'] <- 'Somatic SNV/indel'
      
      max_score_origin[max_score_origin!='known'] <- ' '
      max_score_origin[max_score_origin=='known'] <- 'Known pathogenic'
    })
    
    df_heatmap$max_score_origin <- factor(df_heatmap$max_score_origin, c('Known pathogenic','Somatic SNV/indel',' '))
    
    #--------- Set colors ---------#
    col_pal <- list()
    
    col_pal$eff <- (function(){
      ## Define colors for highly pathogenic effs; generate the rest
      colors_manual <- (function(){
        v <- structure(eff.metadata$color_h_s, names=eff.metadata$ann_h_s)
        v <- v[ !duplicated(names(v)) ]
        v <- v[v!='auto']
        return(v)
      })()
      
      remain_eff <- eff_order[!(eff_order %in% names(colors_manual))]
      
      if(is.null(auto.color)){
        colors_auto <- (function(){
          pal_name <- 'YlOrBr'
          if(length(remain_eff)==0){ return(NULL) }
          
          n_colors <- length(remain_eff)
          max_pal_colors <- brewer.pal.info[rownames(brewer.pal.info)==pal_name,'maxcolors']
          
          if(n_colors <= max_pal_colors){
            out <- suppressWarnings( brewer.pal(n_colors, pal_name)[1:n_colors] ) ## Prevent min n colors error
          } else {
            out <- colorRampPalette(brewer.pal(max_pal_colors, pal_name))( n_colors )
          }
          
          names(out) <- remain_eff
          
          return(out)
        })()
      } else {
        colors_auto <- structure(rep(auto.color, length(remain_eff)), names=remain_eff)
      }
      
      out <- c(colors_manual, colors_auto)
      out <- out[eff_order]
      return(out)
    })()
    
    ## Meta annotation
    ## Make color vector for scale_color_manual() hack
    col_pal$booleans <- c('Known pathogenic'='red','Somatic SNV/indel'='black',' '=NA)
    
    df_heatmap$a_origin.color <- ifelse(df_heatmap$a_origin=='Somatic SNV/indel',col_pal$booleans['Somatic SNV/indel'],NA)
    df_heatmap$max_score_origin.color <- ifelse(df_heatmap$max_score_origin=='Known pathogenic',col_pal$booleans['Known pathogenic'],NA)
    
    #--------- Heatmap ---------#
    ## Create axis indexes
    df_heatmap$index_y <- as.integer(df_heatmap$allele)
    df_heatmap$index_x <- as.integer(df_heatmap$sample)
    
    ## Use gene name for the corresponding allele columns
    gene_labels <- unique(gsub('_\\d','',levels(df_heatmap$allele)))
    gene_breaks <- seq(1, 2*length(gene_labels), by=2) + 0.5
    
    heatmap <- ggplot(df_heatmap, aes(x=index_x,y=index_y)) +
      
      ## Main tiles (eff)
      geom_tile(aes(fill=eff), alpha=df_heatmap$alpha, color='grey', size=0.1) +
      scale_x_continuous(
        name=sprintf('Patient index (n=%s)',length(unique(df_heatmap$sample))),
        #name='Sample index',
        expand=c(0,0)
      ) +
      
      scale_fill_manual(
        name='', ## Start name with space to ensure legend is first
        values=col_pal$eff,
        breaks=names(col_pal$eff)[ tolower(names(col_pal$eff)) != 'none' ], ## Remove none colors from legend
        guide=guide_legend( 
          label.theme=element_text(size=8), 
          title.theme=element_blank()
          #title.theme=element_text(size=0.1) 
        )
      ) + 
      
      ## Metadata points
      ## Color is same as fill. Then use aes(color) with scale_color_manual() to hack in the legend
      # a_origin
      geom_tile(
        aes(color=a_origin), fill=df_heatmap$a_origin.color, alpha=df_heatmap$alpha,
        height=0.2, position = position_nudge(y=-0.1)
      ) +
      
      # max_score_origin
      geom_tile(
        aes(color=max_score_origin), fill=df_heatmap$max_score_origin.color, alpha=df_heatmap$alpha, 
        height=0.2, position = position_nudge(y=0.1)
      ) +
      
      scale_color_manual(
        name=' ',
        breaks=names(col_pal$booleans),
        values=alpha(col_pal$booleans,0), ## Make alpha 0 so that colors don't bleed into neighboring tiles
        guide=guide_legend(override.aes=list(
          fill=c(col_pal$booleans[1:2],'white'), color=NA),
          label.theme=element_text(size=8),
          title.theme=element_blank()
          #title.theme=element_text(size=0.1)
        )
      ) +
      
      scale_y_continuous(
        breaks=gene_breaks, labels=gene_labels, trans='reverse',
        expand=c(0,0)#, sec.axis=dup_axis()
      ) +
      
      geom_hline(yintercept=seq(1,length(unique(df_heatmap$index_y)),by=2) - 0.5, size=0.25) + ## Draw gene dividers
      
      theme_bw() +
      theme(
        plot.title=element_text(hjust=0.5),
        
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=9),
        #axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        #panel.border = element_rect(fill=NA),
        panel.grid = element_blank(),
        legend.box='vertical',
        legend.spacing.x=unit(0.1,'cm'),
        legend.key.size = unit(0.8,"line"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    
    ## Cluster boundaries
    if(!is.null(clusters)){
      heatmap <- heatmap +  
        geom_vline(xintercept=cumsum(rle(clusters)$lengths) + 0.5, size=0.25)
    }
    
    ## Add side title
    if(!is.null(title)){ heatmap <- heatmap + ylab(title) + theme(axis.title.y.right=element_blank()) }
    else { heatmap <- heatmap + theme(axis.title.y=element_blank()) }
    
    return(heatmap)
  }
  
  #plotDiplotypeHeatmap(l_m_diplotypes, rank_order_clust$clusters)
  
  #========= Plot clusters =========#
  plotClusters <- function(clusters, palette='Set2', hjust=0.2, vjust=0.5, angle=0, ...){
    
    #clusters=sample_clusters
    
    df_clusters <- data.frame(
      sample=factor(names(clusters),unique(names(clusters))),
      cluster=factor(clusters, unique(clusters))
    )
    
    breaks <- c(0, cumsum(rle(clusters)$lengths))
    midpoints <- breaks[-length(breaks)] + diff(breaks)/2
    
    cluster_tiles <- ggplot(df_clusters, aes(y=1, x=sample)) + 
      geom_tile(aes(fill=cluster), show.legend=F) + 
      geom_vline(xintercept=breaks + 0.5, size=0.25) +
      scale_y_continuous(expand=c(0,0)) +
      annotate(
        'text', y=1, x=midpoints, label=levels(df_clusters$cluster), 
        hjust=hjust, vjust=vjust, angle=angle, ...
      ) +
      scale_fill_brewer(palette=palette) +
      theme_void() +
      theme(
        panel.border=element_rect(fill=NA),
        plot.margin = unit(c(0, 1, 2.5, 1), "pt")
      )
    
    return(cluster_tiles)
  }
  
  #plotClusters(rank_order_clust$clusters, angle=90)
  
  #========= Plot sigs and HRD score =========#
  #--------- Prep data ---------#
  sigs_split <- (function(){
    
    sigs <- pred[,colnames(pred) %in% chord_features]
    rownames(sigs) <- pred$sample
    
    hrd_samples <- pred_hrd$sample
    sigs <- sigs[rownames(sigs) %in% hrd_samples, ]
    
    #rank_order_clust$clusters
    
    sigs_split1 <- list(
      snv=sigs[,grep('[A-Z][.][A-Z]', colnames(sigs))],
      indel=sigs[,grep('[a-z][.][a-z]', colnames(sigs))],
      sv=sigs[,grep('[A-Z]_1+', colnames(sigs))]
    )
    
    ## Make selection
    sigs_sel <- list(
      snv=NA,
      indel=NA,
      sv=c('DUP_1e03_1e04_bp','DUP_1e04_1e05_bp')
    )
    
    sigs_split2 <- lapply(names(sigs_split1), function(i){
      
      #i='sv'
      
      sel <- sigs_sel[[i]]
      df <- sigs_split1[[i]]
      
      if(!anyNA(sel)){
        df <- df[,sel]
      }
      
      df <- cbind(sample=rownames(df), df)
      rownames(df) <- NULL
      #df$cluster <- rank_order_clust$clusters[match(df$sample, names(rank_order_clust$clusters))]
      df <- df[match(rownames(l_m_diplotypes$eff),df$sample),]
      return(df)
    })
    
    names(sigs_split2) <- names(sigs_split1)
    
    return(sigs_split2)
  })()
  
  pred_hrd_pd <- (function(){
    df <- pred_hrd[,c('sample','BRCA2','BRCA1')]
    colnames(df)[2:3] <- paste0('P(',colnames(df)[2:3],'-type HRD)')
    
    df <- df[match(rownames(l_m_diplotypes$eff),df$sample),]
    return(df)
  })()
  
  #--------- Main ---------#
  plotHeatmapDefault <- function(
    df, clusters, 
    palette='YlGnBu', trans='identity',title=NULL
  ){
    # df=sigs_split$sv
    # clusters=rank_order_clust$clusters
    
    df$cluster <- clusters[match(df$sample, names(clusters))]
    
    df <- df[match(names(clusters), df$sample),]
    
    df$sample <- factor(df$sample, levels=unique(df$sample))
    df <- melt(df, c('sample','cluster'))
    
    heatmap <- ggplot(df, aes(y=as.integer(variable), x=sample)) +
      geom_tile(aes(fill=value), color='grey', size=0.1) +
      geom_vline(xintercept=which(!duplicated(df$cluster)) - 0.5, size=0.25) + ## Cluster boundaries
      scale_fill_distiller(
        palette=palette, direction=1, trans=trans, limits=c(0, max(df$value)),
        guide=guide_colorbar(
          title=element_blank(),
          direction='horizontal', title.position='top',
          label.position='bottom', label.hjust=0.3, label.vjust=0.5, 
          label.theme=element_text(angle=60, size=8),
          frame.colour='black', ticks.colour='black'
        )
      ) +
      
      scale_y_continuous(
        breaks=unique(as.integer(df$variable)), labels=levels(df$variable),
        expand=c(0,0), trans='reverse'
      ) +
      scale_x_discrete(expand=c(0,0)) +
      
      theme_bw() +
      theme(
        panel.background=element_blank(),
        panel.border=element_rect(fill=NA),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        
        legend.key.height=unit(0.4,'cm'),
        
        plot.title = element_text(hjust=0.5, size=12),
        plot.margin = unit(c(0,0,0,0), 'pt')
      )
    
    ## Add side title
    if(!is.null(title)){ heatmap <- heatmap + ylab(title) + theme(axis.title.y.right=element_blank()) }
    else { heatmap <- heatmap + theme(axis.title.y=element_blank()) }
    
    return(heatmap)
  }
  
  #========= Combine plots =========#
  message('Plotting diplotypes...')
  ## Diplotypes
  diplotype_heatmap <- plotDiplotypeHeatmap(
    l_m_diplotypes, 
    clusters=rank_order_clust$clusters, 
    title='Genotypes'
  )
  
  diplotype_heatmap_extended <- plotDiplotypeHeatmap(
    l_m_diplotypes_extended,
    clusters=rank_order_clust$clusters,
    title='Genotypes'
  )

  ## Cluster tiles
  cluster_tiles <- plotClusters(rank_order_clust$clusters)
  
  ## CHORD features
  counter <- 0
  sig_heatmap <- lapply(sigs_split, function(i){ 
    counter <<- counter+1
    plot_title <- toupper(names(sigs_split)[counter])
    plotHeatmapDefault(i, clusters=rank_order_clust$clusters, title=plot_title) 
  })
  names(sig_heatmap) <- names(sigs_split)
  # sig_heatmap_combined <- plot_grid(plotlist = sig_heatmap, nrow=1, axis='tblr', align='h')
  
  ## CHORD predictions
  chord_heatmap <- plotHeatmapDefault(
    pred_hrd_pd, clusters=rank_order_clust$clusters,
    palette='Reds', title='CHORD'
  )
  
  #========= Final plot =========#
  #--------- Main ---------#
  plots <- list(
    ggplot() + theme_void(), ## Add top padding
    chord_heatmap,
    sig_heatmap$sv,
    cluster_tiles,
    diplotype_heatmap
  )
  rel_heights <- c(
    0.2,
    chord=ncol(pred_hrd_pd)-1,
    sv=ncol(sigs_split$sv)-1,
    cluster=0.7,
    diplotypes=ncol(l_m_diplotypes$eff)+1
  )
  
  plots_combined <- plot_grid(
    plotlist = plots,
    ncol=1, axis='tblr', align='v',rel_heights=rel_heights
  )
  
  png(paste0(plot_out_dir,'/overview_hrd_samples.png'),10, 4.5,'in', res=1000)
  grid.draw(plots_combined)
  dev.off()
  
  #--------- Extended ---------#
  plots_extended <- list(
    ggplot() + theme_void(), ## Add top padding
    chord_heatmap,
    sig_heatmap$sv,
    cluster_tiles,
    diplotype_heatmap_extended
  )
  
  # ## HMF
  # rel_heights_extended <- c(
  #   0.2,
  #   chord=ncol(pred_hrd_pd)-1,
  #   sv=ncol(sigs_split$sv)-1,
  #   cluster=0.7,
  #   diplotypes=ncol(l_m_diplotypes_extended$eff) * 0.5
  # )
  # 
  # plots_combined_extended <- plot_grid(
  #   plotlist = plots_extended,
  #   ncol=1, axis='tblr', align='v',rel_heights=rel_heights_extended
  # )
  # 
  # png(paste0(plot_out_dir,'/overview_hrd_samples_extended.png'),10, 7,'in', res=600)
  # grid.draw(plots_combined_extended)
  # dev.off()
  
  ## PCAWG
  rel_heights_extended <- c(
    0.2,
    chord=ncol(pred_hrd_pd)-1,
    sv=ncol(sigs_split$sv)-1,
    cluster=0.7,
    diplotypes=ncol(l_m_diplotypes_extended$eff)*0.5
  )
  
  plots_combined_extended <- plot_grid(
    plotlist = plots_extended,
    ncol=1, axis='tblr', align='v',rel_heights=rel_heights_extended
  )
  
  png(paste0(plot_out_dir,'/overview_hrd_samples_extended.png'),10, 6.5,'in', res=600)
  grid.draw(plots_combined_extended)
  dev.off()
  
  
  
  ####################################################################################################
  # Plotting (diplotypes by cancer type)                                                             #
  ####################################################################################################
  
  message('Plotting diplotypes per cancer type...')
  
  #metadata <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/metadata/DR-047-update1_metadata.txt'))
  metadata_ss <- metadata[match(names(rank_order_clust$clusters), metadata$sample),c('sample','cancer_type')]
  
  cancer_types_sel <- c('Ovary','Pancreas','Breast','Prostate','Biliary','Urinary tract')
  sample_groups <- lapply(cancer_types_sel, function(i){
    metadata_ss[metadata_ss$cancer_type==i,'sample']
  }); names(sample_groups) <- cancer_types_sel
  
  sample_groups[['Other']] <- metadata_ss[!(metadata_ss$cancer_type %in% cancer_types_sel),'sample']
  
  sample_clusters <- rep(names(sample_groups), sapply(sample_groups, length, USE.NAMES=F))
  names(sample_clusters) <- unlist(sample_groups, use.names=F)
  
  l_m_diplotypes2 <- (function(){
    subsetDiplotypeSamples <- function(l.m.diplotypes, samples){
      lapply(l.m.diplotypes, function(i){
        i[match(samples,rownames(i)),]
      })
    }
    
    l_m_ss <- lapply(sample_groups, function(i){
      #i=sample_groups[[1]]
      subsetDiplotypeSamples(l_m_diplotypes,i)
    })
    
    datatype_names <- names(l_m_diplotypes)
    
    out <- lapply(datatype_names, function(i){
      m <- do.call(rbind, lapply(l_m_ss, `[[`, i))
      rownames(m) <- gsub('^.+[.]','',rownames(m))
      return(as.data.frame(m))
    })
    
    names(out) <- datatype_names
    
    return(out)
  })()
  
  chord_heatmap2 <- plotHeatmapDefault(pred_hrd_pd, clusters=sample_clusters,palette='Reds', title='CHORD')
  sig_heatmap2 <- plotHeatmapDefault(sigs_split$sv, clusters=sample_clusters, title='SV')
  
  cluster_tiles2 <- plotClusters(sample_clusters, palette='Pastel1', hjust=0.5, angle=90, size=2.5)
  
  diplotype_heatmap2 <- plotDiplotypeHeatmap(
    l_m_diplotypes2, title='Genotypes',
    clusters=structure(rle2Clusters(rle(sample_clusters)), names=names(sample_clusters))
  )
  
  plots2 <- list(
    ggplot() + theme_void(), ## Add top padding
    chord_heatmap2,
    sig_heatmap2,
    cluster_tiles2,
    diplotype_heatmap2
  )
  rel_heights2 <- c(
    0.2,
    chord=ncol(pred_hrd_pd)-1,
    sv=ncol(sigs_split$sv)-1,
    cluster=2,
    diplotypes=ncol(l_m_diplotypes2$eff)+1
  )
  
  plots_combined2 <- plot_grid(
    plotlist=plots2, ncol=1, axis='tblr', align='v',
    rel_heights=rel_heights2
  )
  
  png(paste0(plot_out_dir,'/overview_hrd_samples_by_cancer_type.png'),10, 5,'in', res=1000)
  grid.draw(plots_combined2)
  dev.off()
}







