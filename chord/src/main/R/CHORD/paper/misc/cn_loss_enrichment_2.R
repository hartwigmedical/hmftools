library(ggplot2)
library(ggrepel)
library(cowplot)
library(openxlsx)
# library(reshape2)
# library(cowplot)#; theme_set(theme_grey())
# library(grid)
library(RColorBrewer)

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

#========= Load data =========#
#--------- CHORD predictions ---------#
metadata <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/metadata_for_analysis.txt'))
pancancer_analysis_samples <- metadata$sample[ metadata$used_for_analysis ]

pred <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/pred_analysis_samples.txt'))
hrd_samples <- pred[pred$hrd >= 0.5,'sample']
hrp_samples <- pred[pred$hrd < 0.5,'sample']

genes_bed <- read.delim(paste0(base_dir,'/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/genes_chord_paper.bed'))
colnames(genes_bed)[1] <- 'chrom'

#--------- Monoallelic loss samples ---------#
rank_order_clust <- readRDS(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/rank_order_clust.rds'))
#biallel_loss_hrd_samples <- names(rank_order_clust$clusters)[!(rank_order_clust$clusters %in% c(4,6))]
# deep_del_samples <- rownames(rank_order_clust$df_ranked)[
#    apply(rank_order_clust$df_ranked,1,function(i){ any(i==10.6) })
# ]
   

#--------- Diplotypes ---------#
diplotypes <- (function(){
   path <- '/Users/lnguyen/Documents/R_cache/gene_diplotypes_with_amps_HMF_DR010_DR047.txt.gz'
   
   if(!file.exists(path)){
      file.copy(
         paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/gene_diplotypes_with_amps_HMF_PCAWG.txt.gz'),
         path
      )
   }
   
   read.delim(path)
})()

#========= Functions =========#
## Do fisher test
fisherTestCustom <- function(df, by.col='has_cn_loss'){
   
   ## Prep data
   df_ss <- subset(
      df,
      ensembl_gene_id %in% genes_bed$ensembl_gene_id,
      select=c(sample, hgnc_symbol, ensembl_gene_id, a1.eff)
   )
   df_ss$is_hrd <- df_ss$sample %in% hrd_samples
   
   df_ss$has_cn_loss <- df_ss$a1.eff %in% c('deep_deletion','loh')
   df_ss$has_loh <- df_ss$a1.eff == 'loh'
   
   uniq_samples <- unique(df_ss$sample)
   
   n_hrd <- sum(uniq_samples %in% hrd_samples)
   n_hrp <- sum(uniq_samples %in% hrp_samples)
   
   main <- function(df){
      #df=subset(df_ss, ensembl_gene_id=='ENSG00000012048')
      #df[order(df$sample),]
      
      conting <- rbind(
         c(
            sum(df$is_hrd & df[[by.col]]), 
            sum(!df$is_hrd & df[[by.col]])
         ),
         
         c(
            n_hrd, 
            n_hrp
         )
      )
      
      data.frame(
         n_pos=conting[1,1],
         n_hrd,
         n_neg=conting[1,2],
         n_hrp,
         pvalue=fisher.test(conting, alternative='greater')$p.value
      )
   }
   
   l <- split(df_ss[,c('is_hrd',by.col)], df_ss$ensembl_gene_id)
   out <- do.call(rbind, lapply(l,main) )
   out$qvalue <- p.adjust(out$pvalue,'hochberg')
   
   out <- cbind(ensembl_gene_id=rownames(out),out)
   
   out <- cbind(
      chrom=genes_bed[match(out$ensembl_gene_id, genes_bed$ensembl_gene_id),'chrom'],
      start=genes_bed[match(out$ensembl_gene_id, genes_bed$ensembl_gene_id),'start'],
      end=genes_bed[match(out$ensembl_gene_id, genes_bed$ensembl_gene_id),'end'],
      hgnc_symbol=genes_bed[match(out$ensembl_gene_id, genes_bed$ensembl_gene_id),'hgnc_symbol'],
      out
   )
   rownames(out) <- NULL
   
   out <- out[order(out$pvalue),]
   out$rank <- 1:nrow(out)
   
   #out[out$hgnc_symbol %in% c('BRCA2','BRCA1','RAD51C','PALB2'),]
   return(out)
}

## Plot
plotEnrichment <- function(df, title=NULL, sel.genes=c('BRCA1','BRCA2','RAD51C','PALB2')){
   
   #df=cn_loss_enrichment$pc_analysis_samples_cluster_4_6
   ## Add plotting data
   
   df$label <- (function(){
      v <- paste0(
         df$hgnc_symbol,'; rank = ',df$rank,
         '\nq = ', signif(df$qvalue,2),
         '\nHRD: ',df$n_pos,'/',df$n_hrd,', HRP: ',df$n_neg,'/',df$n_hrp
      )
      v[!(df$hgnc_symbol %in% sel.genes)] <- NA
      return(v)
   })()
   
   df$label_short <- (function(){
      v <- df$hgnc_symbol
      v[!(df$hgnc_symbol %in% sel.genes)] <- NA
      return(v)
   })()
   
   df$is_selected_gene <- !is.na(df$label)
   
   df$chrom <- factor(df$chrom,c(1:22,'X'))
   #brewer.pal.info
   #display.brewer.all()
   chrom_colors <- c(
      brewer.pal(8,'Dark2'),
      brewer.pal(9,'Set1'),
      brewer.pal(12,'Set3')
   )
   
   p <- ggplot(df, aes(x=rank, y=-log10(qvalue), fill=chrom, size=is_selected_gene)) + 
      
      geom_bar(stat='identity', width=1) +
      geom_point(shape=21) +
      geom_point(data=subset(df, is_selected_gene), shape=21) +
      scale_fill_manual(values=chrom_colors) +
      scale_size_manual(values=c('TRUE'=3.5,'FALSE'=1)) +
      guides(size=FALSE) +
      
      geom_label_repel(aes(label=label),nudge_x=200, nudge_y=2, fill='white', size=2.5) +
      xlab(paste0(nrow(df),' genes')) +
      
      theme_bw() +
      theme(
         plot.title=element_text(hjust=0.5),
         plot.subtitle=element_text(hjust=0.5)
      )
   
   if(!is.null(title)){ p <- p + ggtitle(title) }
   
   return(p)
}


#========= Exec =========#
# cn_loss_enrichment <- (function(){
#    df <- subset(
#       diplotypes, 
#       sample %in% pancancer_analysis_samples 
#       & !(sample %in% deep_del_samples)
#    )
#    fisherTestCustom(df, by.col='has_loh')
# })()
# write.xlsx(
#    cn_loss_enrichment, 
#    paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/data/cn_loss_enrichment.xlsx')
# )


p_cn_loss_enrichment <- plotEnrichment(cn_loss_enrichment)
pdf(paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/plots/cn_loss_enrichment.pdf'),8,5)
plot(p_cn_loss_enrichment)
dev.off()


# plotEnrichmentPerCancerType <- function(df, by.col){
#    #df=subset(diplotypes, sample %in% pancancer_analysis_samples & !(sample %in% biallel_loss_hrd_samples))
#    
#    df$cancer_type <- metadata[match(df$sample, metadata$sample),'cancer_type']
#    
#    sel_cancer_types <- c('Ovary','Breast','Prostate','Pancreas')
#    df$cancer_type <- factor(df$cancer_type, sel_cancer_types)
#    
#    l_pd <- lapply(split(df, df$cancer_type), function(i){
#       fisherTestCustom(i, by.col=by.col)
#    })
#    
#    l_p <- lapply(names(l_pd), function(i){
#       plotEnrichment(l_pd[[i]], title=i)
#    })
#    
#    return(l_p)
# }
# 
# p_cn_loss_enrichment_per_ct <- cn_loss_enrichment <- (function(){
#    df <- subset(
#       diplotypes, 
#       sample %in% pancancer_analysis_samples 
#       & !(sample %in% deep_del_samples)
#    )
#    plotEnrichmentPerCancerType(df, by.col='has_loh')
# })()














p_cn_loss_enrichment <- lapply(names(cn_loss_enrichment), function(i){
   plotEnrichment(cn_loss_enrichment[[i]], title=i)
})

pdf(paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/plots/cn_loss_enrichment.pdf'),6,5)
for(i in p_cn_loss_enrichment){ plot(i) }
dev.off()

write.xlsx(
   cn_loss_enrichment, 
   paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/data/cn_loss_enrichment.xlsx')
)


#========= Fisher test for main HRD cancer types =========#
plotEnrichmentPerCancerType <- function(df){
   #df=subset(diplotypes, sample %in% pancancer_analysis_samples & !(sample %in% biallel_loss_hrd_samples))
   
   df$cancer_type <- metadata[match(df$sample, metadata$sample_id),'primary_tumor_location']
   
   sel_cancer_types <- c('Ovary','Breast','Prostate','Pancreas','Biliary','Urinary tract')
   df$cancer_type <- factor(df$cancer_type, sel_cancer_types)
   
   cn_loss_enrichment_per_ct <- lapply(split(df, df$cancer_type), function(i){
      fisherTestCustom(i, by.col='has_loh')
   })
   
   p_cn_loss_enrichment_per_ct <- lapply(names(cn_loss_enrichment_per_ct), function(i){
      plotEnrichment(cn_loss_enrichment_per_ct[[i]], title=i)
   })
   
   return(p_cn_loss_enrichment_per_ct)
}

p_cn_loss_enrichment_per_ct <- plotEnrichmentPerCancerType(
   subset(diplotypes, sample %in% pancancer_analysis_samples & !(sample %in% biallel_loss_hrd_samples))
)


pdf(paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/plots/cn_loss_enrichment_per_ct_2.pdf'),18,10)
plot_grid(
   plotlist=p_cn_loss_enrichment_per_ct,
   ncol=3, axis='tblr', align='hv'
)
dev.off()


#========= Plot freq of gene defs across genome =========#
cn_loss_per_bin <- readRDS(paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/data/cn_loss_per_bin_hrd_samples.rds'))
cn_loss_per_bin <- cn_loss_per_bin[,colnames(cn_loss_per_bin) %in% pancancer_analysis_samples]

plotFreqCnLoss <- function(fisher_out, cn_loss_per_bin, plot_title=NULL){
   
   #--------- Prep CN loss profiles ---------#
   cn_loss_per_bin_sum <- as.data.frame(do.call(rbind, strsplit(rownames(cn_loss_per_bin),':')))
   colnames(cn_loss_per_bin_sum) <- c('chrom','start')
   cn_loss_per_bin_sum$start <- as.integer(cn_loss_per_bin_sum$start)
   cn_loss_per_bin_sum$chrom <- factor(cn_loss_per_bin_sum$chrom, unique(cn_loss_per_bin_sum$chrom))
   
   cn_loss_per_bin_sum$counts <- rowSums(cn_loss_per_bin)
   
   cn_loss_per_bin_sum$chrom <- factor(cn_loss_per_bin_sum$chrom, c(1:22,'X'))
   
   total_samples <- ncol(cn_loss_per_bin)
   
   #--------- Prep fisher test output ---------#
   #fisher_out <- cn_loss_enrichment[[2]]
   fisher_out$chrom <- factor(fisher_out$chrom, c(1:22,'X'))
   
   #--------- Plot ---------#
   p <- ggplot(cn_loss_per_bin_sum, aes(x=start,y=counts)) +
      facet_grid(.~chrom, scales='free_x') +
      geom_line() +
      
      geom_point(data=fisher_out, aes(fill=chrom, x=start, y=n_pos), shape=21) +
      guides(fill=F) + 
      geom_label_repel(
         data=fisher_out, aes(x=start, y=n_pos, label=label_short), 
         nudge_y=25, fill='white', color='red', size=3
      ) +
      
      xlab('Genome position') +
      ylab(sprintf('n HRD samples with CN loss\n(total no. samples = %s)', total_samples)) +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank()
      )
   
   if(!is.null(plot_title)){
      p <- p + ggtitle(plot_title)
   }
   
   return(p)
}

#--------- Pancancer ---------#
p_freq_cn_loss <- plotFreqCnLoss(
   cn_loss_enrichment$pancancer_analysis_samples, 
   cn_loss_per_bin
)

pdf(paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/plots/n_hrd_cn_loss.pdf'),20,5)
plot(p_freq_cn_loss)
dev.off()

#--------- Per cancer type ---------#
cancer_type_samples <- split(metadata$sample_id, metadata$primary_tumor_location)
cancer_type_samples <- cancer_type_samples[sel_cancer_types]

cn_loss_per_bin_per_ct <- lapply(cancer_type_samples, function(i){
   cn_loss_per_bin[,colnames(cn_loss_per_bin) %in% i]
})

p_cn_loss_per_bin_per_ct <- lapply(sel_cancer_types, function(i){
   #i='Ovary'
   df1 <- cn_loss_enrichment_per_ct[[i]]
   df2 <- cn_loss_per_bin[,colnames(cn_loss_per_bin) %in% cancer_type_samples[[i]] ]
   
   plotFreqCnLoss(df1, df2, plot_title=i)
})

pdf(paste0(base_dir,'/CHORDv2/analysis/cn_loss_enrichment/plots/n_hrd_cn_loss_per_ct.pdf'),20,5)
for(i in p_cn_loss_per_bin_per_ct){ plot(i) }
dev.off()








