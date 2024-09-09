library(reshape2)
library(grid)
library(ggplot2)
library(ggrepel)
library(cowplot)#; theme_set(theme_grey())

#library(dndscv)

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

#========= Load data =========#
#--------- Other data ---------#
genes_bed <- read.delim(paste0(base_dir,'/CHORD/scripts_main/hmfGeneAnnotation/inst/misc/genes_chord_paper.bed'))
colnames(genes_bed)[1] <- 'chrom'

#--------- CHORD predictions ---------#
#metadata_hmf <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))

l_metadata <- list()
l_metadata$hmf <- read.delim(paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'))
l_metadata$pcawg <-  read.delim(paste0(base_dir,'/datasets/PCAWG_2020/metadata/pcawg_metadata_ann.txt.gz'))

## Select samples for analysis here
pred <- (function(){
   l <- readRDS(paste0(base_dir,'/CHORDv2/training/models/2.02_mergedDelMh2-5/plots/preds.rds'))
   l <- l[c('HMF','PCAWG')]
   l$PCAWG$response <- l$PCAWG$response_umcu
   l$PCAWG$response_umcu <- NULL
   
   common_cols <- Reduce(intersect, lapply(l,colnames))
   
   l <- lapply(l, `[`, common_cols)
   
   counter <- 0
   l <- lapply(l, function(i){
      counter <<- counter+1
      i$cohort <- names(l)[counter]
      return(i)
   })
   
   df <- do.call(rbind, l)
   rownames(df) <- NULL
   
   ## Remove HMF tumors from same patient
   df <- df[
      !(df$sample %in% subset(l_metadata$hmf, !max_purity_biopsy, sample, drop=T))
   ,]
   
   ## Keep samples that pass CHORD QC
   ## Keep PCAWG samples with genetic annotation
   df <- df[
      df$hr_status != 'cannot_be_determined' & df$hrd_type != 'cannot_be_determined'
      & !is.na(df$response)
   ,]

   return(df)
})()

hrd_samples <- pred[pred$hrd>=0.5,'sample']

write.table(
   pred, paste0(script_dir,'/data/pred_analysis_samples.txt'),
   sep='\t',quote=F,row.names=F
)

## Save formatted metadata for downstream analysis
(function(){
   df1 <- l_metadata$hmf[,c('sample','hmf_id','primary_tumor_location')]
   colnames(df1)[colnames(df1)=='primary_tumor_location'] <- 'cancer_type'
   df1$cohort <- 'HMF'

   df2 <- l_metadata$pcawg[,c('sample','cancer_type')]
   df2$hmf_id <- NA
   df2$cohort <- 'PCAWG'
   
   df <- rbind(df1,df2)
   
   df$used_for_analysis <- df$sample %in% pred$sample
   
   write.table(
      df, paste0(script_dir,'/data/metadata_for_analysis.txt'),
      sep='\t',quote=F,row.names=F
   )
})()


#--------- Read and cache diplotype table ---------#
diplotypes_raw <- (function(){
   path <- '/Users/lnguyen/Documents/R_cache/gene_diplotypes_max_hmf_pcawg.txt.gz'
   if(!file.exists(path)){
      file.copy(
         paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/gene_diplotypes_with_amps_HMF_PCAWG.txt.gz'),
         path
      )
   }
   read.delim(path)
})()

diplotypes <- diplotypes_raw[diplotypes_raw$sample %in% pred$sample,]


####################################################################################################
# Filtering high frequency variants                                                                #
####################################################################################################

#========= Set AF cutoff and make variant PON ==========#
calcUniqVariantFreq <- function(diplotypes, mode='germ'){
   
   sel_cols <- list(
      common=c('hgnc_symbol'),
      #allele,
      per_allele=c('hgvs_c','max_score','max_score_origin')
   )
   
   getColSelection <- function(allele='a2'){
      sel_cols$per_allele <- paste0(allele,'.',sel_cols$per_allele)
      return(unlist(sel_cols, use.names=F))
   }
   
   ## cnv_germ/cnv_som should contain all possible germ/som variants. No need to scan germ_som rows
   if(mode=='germ'){
      col_names <- getColSelection('a2')
      diplotypes_ss <- diplotypes[diplotypes$diplotype_origin=='cnv_germ',col_names]
   } else if(mode=='som'){
      col_names <- getColSelection('a2')
      diplotypes_ss <- diplotypes[diplotypes$diplotype_origin=='cnv_som',col_names]
   }
   #diplotypes_ss[is.na(diplotypes_ss)] <- 'none'
   #head(diplotypes_ss)
   
   message("Converting variant hgnc_symbol and hgvs_c to string...")
   #diplotypes_ss_string <- with(diplotypes_ss,{ paste0(hgnc_symbol,':',a2.hgvs_c) })
   diplotypes_ss_string <- do.call(function(...){ paste(..., sep=':') }, diplotypes_ss)
   
   message("Counting unique occurrences...")
   v <- table(diplotypes_ss_string)
   #names(v)[grepl('FANCF:c.*1338dupA',names(v))]
   
   ## Reformat table as data frame
   df <- as.data.frame(do.call(rbind,strsplit(names(v),':')))
   colnames(df) <- unlist(sel_cols)
   
   df$freq <- as.integer(v)
   df$freq_prop <- df$freq/length(unique(diplotypes$sample))
   
   ## Remove row where hgvs_c is 0. These used to be NA values
   df <- df[df$hgvs_c!='none',]
   
   ## Assign rank
   df <- df[order(df$freq,decreasing=T),]
   df <- cbind(rank=1:nrow(df),df)
   
   ## Annotations
   df$is_known_patho <- df$max_score==5 & grepl('clinvar|enigma',df$max_score_origin)
   
   # ## Mark pon
   # max_freq <- (function(){
   #    df_known_patho <- df[df$is_known_patho,]
   #    round(median(df_known_patho[n.top.known.patho,'freq']))
   # })()
   # df$is_in_pon <- df$freq > max_freq
   # 
   # ## Store pon cutoff as class
   # class(df) <- c(
   #    pon_cutoff=max_freq, 
   #    origin=if(mode=='germ'){'germline'}else if(mode=='som'){'somatic'}, 
   #    data_type='data.frame'
   # )
   
   return(df)
}

## The cnv_germ subset includes all germline variants found in the germ_som subset. Same for cnv_som.
pon_germ <- calcUniqVariantFreq(diplotypes,'germ')
max_freq <- with(pon_germ,{
   max(freq[is_known_patho & freq <= 100])
})
pon_germ$is_in_pon <- pon_germ$freq > max_freq

# write.table(
#    pon_germ, paste0(script_dir,'/data/pon_germ.txt'),
#    sep='\t', quote=F, row.names=F
# )

plotPon <- function(df, pon.cutoff, n.samples=length(unique(diplotypes$sample))){
   # df=pon_germ
   # pon.cutoff=max_freq
   
   pon_cutoff <- pon.cutoff/n.samples
   
   df$freq_prop <- df$freq / n.samples
   df$is_known_patho <- factor(df$is_known_patho, levels=c('TRUE','FALSE'))
   col_pal <- c('TRUE'='red','FALSE'='lightgrey')
   
   ggplot() +
      geom_point(df, mapping=aes(rank, freq_prop, color=is_known_patho)) + 
      geom_point(subset(df,is_known_patho=='TRUE'), mapping=aes(rank, freq_prop, color=is_known_patho)) +
      geom_hline(yintercept=pon_cutoff, linetype='dashed') +
      
      scale_y_log10(name='Frequency', labels=function(x){ paste0(x*100, '%') }) +
      scale_color_manual(values=col_pal, breaks='TRUE', labels='Known pathogenic \nvariant in ClinVar') +
      labs(title='Frequency of unique SNVs/indels',x='Variant rank') +

      annotate(
         'text', y=pon_cutoff, x=max(df$rank)/2, vjust=-1,
         label=paste0('High frequency: >', format(round(100*pon_cutoff,2), nsmall=2), '%')
      ) +
      
      theme_bw() +
      theme(
         plot.title=element_text(hjust=0.5, size=11),
         legend.position=c(0.8,0.8),
         legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.background=element_rect(color='black', size=0.25)
      )
}
pon_plot <- plotPon(pon_germ, pon.cutoff=max_freq)

## Export
# pdf(paste0(script_dir,'/plots/pon_germ.pdf'), 5,5)
# grid.draw(pon_plot)
# dev.off()

#========= Apply PON filter ==========#
applyPonFilter <- function(diplotypes, pon.germ=NULL, pon.som=NULL){
   # head(diplotypes)
   # pon.germ=pon_germ
   # pon.som=NULL
   
   mkBlacklist <- function(df){
      df_ss <- df[df$is_in_pon, c('hgnc_symbol','hgvs_c')]
      apply(df_ss,1,function(i){ paste(i,collapse=':') })
   }
   
   blacklists <- list()
   if(!is.null(pon.germ)){
      blacklists$germ <- mkBlacklist(pon.germ)
   }
   
   if(!is.null(pon.som)){
      blacklists$som <- mkBlacklist(pon.som)
   }
   
   message("Splitting diplotypes table by diplotype_origin...")
   diplotype_origins <- c('cnv_cnv','cnv_germ','cnv_som','som_som','germ_som')
   diplotypes_split <- lapply(diplotype_origins, function(i){
      #i <- 'cnv_germ'
      df_ss <- diplotypes[diplotypes$diplotype_origin==i, ]
      
      if(nrow(df_ss)==0){ return(NULL) }
      
      ## Initiate in_pon columns here so that rbind doesn't break later
      df_ss$a1.in_pon <- F
      df_ss$a2.in_pon <- F
      
      ## Store old scores
      df_ss$a1.max_score_old <- df_ss$a1.max_score
      df_ss$a2.max_score_old <- df_ss$a2.max_score
      
      # df_ss$a1.hgvs_c[is.na(df_ss$a1.hgvs_c),] <- 'none'
      # df_ss$a2.hgvs_c[is.na(df_ss$a2.hgvs_c),] <- 'none'
      
      return(df_ss)
   }); names(diplotypes_split) <- diplotype_origins
   
   adjustMaxScores <- function(df, allele, blacklist){
      # df=diplotypes_split$cnv_germ
      # allele='a2'
      # blacklist=blacklists$germ
      
      sel_cols <- c('hgnc_symbol', paste0(allele,'.hgvs_c'))
      #variant_strings <- apply(df[sel_cols],1,function(j){ paste(j,collapse=':') })
      variant_strings <- do.call(function(...){ paste(..., sep=':') }, df[sel_cols])
      
      ## Set scores to zero for common variants
      in_pon <- variant_strings %in% blacklist
      df[in_pon,paste0(allele,'.max_score')] <- 0
      df[in_pon,paste0(allele,'.in_pon')] <- T
      
      return(df)
   }
   
   if(!is.null(pon.germ)){
      message("Applying PON filter to germline variants...")
      diplotypes_split$cnv_germ <- adjustMaxScores(diplotypes_split$cnv_germ,'a2',blacklists$germ)
      diplotypes_split$germ_som <- adjustMaxScores(diplotypes_split$germ_som,'a1',blacklists$germ)
   }

   if(!is.null(pon.som)){
      message("Applying PON filter to somatic variants...")
      diplotypes_split$cnv_som <- adjustMaxScores(diplotypes_split$cnv_som,'a2',blacklists$som)
      diplotypes_split$germ_som <- adjustMaxScores(diplotypes_split$germ_som,'a2',blacklists$som)
   }
   
   message("Merging back split diplotypes table...")
   do.call(rbind, diplotypes_split)
}

detIsDef <- function(
   diplotypes,
   min.allele.scores=list(
      cnv_cnv=c(5,5), 
      cnv_som=c(5,3), 
      cnv_germ=c(5,4), 
      som_som=c(5,5), 
      germ_som=c(5,5)
   )
){
   
   counter <- 0
   pb <- txtProgressBar(max=nrow(diplotypes), style=3)
   
   diplotypes$is_def <- with(diplotypes, {
      unlist(Map(function(diplotype_origin, a1.max_score, a2.max_score){
         
         counter <<- counter + 1
         setTxtProgressBar(pb, counter)
         
         cutoffs <- min.allele.scores[[diplotype_origin]]
         a1.max_score >= cutoffs[1] & a2.max_score >= cutoffs[2]
      }, diplotype_origin, a1.max_score, a2.max_score, USE.NAMES=F))
   })
   
   return(diplotypes)
}

diplotypes_filt <- applyPonFilter(diplotypes, pon_germ)
diplotypes_filt <- detIsDef(diplotypes_filt)

path_diplotypes_filt <- paste0(script_dir,'/data/diplotypes_filt.txt.gz')
if(!file.exists(path_diplotypes_filt)){
   write.table(diplotypes_filt, gzfile(path_diplotypes_filt), sep='\t', quote=F, row.names=F)
}




####################################################################################################
# Determining important HR genes                                                                   #
####################################################################################################

#========= Fisher test =========#
fisherTestCustom <- function(
   diplotypes, pred, genes=NULL, 
   chord.cutoff=0.5, 
   fisher.cutoff=0.05, use.pvalue=F,
   min.freq.p=5, 
   rm.likely.non.hr.genes=F
){
   if(is.null(genes)){
      genes <- unique(diplotypes$ensembl_gene_id)
   }
   
   message("Preparing input table...")
   #diplotypes_ss <- diplotypes[,c('sample','ensembl_gene_id','diplotype_origin','a1.max_score','a2.max_score')]
   diplotypes_ss <- diplotypes
   diplotypes_ss$hrd <- pred[match(diplotypes_ss$sample, pred$sample),'hrd']
   
   diplotypes_ss$is_hrd <- diplotypes_ss$hrd >= chord.cutoff
   
   n_samples <- length(unique(diplotypes_ss$sample))
   
   n_hrd_tmp <- unique(diplotypes_ss[c('sample','hrd')])
   n_hrd <- sum(n_hrd_tmp$hrd >= chord.cutoff)
   rm(n_hrd_tmp)
   
   n_hrp <- n_samples-n_hrd
   
   fisherTestByGene <- function(ensembl_gene_id){
      #hgnc_symbol='BRCA2'
      
      df <- diplotypes_ss[diplotypes_ss$ensembl_gene_id==ensembl_gene_id,]
      
      p <- sum(df$is_def & df$is_hrd)
      n <- sum(df$is_def & !df$is_hrd)
      
      m <- rbind(
         c(p, n),
         c(n_hrd, n_hrp)
      )
      
      fisher <- fisher.test(m, alternative='greater')$p.value
      
      out <- data.frame(
         ensembl_gene_id,
         p, p_total=n_hrd, p_rel=p/n_hrd, 
         n, n_total=n_hrp, n_rel=n/n_hrp,
         pvalue=fisher
      )
      out$p_rel[is.na(out$p_rel)] <- 0
      return(out)
   }
   
   ## Main
   counter <- 0
   n_genes <- length(genes)

   message("Performing Fisher's exact test...")
   pb <- txtProgressBar(max=n_genes, style=3)
   freqs <- do.call(rbind, lapply(genes, function(i){
      counter <<- counter + 1
      #i='BRCA2'
      setTxtProgressBar(pb, counter)
      fisherTestByGene(i)
   }))
   freqs <- freqs[order(freqs$pvalue),]
   freqs$qvalue <- p.adjust(freqs$pvalue,'hochberg')

   if(!use.pvalue){
      freqs$keep <- as.integer(freqs$qvalue < fisher.cutoff & freqs$p >= min.freq.p )
   } else {
      freqs$keep <- as.integer(freqs$pvalue < fisher.cutoff & freqs$p >= min.freq.p)
   }
   
   
   if(rm.likely.non.hr.genes){ 
      message(sprintf('Removing samples under with q-values under %s...',fisher.cutoff))
      freqs <- freqs[freqs$keep==T,] 
   }
   
   ## Attach metadata
   class(freqs) <- c(
      cutoff=fisher.cutoff,
      # n_hrd=n_hrd,
      # n_hrp=n_hrp,
      # n_samples=n_samples,
      data_type='data.frame'
   )
   
   return(freqs)
}

path_hr_gene_likelihood <- paste0(script_dir,'/data/hr_gene_likelihood.txt')

if(file.exists(path_hr_gene_likelihood)){
   hr_gene_likelihood <- read.delim(path_hr_gene_likelihood)
} else {
   hr_gene_likelihood <- fisherTestCustom(diplotypes_filt, pred)
   hr_gene_likelihood <- merge(genes_bed, hr_gene_likelihood, by='ensembl_gene_id')
   hr_gene_likelihood <- hr_gene_likelihood[order(hr_gene_likelihood$pvalue),]
   
   #write.table(hr_gene_likelihood, path_hr_gene_likelihood, sep='\t',row.names=F,quote=F)
}

#========= Plot =========#
plotFisherTest <- function(
   df, fisher.cutoff=0.05, min.freq.p=5,
   signif.cutoffs=10^-c(0:5), legend.sci.nota.cutoff=0.01,
   plot.title=NULL, x.lab=NULL, y.lab=NULL, rel.freq.on.xy=F, legend.in.plot=F,
   xlims=c(NA,NA), ylims=c(NA,NA),
   order.rows.by.gene.name=T
){
   #df=hr_gene_likelihood
   # df=cn_eff_coocc[[1]]
   # signif.cutoffs=10^-c(0,100,200,300)
   # fisher.cutoff=10^-100
   
   ## Transform
   df <- within(df,{
      pvalue_log <- -log10(pvalue)
      qvalue_log <- -log10(qvalue)
   })
   
   ## Labels
   df <- within(df,{
      label <- hgnc_symbol
      label[!(qvalue < fisher.cutoff)] <- ''
      label[!(p >= min.freq.p)] <- ''
   })
   
   ## Signif bins
   signif.cutoffs <- sort(signif.cutoffs, decreasing=T)
   cutoffs <- -log10(signif.cutoffs)
   
   df$signif_bin <- .bincode(df$qvalue_log, breaks=c(cutoffs, Inf))
   df$signif_bin[is.na(df$signif_bin)] <- 1
   
   #signif_bin_labels <- signif.cutoffs[df$signif_bin]
   signif_bin_labels <- sapply(format(signif.cutoffs,scientific=T), function(i){
      #i="1e-02"
      i_numeric <- as.numeric(i)
      if(i_numeric >= legend.sci.nota.cutoff){
         ## Regular notation
         out <- format(i_numeric, scientific=F, drop0trailing=T)
      } else {
         ## Exponent notation
         out <- sub('1e[+-]','10^-',i)
         out <- sub('[-]0+','-',out)
      }
      return(out)
   }, USE.NAMES=F)
   signif_bin_labels <- paste0('q<',signif_bin_labels)
   signif_bin_labels <- parse(text=signif_bin_labels)
   
   ##
   if(order.rows.by.gene.name){
      df <- df[order(df$hgnc_symbol),]
   }
   
   df$index <- 1:nrow(df)

   superscriptNotation <- function(x){
      #x=c(0, 10, 20)
      out <- paste0('10^-',x)
      out <- gsub('^10\\^-0$', 1, out) ## Change log10(0) to 1
      parse(text=out)
   }

   ## Main
   if(!rel.freq.on.xy){
      plot <- ggplot(df, aes(index, qvalue_log, label=label)) +
         geom_hline(yintercept=-log10(fisher.cutoff),linetype='dotted') +
         geom_point(aes(size=signif_bin, fill=signif_bin), shape=21, stroke=0.3) +
         geom_text_repel(
            point.padding=0.25, nudge_y=0.6, min.segment.length=1, size=3.2, 
            segment.size=0.25, segment.alpha=0.6
         ) +
         
         scale_x_continuous(limits=xlims) +
         scale_y_continuous(limits=ylims, labels=superscriptNotation) +
         
         scale_fill_distiller(labels=signif_bin_labels,name='Fisher signif.', palette='Spectral') +
         scale_size_continuous(labels=signif_bin_labels,name='Fisher signif.', range=c(0.1,3)) +
         
         labs(x='Gene name',y='q-value fisher') +
         guides(size=guide_legend(reverse=T), fill=guide_legend(reverse=T)) +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            axis.text.y=element_text(size=10),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.text=element_text(hjust=0)
         )
   } else {
      plot <- ggplot(df, aes(n_rel, p_rel, label=label)) +
         geom_abline(slope=1, intercept=0, linetype='dotted', size=0.25) +
         geom_point(aes(size=signif_bin, fill=signif_bin), shape=21, stroke=0.3) +
         geom_text_repel(
            point.padding=0.25, nudge_x=0.03, min.segment.length=0.7, size=3.2, 
            segment.size=0.25, segment.alpha=0.6
         ) +
         
         scale_x_continuous(limits=xlims, breaks=seq(0,1,0.1), labels=scales::percent_format(accuracy=1)) +
         scale_y_continuous(limits=ylims, breaks=seq(0,1,0.1), labels=scales::percent_format(accuracy=1)) +
         
         scale_fill_distiller(labels=signif_bin_labels,name='Fisher signif.', palette='Spectral') +
         scale_size_continuous(labels=signif_bin_labels,name='Fisher signif.', range=c(0.1,3)) +
         guides(size=guide_legend(reverse=T), fill=guide_legend(reverse=T)) +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            axis.title=element_text(size=10),
            legend.text=element_text(hjust=0)
         )
   }
   
   if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
   if(!is.null(x.lab)){ plot <- plot + xlab(x.lab) }
   if(!is.null(y.lab)){ plot <- plot + ylab(y.lab) }
   
   if(legend.in.plot){
      plot <- plot + theme(
         legend.position=c(0.9, 0.9),
         legend.justification=c(1,1),
         legend.box.background=element_rect(fill=NA, color='black')
      )
   }

   return(plot)
}

## Export
plotFisherTestWrapper <- function(data, export.path, ...){
   #data=subset(hr_gene_likelihood[nrow(hr_gene_likelihood):1,], keep!=-1)
   p <- plotFisherTest(
      data,
      order.rows.by.gene.name=F,
      #plot.title='Gene deficiency enrichment in HRD samples',
      rel.freq.on.xy=T, y.lab='% HRD patients with gene deficiency', x.lab='% HRP patients with gene deficiency',
      ylims=c(0, 0.45), xlims=c(0, 0.3)
   )
   
   # guide_legend_custom <- guide_legend(
   #    title.position='bottom', title.hjust=0.5, 
   #    label.position='bottom', label.theme=element_text(angle=90, hjust=1, vjust=0.5, size=9),
   #    nrow=1
   # )
   # 
   # p <- p + 
   #    theme(legend.position='bottom') +
   #    guides(size=guide_legend_custom, fill=guide_legend_custom)
   
   pdf(export.path, 4.75, 3.75)
   grid.draw(p)
   dev.off()
}

## Blacklist STARD13, NF1 
hr_gene_likelihood$keep[hr_gene_likelihood$hgnc_symbol %in% c('STARD13','NF1')] <- -1

plotFisherTestWrapper(
   subset(hr_gene_likelihood, keep!=-1),
   paste0(script_dir,'/plots/hr_gene_likelihood_fisher_filt.pdf')
)

plotFisherTestWrapper(
   hr_gene_likelihood,
   paste0(script_dir,'//plots/hr_gene_likelihood_fisher.pdf')
)


#========= Co-occurrence of LOH of BRCA2/STARD13 and BRCA1/NF1 =========#
#--------- Fisher ---------#
mkCnvCooccMat <- function(diplotypes, gene1, gene2){
   # diplotypes=diplotypes_filt
   # gene1='BRCA2'
   # gene2='STARD13'
   
   genes <- c(gene1,gene2)
   
   l <- lapply(genes, function(i){
      df <- diplotypes[diplotypes$hgnc_symbol==i,c('sample','a1.eff')]
      colnames(df)[2] <- i
      return(df)
   })
   
   df <- Reduce(function(x,y){ merge(x,y,by='sample',all=T) }, l)
   df[is.na(df)] <- 'none'
   rownames(df) <- df$sample; df$sample <- NULL
   
   cn_effects <- c('deep_deletion','loh')

   df[,1] <- as.integer(df[,1] %in% cn_effects)
   df[,2] <- as.integer(df[,2] %in% cn_effects)
   
   df <- df[rev(do.call(order, df)),]
   
   return(df)
}
#mkCnvCooccMat(diplotypes_filt, 'BRCA2', 'STARD13')

fisherTest1VsAll <- function(diplotypes, constant.gene, fisher.cutoff=0.1){
   #constant.gene='BRCA2'
   
   genes <- unique(diplotypes$hgnc_symbol)
   genes <- genes[genes!=constant.gene]
   
   pb <- txtProgressBar(max=length(genes), style=3)
   counter<-0
   l <- lapply(genes, function(i){
      counter <<- counter+1
      setTxtProgressBar(pb, counter)
      
      #i='STARD13'
      df <- mkCnvCooccMat(diplotypes, constant.gene, i)
      
      p <- sum(df[[1]]==1 & df[[2]]==1)
      p_total <- sum(df[[1]]==1)
      
      n <- sum(df[[1]]==0 & df[[2]]==1)
      n_total <- sum(df[[1]]==0)
      
      tab <- rbind(
         c(p, n),
         c(p_total, n_total),
         deparse.level=0
      )
      
      pvalue <- fisher.test(tab, alternative='greater')$p.value
      
      data.frame(
         hgnc_symbol=i,
         p, p_total, p_rel=p/p_total,
         n, n_total, n_rel=n/n_total,
         pvalue
      )
   })
   
   df <- do.call(rbind, l)
   df <- df[order(df$pvalue),]
   df$qvalue <- p.adjust(df$pvalue,'hochberg')
   
   df$keep <- as.integer(df$qvalue < fisher.cutoff)
   
   class(df) <- c(
      cutoff=fisher.cutoff,
      # n_hrd=n_hrd,
      # n_hrp=n_hrp,
      # n_samples=n_samples,
      data_type='data.frame'
   )
   
   return(df)
}

#--------- Main ---------#
gene_chroms <- c('17'='BRCA1','13'='BRCA2')

path_cnv_coocc <- paste0(script_dir,'/data/cnv_coocc.rds')
if(file.exists(path_cnv_coocc)){
   cnv_coocc <- readRDS(path_cnv_coocc)
} else {
   cnv_coocc <- lapply(gene_chroms, function(i){
      message('\nDetermining 1 vs. all co-occurrences of CNVs for: ', i)
      fisherTest1VsAll(diplotypes_filt, i)
   })
   names(cnv_coocc) <- gene_chroms
   saveRDS(cnv_coocc, path_cnv_coocc)
}

#--------- Mark Chr 13/17 and plot ---------#
detIsInChrom <- function(genes, bed, target.chrom){
   # genes=df$hgnc_symbol
   # bed=genes_bed
   # target.chrom='17'
   
   chrom <- bed$chrom[ match(genes, bed$hgnc_symbol) ]
   chrom[is.na(chrom)] <- '0'
   as.integer(chrom==target.chrom)
}

plotFisherTestCnvCoocc <- function(
   df, fisher.cutoff=0.1, 
   signif.cutoffs=10^-c(0:5), legend.sci.nota.cutoff=0.01,
   plot.title=NULL, x.lab=NULL, y.lab=NULL, rel.freq.on.xy=F, legend.in.plot=F
){
   #df=hr_gene_likelihood
   # df=cn_eff_coocc[[1]]
   # signif.cutoffs=10^-c(0,100,200,300)
   # fisher.cutoff=10^-100
   
   ## Transform
   df <- within(df,{
      pvalue_log <- -log10(pvalue)
      qvalue_log <- -log10(qvalue)
   })
   
   ## Labels
   df <- within(df,{
      label <- hgnc_symbol
      label[!(qvalue < fisher.cutoff)] <- ''
   })
   
   ## Signif bins
   signif.cutoffs <- sort(signif.cutoffs, decreasing=T)
   cutoffs <- -log10(signif.cutoffs)
   
   df$signif_bin <- .bincode(df$qvalue_log, breaks=c(cutoffs, Inf))
   df$signif_bin[is.na(df$signif_bin)] <- 1
   
   #signif_bin_labels <- signif.cutoffs[df$signif_bin]
   signif_bin_labels <- sapply(format(signif.cutoffs,scientific=T), function(i){
      #i="1e-02"
      i_numeric <- as.numeric(i)
      if(i_numeric >= legend.sci.nota.cutoff){
         ## Regular notation
         out <- format(i_numeric, scientific=F, drop0trailing=T)
      } else {
         ## Exponent notation
         out <- sub('1e[+-]','10^-',i)
         out <- sub('[-]0+','-',out)
      }
      return(out)
   }, USE.NAMES=F)
   signif_bin_labels <- paste0('q<',signif_bin_labels)
   signif_bin_labels <- parse(text=signif_bin_labels)
   
   ##
   df <- df[order(df$hgnc_symbol),]
   df$index <- 1:nrow(df)
   
   superscriptNotation <- function(x){
      #x=c(0, 10, 20)
      out <- paste0('10^-',x)
      out <- gsub('^10\\^-0$', 1, out) ## Change log10(0) to 1
      parse(text=out)
   }
   
   ## Main
   if(!rel.freq.on.xy){
      plot <- ggplot(df, aes(index, qvalue_log, label=label, color=is_in_chrom)) +
         geom_hline(yintercept=-log10(fisher.cutoff),linetype='dotted') +
         geom_point(aes(size=signif_bin, fill=signif_bin), shape=21, stroke=0.3) +
         geom_text_repel(
            point.padding=0.25, nudge_y=0.6, min.segment.length=1, size=3, 
            segment.size=0.25, segment.alpha=0.6
         ) +
         
         scale_x_continuous() +
         scale_y_continuous(labels=superscriptNotation) +
         
         scale_fill_distiller(labels=signif_bin_labels,name='Fisher signif.', palette='Spectral') +
         scale_size_continuous(labels=signif_bin_labels,name='Fisher signif.', range=c(0.5,5)) +
         
         labs(x='Gene name',y='q-value fisher') +
         guides(
            size=guide_legend(reverse=T, override.aes=list(color='black')), 
            fill=guide_legend(reverse=T),
            color=F
         ) +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            axis.text.y=element_text(size=10),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.text=element_text(hjust=0)
         )
   } else {
      plot <- ggplot(df, aes(n_rel, p_rel, label=label, color=is_in_chrom)) +
         geom_abline(slope=1, intercept=0, linetype='dotted', size=0.25) +
         geom_point(aes(size=signif_bin, fill=signif_bin), shape=21, stroke=0.3) +
         geom_text_repel(
            point.padding=0.25, nudge_x=0.03, min.segment.length=0.7, size=2.5, 
            segment.size=0.25, segment.alpha=0.6
         ) +
         
         scale_x_continuous(breaks=seq(0,1,0.1), labels=function(x){ paste0(x*100, '%') }) +
         scale_y_continuous(breaks=seq(0,1,0.1), labels=function(x){ paste0(x*100, '%') }) +
         
         scale_fill_distiller(labels=signif_bin_labels,name='Fisher signif.', palette='Spectral') +
         scale_size_continuous(labels=signif_bin_labels,name='Fisher signif.', range=c(0.25,3)) +
         guides(
            size=guide_legend(reverse=T, override.aes=list(color='black')), 
            fill=guide_legend(reverse=T),
            color=F
         ) +
         
         theme_bw() +
         theme(
            plot.title=element_text(hjust=0.5, size=11),
            axis.title=element_text(size=10),
            legend.text=element_text(hjust=0)
         )
   }
   
   if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
   if(!is.null(x.lab)){ plot <- plot + xlab(x.lab) }
   if(!is.null(y.lab)){ plot <- plot + ylab(y.lab) }
   
   if(legend.in.plot){
      plot <- plot + theme(
         legend.position=c(0.9, 0.9),
         legend.justification=c(1,1),
         legend.box.background=element_rect(fill=NA, color='black')
      )
   }
   
   return(plot)
}

#plotFisherTest(cnv_coocc$BRCA1, rel.freq.on.xy=T, signif.cutoffs=10^-c(0,100,200,300), fisher.cutoff=10^-150)

plots_cnv_coocc <- lapply(1:length(gene_chroms), function(i){
   gene=gene_chroms[[i]]
   chrom=names(gene_chroms[i])
   
   df <- cnv_coocc[[gene]]
   df$is_in_chrom <- detIsInChrom(df$hgnc_symbol, genes_bed, chrom)
   df$is_in_chrom <- factor(df$is_in_chrom, levels='1','0')
   
   # plot <- plotFisherTestCnvCoocc(
   #    df, x.lab='Gene2 name',
   #    signif.cutoffs=10^-c(0,100,200,300), fisher.cutoff=10^-150
   # )
   
   plot <- plotFisherTestCnvCoocc(
      df, signif.cutoffs=10^-c(0,100,200,300), fisher.cutoff=10^-150,
      rel.freq.on.xy=T,
      y.lab=paste0('Freq. of chromosomal alteration occurring\nin both Gene2 and ',gene),
      x.lab=paste0('Freq. of chromosomal alteration occurring\nin Gene2 but not ',gene)
   )
   
   plot <- plot + 
      ggtitle(
         sprintf('Co-occurrence of chromosomal alterations (Deep deletion or LOH): %s~Gene2', gene),
         sprintf('(Red: Gene2 on Chr. %s)', chrom)
      ) +
      theme(
         plot.title=element_text(size=11, hjust=0.5),
         plot.subtitle=element_text(size=9, hjust=0.5)
      )
   
   return(plot)
})

pdf(paste0(script_dir,'/plots/cnv_coocc_fisher.pdf'), 10, 6)
for(i in plots_cnv_coocc){ 
   suppressWarnings({ grid.draw(i) })
}
dev.off()


####################################################################################################
# Export final output                                                                              #
####################################################################################################
#--------- Fisher table ---------#
write.table(
   hr_gene_likelihood, paste0(script_dir,'/data/hr_gene_likelihood.txt'),
   sep='\t',row.names=F,quote=F
)

hr_gene_likelihood_exp <- hr_gene_likelihood
colnames(hr_gene_likelihood_exp)[
   colnames(hr_gene_likelihood_exp) %in% c('p','p_total','p_rel','n','n_total','n_rel')
   ] <- c('n_hrd_def','n_hrd','n_hrd_def_rel','n_hrp_def','n_hrp','n_hrp_def_rel')

write.table(hr_gene_likelihood_exp, paste0(script_dir,'/data/hr_gene_likelihood_exp.txt'), sep='\t',row.names=F,quote=F)

#--------- Diplotypes subset ---------#
hr_genes_sel <- hr_gene_likelihood[hr_gene_likelihood$keep==1,'hgnc_symbol']
diplotypes_hrGenes  <- (function(){
   df <- diplotypes_filt
   df <- df[df$hgnc_symbol %in% hr_genes_sel,]
   return(df)
})()
write.table(
   diplotypes_hrGenes, gzfile(paste0(script_dir,'/data/diplotypes_hrGenes.txt.gz')),
   sep='\t',row.names=F,quote=F
)

diplotypes_hrd_hrGenes  <- (function(){
   df <- diplotypes_filt
   df <- df[df$sample %in% hrd_samples & df$hgnc_symbol %in% hr_genes_sel,]
   return(df)
})()
write.table(
   diplotypes_hrd_hrGenes, gzfile(paste0(script_dir,'/data/diplotypes_hrd_hrGenes.txt.gz')),
   sep='\t',row.names=F,quote=F
)


















