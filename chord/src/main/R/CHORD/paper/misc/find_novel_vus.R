library(reshape2)

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
## Metadata
metadata <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/metadata_for_analysis.txt'))

## Preds
pred <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/pred_analysis_samples.txt'))
hrp_samples <- pred[pred$hrd < 0.5,'sample']

## Diplotypes
diplotypes_hrd_hrGenes <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/diplotypes_hrd_hrGenes.txt.gz'))

diplotypes_raw <- read.delim(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/gene_diplotypes_with_amps_HMF_PCAWG.txt.gz'))
#diplotypes_raw <- read.delim('/Users/lnguyen/Documents/R_cache/gene_diplotypes_max_hmf_pcawg.txt.gz')


## 
l_m_diplotypes <- readRDS(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/l_m_diplotypes.rds'))
rank_order_clust <- readRDS(paste0(base_dir,'/CHORDv2/analysis/pancancer_overview/data/rank_order_clust.rds'))

#========= Narrow down VUS's =========#
splitDiplotypeMatrix <- function(m){
   l <- list(
      a1=m[,grep('_1$',colnames(m))],
      a2=m[,grep('_2$',colnames(m))]
   )
}

meltDiplotypeMatrix <- function(l_m_diplotypes){
   meltDiplotypeMatrixAllele <- function(allele='a2'){
      counter <- 0
      l <- lapply(l_m_diplotypes, function(i){
         counter <<- counter + 1
         df <- as.data.frame(splitDiplotypeMatrix(i)[[allele]])
         df <- cbind(sample=rownames(df), df)
         df_melt <- melt(df, 'sample')
         
         df_melt$variable <- sapply(strsplit(as.character(df_melt$variable),'_'),`[`,1)
         colnames(df_melt)[2] <- 'hgnc_symbol'
         colnames(df_melt)[3] <- names(l_m_diplotypes)[counter]
         return(df_melt)
      })
      
      Reduce(function(x,y){ merge(x,y,sort=F) },l)
   }
   
   list(
      a1=meltDiplotypeMatrixAllele('a1'),
      a2=meltDiplotypeMatrixAllele('a2')
   )
}

getDiplotypesWithNovelVus <- function(l_m_diplotypes, rank_order_clust, diplotypes){
   
   #--------- Make filters ---------#
   ## Convert diplotype matrix back to diplotype table
   l_m_diplotypes_melt <- meltDiplotypeMatrix(l_m_diplotypes)
   
   df_selection <- data.frame(
      sample=l_m_diplotypes_melt$a1$sample,
      hgnc_symbol=l_m_diplotypes_melt$a1$hgnc_symbol,
      a1.eff=l_m_diplotypes_melt$a1$eff,
      a2.eff=l_m_diplotypes_melt$a2$eff,
      a1.max_score=l_m_diplotypes_melt$a1$score,
      a2.max_score=l_m_diplotypes_melt$a2$score,
      a2.max_score_origin=l_m_diplotypes_melt$a2$max_score_origin
   )
   
   df_selection$has_new_pathogenic <- with(df_selection,{
      # a1.max_score >= 5 & a2.max_score >= 3 & 
      #    (a2.max_score < 5 | (a2=='frameshift' & a2.max_score_origin!='known'))
      a1.max_score >= 5 & a2.max_score >= 3 & a2.max_score < 5
   })
   
   ## Get cluster per sample
   sample_clusters <- rank_order_clust$clusters
   cluster_names <- rank_order_clust$cluster_names
   
   for(i in cluster_names){
      sample_clusters[sample_clusters==i] <- names(cluster_names)[i]
   }
   
   df_selection$cluster <- sample_clusters
   
   ## Subset for potentially new pathogenic variants
   df_selection <- df_selection[df_selection$has_new_pathogenic,]
   
   
   #--------- Get diplotypes ---------#
   diplotypes_ss <- do.call(rbind,lapply(1:nrow(df_selection), function(i){
      row <- df_selection[i,]
      
      
      variant_type_selection <- row$a2.eff
      
      if(variant_type_selection=='frameshift'){ variant_type_selection <- 'frameshift_variant' }
      else if(variant_type_selection=='essential_splice_variant'){ variant_type_selection <- c('splice_acceptor_variant','splice_donor_variant') }
      else if(variant_type_selection=='nonsense'){ variant_type_selection <- c('start_lost','stop_gained','stop_lost') }
      else if(variant_type_selection=='missense'){ variant_type_selection <- 'missense_variant' }
      
      diplotypes_raw[
         diplotypes_raw$sample==row$sample 
         & diplotypes_raw$hgnc_symbol==row$hgnc_symbol
         & diplotypes_raw$a2.eff %in% variant_type_selection
         & !(diplotypes_raw$diplotype_origin %in% c('germ_som','som_som'))
      ,]
      
   }))
   
   
   
   diplotypes_ss$cluster <- sample_clusters[ match(diplotypes_ss$sample,names(sample_clusters)) ]
   diplotypes_ss <- diplotypes_ss[diplotypes_ss$cluster==diplotypes_ss$hgnc_symbol,]
      
   return(diplotypes_ss)
}

diplotypes_vus <- getDiplotypesWithNovelVus(l_m_diplotypes, rank_order_clust, diplotypes_hrd_hrGenes)


#========= Get variants from raw data using diplotypes_vus as a guide =========#
#gene_ann_dir <- paste0(base_dir,'/datasets/HMF_DR010_DR047/gene_ann/')

mut_profile_paths <- read.table(paste0(base_dir,'/CHORDv2/analysis/find_novel_vus/data/mut_profile_paths.txt'), header=F)

colnames(mut_profile_paths)[1] <- 'path'
mut_profile_paths$sample <- sapply(strsplit(mut_profile_paths$path,'/'),`[[`,10)

mut_profile_paths$origin <- gsub(
   '.txt.gz$','',
   sapply(strsplit(basename(mut_profile_paths$path),'_'),`[[`,3)
)

if(dir.exists('/Users/lnguyen/')){
   mut_profile_paths$path <- paste0('/Users/lnguyen/',mut_profile_paths$path)
}

getVusFromDiplotype <- function(
   diplotypes,
   #ann.data.dir=gene_ann_dir,
   mut.profile.paths=mut_profile_paths
){
   # diplotypes=diplotypes_vus
   # ann.data.dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/HMF_DR010_DR047/vcf_subset/'
   
   #--------- Search through raw data ---------#
   l <- with(diplotypes,{
      Map(function(sample, ensembl_gene_id, diplotype_origin, a2.max_score){
         
         a2_origin <- strsplit(diplotype_origin,'_')[[1]][2]
         #ann_data_path <- paste0(ann.data.dir,'/',sample,'/mut_profiles/mut_profile_',a2_origin,'.txt.gz')
         ann_data_path <- mut.profile.paths[
            mut.profile.paths$sample==sample & mut.profile.paths$origin==a2_origin,
            'path'
         ]
         
         print(ann_data_path)
         
         df <- read.delim(ann_data_path)
         df_ss <- df[df$ensembl_gene_id==ensembl_gene_id & df$max_score==a2.max_score,]
         df_ss$mut_origin <- a2_origin
         
         return(df_ss)
         
      }, sample, ensembl_gene_id, diplotype_origin, a2.max_score)
   })
   
   #--------- Merge and format columns/values ---------#
   out <- do.call(rbind, l)
   sel_cols <- c(
      'chrom'='chrom',
      'pos'='pos',
      'ref'='ref',
      'alt'='alt',
      'hgvs_c'='hgvs_c',
      'snpeff_eff'='variant_type',
      'mut_origin'='mutation_origin',
      'hgnc_symbol'='hgnc_symbol',
      'ensembl_gene_id'='ensembl_gene_id',
      'clinvar_sig'='clinvar_annotation'
   )
   out <- out[names(sel_cols)]
   
   out <- within(out,{
      clinvar_sig[is.na(clinvar_sig)] <- ''
      mut_origin[mut_origin=='germ'] <- 'germline'
      mut_origin[mut_origin=='som'] <- 'somatic'
   })
   
   colnames(out) <- sel_cols
   
   # out$Comment <- ''
   # out$Comment[ grep('benign',out[,'ClinVar annotation'], ignore.case=T) ] <- 'Was highest impact 2nd hit, even though variant is annotated as benign'
   
   out <- cbind(
      sample=sapply(strsplit(rownames(out),'[.]'),`[[`,1),
      out
   )
   rownames(out) <- NULL
   return(out)
}

mut_profile_vus <- getVusFromDiplotype(diplotypes_vus)


#=========  =========#
novel_vus_strings <- with(mut_profile_vus,{ paste(chrom,pos,ref,alt,sep=':') })
#gene_ann_subdirs <- list.dirs(gene_ann_dir,recursive=F,full.names=T)

###
novel_vus_in_hrp_samples_path <- paste0(base_dir,'/CHORDv2/analysis/find_novel_vus/data/novel_vus_in_hrp_samples.rds')

if(file.exists(novel_vus_in_hrp_samples_path)){
   novel_vus_in_hrp_samples <- readRDS(novel_vus_in_hrp_samples_path)
} else {
   progress_counter <- 0
   l <- split(mut_profile_paths,mut_profile_paths$sample)
   novel_vus_in_hrp_samples <- lapply(l, function(i){
      #i=l[[1]]
      
      progress_counter <<- progress_counter + 1
      #if(progress_counter<=5){ break }
      
      sample_name <- i[1,'sample']
      if(!(sample_name %in% hrp_samples)){ return(NULL) }
      
      message('Processing [',progress_counter,']: ', sample_name)
      
      # txt_files <- c(
      #    germ=paste0(i,'/mut_profiles/mut_profile_germ.txt.gz'),
      #    som=paste0(i,'/mut_profiles/mut_profile_som.txt.gz')
      # )
      
      txt_files <- c(
         germ=i[i$origin=='germ','path'],
         som=i[i$origin=='som','path']
      )
      
      counter <- 0
      variants <- do.call(rbind,lapply(txt_files, function(i){
         counter <<- counter + 1
         df <- read.delim(i)
         df <- df[,c('chrom','pos','ref','alt','hgvs_c','snpeff_eff')]
         df$mut_origin <- names(txt_files)[counter]
         return(df)
      }))
      rownames(variants) <- NULL
      
      variant_strings <- with(variants,{ paste(chrom,pos,ref,alt,sep=':') })
      
      variants[variant_strings %in% novel_vus_strings,]
   })
   
   saveRDS(
      novel_vus_in_hrp_samples,
      paste0(base_dir, '/CHORDv2/analysis/find_novel_vus/data/novel_vus_in_hrp_samples.rds')
   )
}

##! None of the novel VUS's were found in HRP samples. No further action necessary
table(
   sapply(novel_vus_in_hrp_samples, function(i){
      if(is.null(i)){ return(0) }
      else { nrow(i) }
   })
)


#========= Format output =========#
mut_profile_vus_export <- (function(){
   df <- mut_profile_vus

   df <- cbind(
      cluster=rank_order_clust$clusters[ match(df$sample, names(rank_order_clust$clusters)) ],
      df
   )

   df <- cbind(
      sample_index=match(df$sample, names(rank_order_clust$clusters)),
      df
   )
   
   sample_names <- df$sample
   df$sample <- metadata$hmf_id[ match(df$sample, metadata$sample) ]
   df$sample[is.na(df$sample)] <- sample_names[is.na(df$sample)] 
   
   df <- df[order(df$sample_index),]

   return(df)
})()

write.table(
   mut_profile_vus_export,
   paste0(base_dir,'/CHORDv2/analysis/find_novel_vus/data/likely_pathogenic_vus.txt'),
   sep='\t',quote=F,row.names=F
)











