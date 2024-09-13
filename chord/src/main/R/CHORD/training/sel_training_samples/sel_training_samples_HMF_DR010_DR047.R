library(mutSigExtractor)

options(stringsAsFactors = F)

#========= Paths =========#
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
metadata <- read.delim(paste0(base_dir,'/datasets/HMF_DR010_DR047/metadata/DR-047-update1_metadata.txt'), check.names=F, stringsAsFactors=T)
colnames(metadata) <- sub('#','',colnames(metadata))

brca_status <- read.delim(paste0(base_dir,'/datasets/HMF_DR010_DR047/scripts/annotate_genes/gene_diplotypes_with_amps_HMF_DR010_DR047_BRCA.txt.gz'))
contexts <- readRDS(paste0(base_dir,'/datasets/HMF_DR010_DR047/matrices/contexts_merged.rds'))


#========= Identify samples with MSI =========#
## Based on the MYSQL MSI query, the sample with MSI with the lowest contrib of indel.rep had
## indel.rep==14087. Furthermore, all samples above this value had MSI
## (100% TP rate). MSI was validated at HMF using a PCR based method (?).
indel_load <- transformContexts(indel=contexts$indel, simplify.types='indel')
indel_rep_load <- rowSums(indel_load[,c('ins.rep','del.rep')])
msi_samples <- names(indel_rep_load)[indel_rep_load>=14000]

#========= Main =========#
detTrainingSet <- function(df, genes=c('BRCA2','BRCA1'), msi.samples){

   # df=brca_status
   # genes=c('BRCA2','BRCA1')
   # msi.samples=msi_samples

   ## Determine pre-conditions
   has_msi <- df$sample %in% msi.samples

   is_def <- df$hit_score >= 10
   is_prof <- unlist(with(df,{
      Map(function(diplotype_origin, a1.max_score, a2.max_score){
         if(diplotype_origin=='cnv_germ'){
            a1.max_score==0 & a2.max_score <= 3

         } else if(diplotype_origin=='cnv_som'){
            a1.max_score==0 & a2.max_score <= 3

         } else if(diplotype_origin %in% c('germ_som','som_som')){
            a1.max_score <= 3 & a2.max_score <= 3
         } else {
            FALSE
         }
      }, diplotype_origin, a1.max_score, a2.max_score, USE.NAMES=F)
   }))

   ## Determine in training set (def)
   in_training_set <- rep(FALSE,nrow(df))
   in_training_set[is_def & !has_msi] <- TRUE ## BRCA1/2 class
   #in_training_set[is_prof_confident] <- TRUE ## none class

   ## Combine info into one df
   out_pre <- cbind(
      df[c('sample','hgnc_symbol','hit_score','a1.eff','a1.max_score','a2.eff','a2.max_score','biallelic_status')],
      has_msi,
      is_def,
      is_prof,
      in_training_set
   )

   #--------- Prof (part 1) ---------#
   l <- lapply(genes, function(i){
      out_pre[out_pre$hgnc_symbol==i,c('sample','is_prof')]
   })
   names(l) <- genes

   df_det_is_prof_all <- Reduce(function(x,y){ merge(x,y,by='sample') },l)
   df_det_is_prof_all$is_prof_all <- apply(df_det_is_prof_all[,-1],1,all)

   #--------- Def + make final output ---------#
   ## Sort by genes order
   out <- do.call(rbind, lapply(genes, function(i){
      out_pre[out_pre$hgnc_symbol==i,]
   }))

   ## Then, greedily resolve multiple hit cases
   out <- do.call(rbind, lapply(unique(out$sample),function(i){
      df_ss <- out[out$sample==i,]
      df_ss[which.max(df_ss$is_def),]
   }))

   ## Indicate which gene was deficient
   out$response <- with(out,{
      unlist(Map(function(is_def, hgnc_symbol){
         if(is_def){ hgnc_symbol }
         else{ 'none' }
      }, is_def, hgnc_symbol))
   })

   insertColAfter <- function(df, position, v, name){
      col_index <- which(position==colnames(df))
      df <- cbind(df[,1:col_index], v, df[,(col_index+1):ncol(df)])
      colnames(df)[col_index+1] <- name
      return(df)
   }

   #--------- Prof (part 2) ---------#
   out <- insertColAfter(
      out, 'is_prof',
      df_det_is_prof_all$is_prof_all[match(out$sample, df_det_is_prof_all$sample)],
      'is_prof_all'
   )

   out$in_training_set[out$is_prof_all] <- TRUE

   return(out)
}

annotation <- detTrainingSet(brca_status, genes=c('BRCA2','BRCA1'), msi_samples)
#table(subset(annotation, in_training_set)$response)

#========= Determine max purity biopsy =========#
metadata_ann <- do.call(rbind,lapply(split(metadata,metadata$patient_id), function(i){
   i$max_purity_biopsy <- FALSE
   i$max_purity_biopsy[which.max(i$tumor_purity)] <- TRUE
   return(i)
}))

#========= Finalize =========#
metadata_ann <- cbind(
   metadata_ann,
   annotation[match(metadata_ann$sample_id, annotation$sample),]
)
##metadata_ann$in_training_set[!metadata_ann$max_purity_biopsy] <- FALSE

write.table(
   metadata_ann, paste0(base_dir,'/CHORDv2/training/scripts/sel_training_samples/HMF_DR010_DR047/metadata_training_samples_ann.txt'),
   sep='\t',quote=F, row.names=F
)

