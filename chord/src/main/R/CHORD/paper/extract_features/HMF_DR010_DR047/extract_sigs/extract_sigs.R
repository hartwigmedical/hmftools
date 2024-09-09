options(stringsAsFactors=F)

#========= Path prefixes =========#
base_dir <- (function(suffix='/Luan_projects/'){
   paths <- list(
      hpc='/hpc/cog_bioinf/cuppen/project_data/',
      mnt='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/',
      local='/Users/lnguyen/Documents/'
   )

   for(i in paths){
      if(dir.exists(i)){
         path <- i
         break
      }
   }

   if(!is.null(suffix)){ path <- paste0(path,suffix) }
   if(!file.exists(path)){ warning('Path does not exist: ',path) }

   return(path)
})()


library(devtools)
load_all(paste0(base_dir,'/CHORD/scripts_main/mutSigExtractor/'))

#========= Main =========#
write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

args <- commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
vcf_som <- args[2]
vcf_sv <- args[3]

out_dir <- paste0(base_dir,'/datasets/HMF_DR010_DR047/matrices/')
mut_types <- c('snv','dbs','indel','sv')
for(i in paste0(out_dir,mut_types,'/')){
   dir.create(i, recursive=T, showWarnings=F)
}


main <- function(sample_name, vcf_som, vcf_sv){

   out_paths <- paste0(out_dir,'/',mut_types,'/',sample_name,'_',mut_types,'.txt') ## subdirs were already created (manually)
   names(out_paths) <- mut_types

   contexts <- list()

   message('## Extracting SNV contexts...')
   if(!file.exists(out_paths['snv'])){
      contexts$snv <- extractSigsSnv(vcf_som, output='contexts', vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$snv, out_paths['snv'])
   }

   message('## Extracting SV contexts...')
   if(!file.exists(out_paths['sv'])){
      contexts$sv <- extractSigsSv(vcf_sv, output='contexts', vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$sv, out_paths['sv'])
   }

   message('## Extracting indel contexts...')
   if(!file.exists(out_paths['indel'])){
      contexts$indel <- extractSigsIndel(vcf_som, vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$indel, out_paths['indel'])
   }

   message('## Extracting DBS contexts...')
   if(!file.exists(out_paths['dbs'])){
      contexts$dbs <- extractSigsDbs(vcf_som, output='contexts', vcf.filter='PASS', sample.name=sample_name)
      write.tsv(contexts$dbs, out_paths['dbs'])
   }

}

main(sample_name, vcf_som, vcf_sv)

#========= Local execution =========#

vcf_paths <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/HMF_DR010_DR047/manifest/manifest.txt')
for(i in 1:nrow(vcf_paths)){
   message('Processing [',i,']: ',vcf_paths$sample[i])
   #i=229
   main(
      vcf_paths$sample[i],
      paste0('/Users/lnguyen/',vcf_paths$som[i]),
      paste0('/Users/lnguyen/',vcf_paths$sv[i])
   )
}
