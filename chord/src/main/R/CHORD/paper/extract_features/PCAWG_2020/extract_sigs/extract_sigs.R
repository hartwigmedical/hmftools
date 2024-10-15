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

library(devtools)
load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mutSigExtractor/'))

#========= Main =========#
args <- commandArgs(trailingOnly=TRUE)
out_dir=args[1]
sample_name=args[2]
vcf_snv=args[3]
vcf_indel=args[4]
bedpe_sv=args[5]

manifest <- read.delim(paste0(base_dir,'/datasets/processed/PCAWG_2020/manifest/manifest_all_som_data.txt'))

bedpeToSvContexts <- function(bedpe.file, sample.name){
   #bedpe.file=paste0('/Users/lnguyen/',manifest[1,'sv'])
   bedpe <- read.delim(bedpe.file)

   df <- data.frame(
      sv_type=bedpe$svclass,
      sv_len=abs(bedpe$start1-bedpe$start2)
   )
   df[grep('INV$',df$sv_type),'sv_type'] <- 'INV'

   extractSigsSv(df=df, output='contexts', sample.name=sample.name)
}


write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

#out_dir <- paste0(base_dir,'/datasets/processed//PCAWG_2020/matrices/')
done_dir <- paste0(out_dir,'/done_files/')
mut_types <- c('snv','dbs','indel','sv')

main <- function(sample_name, vcf_snv, vcf_indel, bedpe_sv, ...){

   # i=1
   # sample_name=manifest[i,'sample_id']
   # vcf_snv=paste0('/Users/lnguyen/',manifest[i,'som_smnv'])
   # vcf_indel=paste0('/Users/lnguyen/',manifest[i,'som_indel'])
   # bedpe_sv=paste0('/Users/lnguyen/',manifest[i,'som_sv'])

   out_paths <- paste0(out_dir,'/',mut_types,'/',sample_name,'_',mut_types,'.txt') ## subdirs were already created (manually)
   names(out_paths) <- mut_types

   contexts <- list()

   message('## Extracting SNV contexts...')
   if(!file.exists(out_paths['snv'])){
      contexts$snv <- extractSigsSnv(vcf_snv, output='contexts', sample.name=sample_name, merge.consecutive=T)
      write.tsv(contexts$snv, out_paths['snv'])
   }
   
   message('## Extracting DBS contexts...')
   if(!file.exists(out_paths['dbs'])){
      contexts$dbs <- extractSigsDbs(vcf_snv, output='contexts', sample.name=sample_name, merge.consecutive=T)
      write.tsv(contexts$dbs, out_paths['dbs'])
   }

   message('## Extracting indel contexts...')
   if(!file.exists(out_paths['indel'])){
      contexts$indel <- extractSigsIndel(vcf_indel, sample.name=sample_name)
      write.tsv(contexts$indel, out_paths['indel'])
   }

   message('## Extracting SV contexts...')
   if(!file.exists(out_paths['sv'])){
      contexts$sv <- bedpeToSvContexts(bedpe_sv, sample.name=sample_name)
      write.tsv(contexts$sv, out_paths['sv'])
   }
   
   if(all(file.exists(out_paths))){
      file.create(paste0(done_dir,'/',sample_name,'.done'))
   }
}

main(
   out_dir=args[1],
   sample_name=args[2],
   vcf_snv=args[3],
   vcf_indel=args[4],
   bedpe_sv=args[5]
)


# counter <- 0
# for(i in 1:nrow(manifest)){
#    counter <- counter + 1
# 
#    sample_name = manifest[i,'sample_id']
#    vcf_snv = paste0('/Users/lnguyen/',manifest[i,'snv_mnv'])
#    vcf_indel = paste0('/Users/lnguyen/',manifest[i,'indel'])
#    bedpe_sv = paste0('/Users/lnguyen/',manifest[i,'sv'])
# 
#    message('[',counter,']: ',sample_name)
# 
#    main(
#       sample_name = sample_name,
#       vcf_snv = vcf_snv,
#       vcf_indel = vcf_indel,
#       bedpe_sv = bedpe_sv,
#    )
# }

#========= Exec =========#











