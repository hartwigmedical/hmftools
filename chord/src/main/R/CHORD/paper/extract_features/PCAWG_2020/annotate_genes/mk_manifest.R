options(stringsAsFactors=F)

#========= Path prefixes =========#
base_dir <- list(
   hpc='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/',
   mnt='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/',
   local='/Users/lnguyen/Documents/Luan_projects/datasets/'
)

for(i in base_dir){
   if(dir.exists(i)){ 
      base_dir <- i 
      break
   }
}

write.tsv <- function(x,file,...){
   write.table(x,file,sep='\t',row.names=F,quote=F,...)
}

#========= Main =========#
path_files <- list.files(paste0(base_dir,'/PCAWG_2020/scripts/annotate_genes/paths/'), full.names=T)

paths <- lapply(path_files, function(i){
   df <- read.table(i)
   colnames(df) <- 'path'
   df$sample_id <- sapply(strsplit(basename(df$path),'[.]'),`[[`,1)
   return(df)
})
names(paths) <- sapply(strsplit(basename(path_files),'[.]'),`[[`,1)

paths <- paths[c('germ_smnv_indel','som_smnv_indel','som_cnv')]

manifest <- Reduce(function(x,y){
   merge(x,y,by='sample_id',all=T)
}, paths)
colnames(manifest)[2:ncol(manifest)] <- names(paths)

manifest <- manifest[order(manifest$sample_id),]

manifest_all_data <- na.exclude(manifest)
write.tsv(manifest_all_data, paste0(base_dir,'/PCAWG_2020/scripts/annotate_genes/manifest_samples_with_all_data.txt'))

