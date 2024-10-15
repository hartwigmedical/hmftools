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

#========= Main =========#
out_dir <- paste0(base_dir,'/datasets/processed/PCAWG_2020/matrices/')
#sub_dirs <- list.dirs(out_dir, recursive=F)
sub_dirs <- paste0(out_dir,'/',c('snv','dbs','indel','sv'))

l_paths <- lapply(sub_dirs, list.files, full.names=T)
names(l_paths) <- basename(sub_dirs)

# l_sample_names <- lapply(l_paths,basename)
# l_sample_names <- lapply(l_sample_names, function(i){
#    sapply(strsplit(i,'_'),`[[`,1)
# })

readContextsFromPaths <- function(paths){

   counter <- 0
   pb <- txtProgressBar(max=length(paths), style=3)
   #paths=l_paths$snv[1:2]
   lapply(paths, function(i){
      counter <<- counter+1
      setTxtProgressBar(pb, counter)
      read.delim(i, check.names=F)
   })
   #message('\n')
}

# ll_contexts <- list()
# ll_contexts$snv <- readContextsFromPaths(l_paths$snv)
# ll_contexts$indel <- readContextsFromPaths(l_paths$indel)
# ll_contexts$sv <- readContextsFromPaths(l_paths$sv)
ll_contexts <- lapply(l_paths, readContextsFromPaths)

## Merge into one matrix per mut type
l_contexts <- lapply(ll_contexts, function(i){
   t(do.call(cbind,i))
})

## Export
saveRDS(l_contexts, paste0(out_dir,'/contexts_merged.rds'))
write.table(
   do.call(cbind,l_contexts), paste0(out_dir,'/contexts_merged.txt'),
   sep='\t',quote=F
)

