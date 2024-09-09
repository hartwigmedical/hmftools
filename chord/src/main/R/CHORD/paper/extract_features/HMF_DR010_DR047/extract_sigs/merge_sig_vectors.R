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

#========= Main =========#
out_dir <- paste0(base_dir,'/datasets/HMF_DR010_DR047/matrices/')
sub_dirs <- list.dirs(out_dir, recursive=F)

l_paths <- lapply(sub_dirs, list.files, full.names=T)
names(l_paths) <- basename(sub_dirs)

# l_sample_names <- lapply(l_paths,basename)
# l_sample_names <- lapply(l_sample_names, function(i){
#    sapply(strsplit(i,'_'),`[[`,1)
# })

readContextsFromPaths <- function(paths){
   counter <- 0
   pb <- txtProgressBar(max=length(paths), style=3)
   l <- lapply(paths, function(i){
      counter <<- counter+1
      setTxtProgressBar(pb, counter)
      read.delim(i)
   })
   message('\n')
   return(l)
}

ll_contexts <- list()
ll_contexts$snv <- readContextsFromPaths(l_paths$snv)
ll_contexts$dbs <- readContextsFromPaths(l_paths$dbs)
ll_contexts$indel <- readContextsFromPaths(l_paths$indel)
ll_contexts$sv <- readContextsFromPaths(l_paths$sv)

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

