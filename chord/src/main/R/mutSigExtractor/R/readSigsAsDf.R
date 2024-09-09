#' Read multiple files outputed my mutSigExtractor into a data frame
#'
#' @param dir A directory from which to read the signatures/contexts
#' @param files Alternatively, a list of file paths can also be specified
#' @param samples.as.rows Samples as rows and signatures/contexts as columns?
#' @param verbose Show file reading progress?
#' @param ... Arguments that can be pass to list.files() and read.table()
#'
#' @return A data frame
#' @export
#'
#' @examples
#' base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Breast_Organoids/matrices/'
#' contexts <- list(
#'    snv = readSigsAsDf(paste0(base_dir,'/snv_contexts')),
#'    indel = readSigsAsDf(paste0(base_dir,'/indel')),
#'    sv = readSigsAsDf(paste0(base_dir,'/sv_contexts'))
#' )

readSigsAsDf <- function(dir = NULL, files = NULL, samples.as.rows = T, verbose = T, ...){
   if(!is.null(dir)){
      files <- list.files(dir, full.names = T, ...)
   }

   if(verbose){ n_files <- length(files) }
   m <- do.call(cbind, lapply(1:n_files, function(i){
      file <- files[i]
      if(verbose){ message('Reading ', i, '/', n_files, ': ', basename(file)) }
      read.table(file, check.names = F, ...)
   }))

   if(samples.as.rows){ m <- t(m) }

   return(as.data.frame(m))
}
