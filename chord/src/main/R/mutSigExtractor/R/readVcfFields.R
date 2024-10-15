#' Read vcf into R as data frame
#'
#' @param vcf.file Path to the vcf file
#' @param fields A character or integer vector indicating which columns to keep.
#'
#' @return A dataframe
#' @export
#'
#' @examples
#' readVcfFields('/path/to/vcf', fields=c('CHROM','POS','REF','ALT'))
#' readVcfFields('/path/to/vcf', fields=c(1,2,4,5))
readVcfFields <- function(vcf.file, fields=NULL){

   ## Scan for the header line
   con  <- file(vcf.file, open = "r")
   line_num <- 0
   while(length(line <- readLines(con, n=1, warn=F)) > 0) {
      line_num <- line_num + 1
      #if(grepl('^#CHROM',line)){ print(line) }
      #if(line_num==100){ break }
      if(!grepl('^##',line)){ break }
   }
   close(con)

   vcf <- read.delim(
      vcf.file, skip=line_num-1,
      check.names=F, stringsAsFactors=F,
      colClasses='character', na.strings=''
   )
   vcf$POS <- as.integer(vcf$POS)

   ## Remove '#' from header line
   colnames(vcf) <- sub('^#','',colnames(vcf))

   ## Select desired columns
   if(!is.null(fields)){
      vcf <- vcf[,fields,drop=F]
   }

   return(vcf)
}



####################################################################################################
#' Get values from INFO field
#'
#' @param v A character vector of the INFO field, with each item being a line of the INFO field
#' @param keys A character vector of the names of the INFO field values to retrieve
#'
#' @return A character vector (if 1 key was provided) or dataframe (if >1 key was provided)
#' containing the key names and corresponding values
#' @export
getInfoValues <- function(v, keys){
   #v=vcf$info
   #keys=c('SUBCL','PURPLE_CN')
   require(stringr)
   main <- function(x, key){
      #key='PURPLE_CN'
      #str_extract(x, paste0(key,'=.+?;'))
      x <- str_extract(x, paste0(key,'=.+?;'))
      x <- str_replace(x,'.*=','')
      str_replace(x,';','')
   }

   if(length(keys)==1){
      return(main(v, keys))
   }

   out <- lapply(keys, function(i){ main(v, i) })
   names(out) <- keys
   as.data.frame(out)
}
