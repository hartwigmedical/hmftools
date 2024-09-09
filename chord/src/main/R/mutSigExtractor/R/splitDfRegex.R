#' Split data frame or matrix using regex on colnames
#'
#' @param df A data frame of matrix
#' @param patterns A vector/list of regular expressions. If arguments are named, or items in the
#' vector/list are named, the items in the output list will also be named accordingly
#' @param return.remainder Return columns which don't match the given regular expressions?
#'
#' @return A list of data frames or matrices
#' @export
#'
#' @examples
#' ## Will return a list containing the corresponding names
#' splitDfRegex(df = contexts, patterns = list(snv = '\\[', indel = '[.]', sv = '[A-Z]{3}'))
#'
#' ## The same output list as above, but without names
#' splitDfRegex(df = contexts, patterns = list('\\[', '[.]', '[A-Z]{3}'))

splitDfRegex <- function(df, patterns, return.remainder=T, names=NULL){
   patterns <- as.list(patterns)

   df <- as.data.frame(df)

   indexes <- lapply(patterns, function(i){ grep(i, colnames(df)) })
   l_split <- lapply(indexes, function(i){ df[i] })

   if(return.remainder){
      cols <- seq(1:ncol(df))
      remainder_indexes <- cols[!(cols %in% unlist(indexes))]

      if(length(remainder_indexes) != 0){
         l_split[['remainder']] <- df[remainder_indexes]
      }
   }

   if(!is.null(names)){
      names(l_split) <- names
   }

   return(l_split)
}



