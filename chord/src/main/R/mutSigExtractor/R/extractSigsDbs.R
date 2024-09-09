#' Extract doublet substitution signatures
#'
#' @description Will output a 1-column matrix containing: (if output=='signatures') the absolute
#' signature contributions (i.e. the number of mutations contributing to each mutational signature),
#' or (if output=='contexts') the mutation contexts, or (if output=='df') a dataframe with each
#' mutation annotated by context
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to
#' vcf.file
#' @param output Output the absolute signature contributions (default, 'signatures'), the
#' DBS contexts ('contexts'), or an annotated bed-like dataframe ('df')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#' this. If none is provided, the basename of the vcf file will be used.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#' @param ref.genome Deprecated. Argument only kept for compatibility.
#' @param verbose Print progress messages?
#'
#' @return A 1-column matrix containing the context counts or signature contributions
#' @export
extractSigsDbs <- function(
   vcf.file=NULL, df=NULL, output='contexts', sample.name=NULL,
   signature.profiles=DBS_SIGNATURE_PROFILES, ref.genome=NULL, verbose=F, ...
){

   ## Init --------------------------------
   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(vcf.file, verbose=verbose, ...)
      #df <- variantsFromVcf(vcf.file, ref.genome=ref.genome, verbose=verbose)
   }

   ## Filter variants
   df <- subsetSmnvs(df, type='dbs', verbose=verbose)

   if(verbose){ message('Initializing SNV signature output vector...') }
   context_counts <- structure(
      rep(0, nrow(DBS_TYPES)),
      names=DBS_TYPES$context
   )

   ## Main --------------------------------
   if(nrow(df)!=0){
      df$context <- paste0(df$ref,'>',df$alt)

      if(verbose){ message('Getting reverse complement contexts...') }
      df$index <- 1:nrow(df) ## Create index to maintain original row order
      df_split <- list(
         'TRUE'=df[df$context %in% DBS_TYPES$context,],
         'FALSE'=df[!(df$context %in% DBS_TYPES$context),] ## Needs reverse complementing
      )

      df_split[['FALSE']] <- within(df_split[['FALSE']],{
         context <- DBS_TYPES$context[ match(context, DBS_TYPES$context_rev_comp) ]
      })

      df <- rbind(df_split[['FALSE']],df_split[['TRUE']]); rm(df_split)
      df <- df[order(df$index),]; df$index <- NULL ## Restore original row order

      ## Check for weird nucleotides
      which_weird_nt <- sort(unique(c(
         grep('[^ACTG>]',df$substitution),
         grep('[^ACTG]',df$tri_context)
      )))

      if(length(which_weird_nt)>0){
         warning(
            length(which_weird_nt),
            ' variants containing nucleotides other than A,T,C,G were removed (rows: ',
            paste(which_weird_nt, collapse=', '), ')'
         )
         df <- df[-which_weird_nt,]
      }

      if(verbose){ message('Counting context occurrences...') }
      df$context <- factor(df$context, names(context_counts))
      tab <- table(df$context)
      context_counts[names(tab)] <- tab
   }

   ## Output --------------------------------
   if(output=='df'){
      if(verbose){ message('Returning annotated bed-like dataframe...') }
      return(df)
   }

   if(output=='contexts'){
      if(verbose){ message('Returning DBS context counts...') }
      out <- as.matrix(context_counts)

   } else if(output=='signatures'){
      if(verbose){ message('Returning DBS signature contributions...') }
      ## Least squares fitting
      out <- fitToSignatures(context_counts, signature.profiles)
      names(out) <- colnames(signature.profiles)
      out <- as.matrix(out)
   }

   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }

   return(out)
}


