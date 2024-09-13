#' Extract single nucleotide variant signatures
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
#' 96-trinucleotide contexts ('contexts'), or an annotated bed-like dataframe ('df')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#' this. If none is provided, the basename of the vcf file will be used.
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#' @param verbose Print progress messages?
#'
#' @return A 1-column matrix containing the context counts or signature contributions
#' @export
extractSigsSnv <- function(
   vcf.file=NULL, df=NULL, output='signatures', sample.name=NULL,
   ref.genome=DEFAULT_GENOME, signature.profiles=SBS_SIGNATURE_PROFILES_V2,
   verbose=F, ...
){

   ## Init --------------------------------
   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(vcf.file, verbose=verbose, ...)
   }

   if(verbose){ message('Initializing SNV signature output vector...') }
   context_counts <- structure(
      rep(0, length(SUBS_CONTEXTS_96)),
      names=SUBS_CONTEXTS_96
   )

   ## Filter variants
   df <- subsetSmnvs(df, type='snv', verbose=verbose)

   ## Main --------------------------------
   if(nrow(df)!=0){ ## Don't process empty vcfs
      ## Convert string to variable name
      if(is.character(ref.genome)){ ref.genome <- eval(parse(text=ref.genome)) }

      if(verbose){ message('Getting SNV trinucleotide contexts...') }
      df$substitution <- paste0(df$ref,'>',df$alt)
      df$tri_context = BSgenome::getSeq(
         x = ref.genome,
         names = df$chrom,
         start = df$pos - 1,
         end = df$pos + 1,
         as.character=T
      )

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
         if(nrow(df)==0){ warning('No variants remained after removing weird nucleotides') }
      }

      if(verbose){ message('Converting trinucleotide contexts to substitution contexts...') }
      ## Get trinucleotide contexts that don't conform to C>N or T>N
      select_opp_types <- which(!(df$substitution %in% SUBSTITUTIONS))

      ## Reverse complement for non-conforming contexts
      df[select_opp_types,'tri_context'] <- IRanges::reverse(
         chartr('ATGC', 'TACG', df[select_opp_types,'tri_context'])
      )

      ## Single nt can simply be complemented (no need to reverse)
      df[select_opp_types,'substitution'] <- chartr('ATGC', 'TACG', df[select_opp_types,'substitution'])

      ## Convert trinucleotide context to substitution context
      df$context <- paste0(
         substr(df$tri_context, 1, 1),
         '[', df$substitution, ']',
         substr(df$tri_context, 3, 3)
      )
      df$context <- factor(df$context, names(context_counts))

      ## Update context_counts
      ## Count context occurrences. Fill found contexts into context_counts
      if(verbose){ message('Counting substitution context occurrences...') }
      context_counts_new <- table(df$context)
      context_counts[names(context_counts_new)] <- context_counts_new
   }

   ## Output --------------------------------
   if(output=='df'){
      if(verbose){ message('Returning annotated bed-like dataframe...') }
      return(df)
   }

   if(output=='contexts'){
      if(verbose){ message('Returning context counts...') }
      out <- as.matrix(context_counts)

   } else if(output=='signatures'){
      if(verbose){ message('Returning absolute signature contributions...') }
      ## Least squares fitting
      out <- fitToSignatures(context_counts, signature.profiles)
      names(out) <- colnames(signature.profiles)
      out <- as.matrix(out)
   }
   if(verbose & !is.data.frame(df)){ warning("Input to extractSigsSnv() contained no variants. Returning dataframe of 0's") }

   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }

   return(out)
}

