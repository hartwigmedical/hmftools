#' Assign a mutational signature to each mutation in a vcf or bed-like dataframe
#'
#' @description This function performs the following steps to assign a signature to each mutation
#' in a vcf or bed-like dataframe:
#' - Extract contexts and count occurrences
#' - Calculate signature contributions by fitting contexts counts to reference signature profiles
#' - Multiply each reference signature profile by the signature contribution vector. This gives
#' the probability of each context to be assigned to a signature
#' - Assign each context a signature based on the maximum probability signature
#' - Assign each mutation a signature based on its context
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to
#' vcf.file
#' @param mode Can be 'snv','sbs','indel' or 'dbs'. `mode` defines which `extractSigs*()` function
#' to use. Additionally, if `signature.profiles` is unspecified, the correct signature profile
#' matrix will be automatically chosen
#' @param output Can be 'df' (a dataframe with chrom, pos, ref, alt, context, assigned_sig, sig_prob),
#' 'df.compact' (a dataframe with context, assigned_sig, sig_prob), 'details' or 'vector' (a vector whose
#' names are assigned_sig, and values are sig_prob).
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are the mutational signatures.
#' @param args.extract.sigs A list of args that can be passed to the `extractSigs*()` functions
#' @param fit.method Can be 'lsq' or 'strict'. Method for fitting context counts to signature
#' profiles. See documentation for `fitToSignatures()` for more details
#' @param args.fit A list of args that can be passed to the `fitToSignatures*()` functions
#' @param min.mut.load Samples with fewer mutations than this value will have no mutations assigned
#' to any signatures (all mutations assigned as NA). If 'auto', one of the following values will be
#' selected based on the mut type (i.e. `mode`): snv=500, dbs=10, indel=50
#' @param min.sig.abs.contrib Mutations attributed to signatures with absolute contribution lower
#' than this value will be assigned NA instead. If 'auto', one of the following values will be
#' selected based on the mut type (i.e. `mode`): snv=100, dbs=5, indel=10
#' @param min.sig.rel.contrib Signatures with low relative contribution (less than this value, e.g.
#' 0.05 for 5 percent) of will be excluded from the assignment.
#' @param verbose Print progress messages?
#'
#' @return See `output`
#' @export
#'
assignSigPerMut <- function(
   vcf.file=NULL, df=NULL, mode='snv', output='df',
   signature.profiles=NULL,
   args.extract.sigs=list(vcf.filter='PASS'),
   fit.method='strict', args.fit=NULL,
   min.mut.load='auto', min.sig.abs.contrib='auto', min.sig.rel.contrib=0.05,
   verbose=F
){

   ## Checks --------------------------------
   mode <- tolower(mode)
   if(!(mode %in% c('snv','sbs','indel','id','dbs'))){
      stop("`mode` must be one of the following: 'snv','sbs','indel','id,'dbs'")
   }
   if(mode=='id'){ mode <- 'indel' }
   if(mode=='sbs'){ mode <- 'snv' }

   if(!(fit.method %in% c('lsq','strict'))){
      stop("`fit.method` must be one of the following: 'lsq','strict'")
   }

   if(!(output %in% c('df','df.compact','details','vector'))){
      stop("`output` must be one of the following: 'df','df.compact','details','vector'")
   }

   ## Default thresholds --------------------------------
   if(!is.null(min.mut.load) && min.mut.load=='auto'){
      min.mut.load <- switch(
         mode,
         snv=500, dbs=10, indel=50
      )
   }

   if(!is.null(min.sig.abs.contrib) && min.sig.abs.contrib=='auto'){
      min.sig.abs.contrib <- switch(
         mode,
         snv=100, dbs=5, indel=10
      )
   }

   ## --------------------------------
   if(verbose){ message('Extracting ',mode,' contexts...') }
   f_extract_sigs <- switch(
      mode,
      snv=extractSigsSnv,
      dbs=extractSigsDbs,
      indel=extractSigsIndel
   )

   ## Oligatory args
   args_extract_sigs <- list(vcf.file=NULL, df=NULL, output='df')
   if(mode=='indel'){ args_extract_sigs <- c(args_extract_sigs, method='PCAWG') }

   ## Optional args
   args_extract_sigs <- c(args_extract_sigs, args.extract.sigs)

   ## Rm duplicate args
   args_extract_sigs <- args_extract_sigs[ !duplicated(names(args_extract_sigs)) ]

   if(!is.null(vcf.file)){
      args_extract_sigs$vcf.file <- vcf.file
   } else if(!is.null(df)){
      args_extract_sigs$df <- df
   } else {
      stop('Input must be specified to `vcf.file` or `df`')
   }

   df <- do.call(f_extract_sigs, args_extract_sigs)
   rm(f_extract_sigs, args_extract_sigs)

   ## --------------------------------
   if(!is.null(signature.profiles)){
      sig_profiles <- signature.profiles
   } else {
      sig_profiles <- switch(
         mode,
         snv=SBS_SIGNATURE_PROFILES_V3,
         dbs=DBS_SIGNATURE_PROFILES,
         indel=INDEL_SIGNATURE_PROFILES
      )
   }

   if(nrow(df)!=0){

      ## Select only relevant columns
      df <- df[,c('chrom','pos','ref','alt','context')]

      if(verbose){ message('Counting mutation contexts...') }
      context_counts_pre <- table(df$context)
      context_counts <- structure(
         rep(0, length(levels(df$context))),
         names=levels(df$context)
      )
      context_counts[names(context_counts_pre)] <- context_counts_pre

      ## --------------------------------
      if(verbose){ message('Fitting context counts to signature profiles...') }

      ## Force context names and order to be the same in sig profile matrix
      if(!all(names(context_counts) %in% rownames(sig_profiles))){
         stop('Extracted contexts do not match `rownames(signature.profiles)`. Maybe the wrong `mode` was selected?')
      }
      sig_profiles <- sig_profiles[names(context_counts),,drop=F]

      f_fit <- switch(
         fit.method,
         lsq=fitToSignatures,
         strict=fitToSignaturesStrict
      )

      args_fit <- c(
         list(mut.context.counts=context_counts, signature.profiles=sig_profiles),
         args.fit
      )
      args_fit <- args_fit[ !duplicated(names(args_fit)) ]

      sig_contrib <- do.call(f_fit, args_fit)
      sig_contrib[is.na(sig_contrib)] <- 0 ## Signatures with all 0 probs leads to NA contributions

      rm(f_fit, args_fit)

      ## --------------------------------
      if(!is.null(min.sig.rel.contrib)){

         if(verbose){ message('Remove signatures with <',min.sig.rel.contrib,' rel. contribution...') }

         sig_contrib_rel <- sig_contrib/sum(sig_contrib)
         sig_contrib_rel[is.na(sig_contrib_rel)] <- 0

         sig_blacklist <- names(sig_contrib_rel)[ sig_contrib_rel < min.sig.rel.contrib ]
         sig_contrib[sig_blacklist] <- 0
      }

      # ## Scale sig contrib to have original total mutational load from the fitting
      # sig_contrib <- sig_contrib * (sum(context_counts) / sum(sig_contrib))

      ## --------------------------------
      if(sum(context_counts)>=min.mut.load){
         if(verbose){ message('Calculating signature probabilities per mutation...') }

         ## Adjust signature context probabilities based on signature contributions in sample
         sig_profiles_sample <- apply(sig_profiles, 1, function(i){ i * sig_contrib })

         ## Convert to context x signature matrix
         ## Avoid converting to vector when there is only one signature
         sig_profiles_sample <- matrix(
            sig_profiles_sample,
            ncol=ncol(sig_profiles), byrow=T,
            dimnames=list(rownames(sig_profiles), colnames(sig_profiles))
         )

         sig_profiles_sample <- sig_profiles_sample / rowSums(sig_profiles_sample)
         sig_profiles_sample[is.na(sig_profiles_sample)] <- 0

         ## Get max probability signature per context
         context_sig_assignment <- data.frame(
            assigned_sig = colnames(sig_profiles_sample)[ max.col(sig_profiles_sample) ],
            sig_prob = apply(sig_profiles_sample,1,max),
            row.names=rownames(sig_profiles_sample)
         )

         unassigned_contexts <- rownames(context_sig_assignment)[context_sig_assignment$sig_prob==0]
         if(length(unassigned_contexts)>0){
            warning(length(unassigned_contexts),' contexts could not be unassigned to a siganture')
            context_sig_assignment$assigned_sig[context_sig_assignment$sig_prob==0] <- NA
         }

         ## Assign each mutation a signature based on its context
         df$assigned_sig <- context_sig_assignment$assigned_sig[df$context]
         df$sig_prob <- context_sig_assignment$sig_prob[df$context]

         if(!is.null(min.sig.abs.contrib)){

            if(verbose){ message('Unassigning signatures with <',min.sig.abs.contrib,' abs. contribution...') }
            low_abs_contrib_sigs <- names(sig_contrib)[ sig_contrib<min.sig.abs.contrib ]
            df$sig_prob[df$assigned_sig %in% low_abs_contrib_sigs] <- NA
            df$assigned_sig[df$assigned_sig %in% low_abs_contrib_sigs] <- NA
         }

      } else {

         if(verbose){ message('Sample had <',min.mut.load,' mutations. No signatures were assigned') }
         df$assigned_sig <- NA
         df$sig_prob <- NA
      }

      ## --------------------------------
      if(verbose){ message('Calculating signature contributions from assignment...') }
      sig_contrib_assigned <- structure(
         rep(0,ncol(sig_profiles)),
         names=colnames(sig_profiles)
      )
      tab <- table(df$assigned_sig)
      sig_contrib_assigned[names(tab)] <- tab
      rm(tab)

   } else {
      if(verbose){ message('Input contains no variants, or all variants were filtered out') }
      df <- data.frame(
         chrom=character(),
         pos=character(),
         ref=character(),
         alt=character(),
         context=character(),
         assigned_sig=character(),
         sig_prob=character()
      )

      sig_contrib <- structure(
         rep(0,ncol(sig_profiles)),
         names=colnames(sig_profiles)
      )

      sig_contrib_assigned <- sig_contrib
   }

   ## --------------------------------
   if(verbose){ message('Returning output...') }

   if(output=='df'){ return(df) }

   if(output=='df.compact'){ return( df[,c('context','assigned_sig','sig_prob')] ) }

   if(output=='details'){
      sig_contrib_summ <- cbind(fit=as.vector(sig_contrib), assigned=as.vector(sig_contrib_assigned))
      rownames(sig_contrib_summ) <- colnames(signature.profiles)
      return(
         list(
            muts=df,
            contrib=sig_contrib_summ
         )
      )

   }

   structure(df$sig_prob, names=df$assigned_sig)
}
