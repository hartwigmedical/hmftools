#' Determine SV type and length
#'
#' @param df A dataframe containing the columns: chrom, pos, ref, alt, filter, id, info
#' @param sv.caller Only applies when mode=='sv'. At this moment supports 'manta', 'gridss', or 'pcawg'.
#' Currently there is no standard to how SVs are reported in vcfs. Therefore, output from different
#' callers will need to be parsed separately.
#' @param verbose Print progress messages?
#'
#' @return A dataframe in the same structure as a bed file with an extra column stating the context
#' of each variant
#' @export
getContextsSv <- function(df, sv.caller='gridss', return.raw=F, verbose=F){

   if(nrow(df)==0){
      return(data.frame())
   }

   #========= Callers that report SVs =========#
   if(sv.caller=='manta'){

      if(verbose){ message('Returning SV length and type...') }

      out <- getInfoValues(df$info,c('SVTYPE','SVLEN'))
      colnames(out) <- c('sv_type','sv_len')
      out$sv_len <- as.numeric(out$sv_len)
      out$sv_len[out$sv_type=='TRA'] <- NA

      return(out)
   }

   #========= Callers that report breakends =========#
   df$info <- NULL

   #--------- Determine breakend pairs ---------#
   if(verbose){ message('Identifying breakend pairs...') }
   ## Pair codes:
   ## 0: unpaired
   ## 1: 5' breakend
   ## 2: 3' breakend
   if(sv.caller=='gridss'){
      id_nchars <- nchar(df$id)

      #df$group_id <- substr(df$id, 1, id_nchars-1)
      df <- df[grepl('^gridss',df$id),]
      df$pair_id <- (function(){
         v <- rep(0,nrow(df))
         v[grepl('o$',df$id)] <- 1
         v[grepl('h$',df$id)] <- 2
         return(v)
      })()

   } else if(sv.caller=='pcawg'){
      split_ids <- strsplit(df$id,'_')
      #df$group_id <- sapply(split_ids,`[[`,1)
      df$pair_id <- as.integer(sapply(split_ids,`[[`,2))
   }

   ## Sort by group id and pair id while (mostly) maintaining original vcf order
   ## 5' breakend is forced to be 1st row
   #df$group_id <- factor(df$group_id, unique(df$group_id))
   #df <- df[order(df$group_id, df$pair_id),]

   #--------- SV type and length ---------#
   if(verbose){ message('Determining SV length and type...') }
   ## Select 5prime pairs
   df_ss <- df[df$pair_id %in% c(1,0),]

   ## Get ALT coordinates
   alt_coord <- regmatches(df_ss$alt, gregexpr('\\d|\\w+:\\d+', df_ss$alt))
   alt_coord <- as.data.frame(do.call(rbind, lapply(alt_coord, function(i){
      if(length(i) == 0){ c(NA,NA) }
      else { unlist(strsplit(i, ':')) }
   })))
   colnames(alt_coord) <- c('chrom_alt','pos_alt')
   df_ss <- cbind(df_ss,alt_coord); rm(alt_coord)
   df_ss$chrom <- gsub('chr','',df_ss$chrom)
   df_ss$chrom_alt <- gsub('chr','',df_ss$chrom_alt)

   df_ss$pos_alt <- as.numeric(as.character(df_ss$pos_alt))
   df_ss$sv_len <- df_ss$pos_alt - df_ss$pos

   ## Decision tree
   df_ss$sv_type <- with(df_ss,{
      do.call(rbind,Map(function(pair_id, chrom, chrom_alt, alt, sv_len){
         if(pair_id==0){ return('SGL') }
         if(chrom != chrom_alt){ return('TRA') }
         if(sv_len==1){ return('INS') }

         if(grepl('\\w+\\[.+\\[', alt)){ return('DEL') }
         if(grepl('\\].+\\]\\w+', alt)){ return('DUP') }

         if(grepl('\\w+\\].+\\]', alt) | grepl('\\[.+\\[\\w+', alt) ){
            return('INV')
         }

         return(NA)

      }, pair_id, chrom, chrom_alt, alt, sv_len, USE.NAMES=F))
   })

   df_ss[df_ss$sv_type  %in% c('TRA','SGL'),'sv_len'] <- NA

   if(!return.raw){
      df_ss[,c('sv_type','sv_len')]
   } else {
      df_ss
   }

}

####################################################################################################
#' Extract structural variant signatures
#'
#' @description Will return a 1-column matrix containing: (if output = 'signatures') the absolute
#' signature contributions (i.e. the number of mutations contributing to each mutational signature),
#' or (if output = 'contexts') the mutation contexts.
#'
#' To elaborate, the 6 SV signatures used are those described in this paper: https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip,
#' in Supplementary.Table.21.Signatures.v3.xlsx. These are derived from mutation contexts composed
#' of SV type/length.
#'
#' Note that the probabilities of the clustered and non-clustered rearrangements in the signature
#' profile have been combined. In other words, whether the rearrangements were
#' clustered/non-clustered were not considered.
#'
#' @param vcf.file Path to the vcf file
#' @param df A dataframe with the columns: sv_type, sv_len. sv_type can be DEL, DUP, INV, TRA, BND.
#' Note that for TRA and BND, sv_len will be ignored. Alternative input option to vcf.file
#' @param output Output the absolute signature contributions (default, 'signatures'), or the SV
#' type/length contexts ('contexts')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#' this. If none is provided, the basename of the vcf file will be used.
#' @param half.tra.counts Divide translocation counts by 2?
#' @param sv.caller Can be 'manta', 'gridss', or 'pcawg'
#' @param sv.len.cutoffs SV length cutoff intervals as a numeric vector.
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are  the mutational signatures.
#' @param verbose Print progress messages?
#'
#' @return A 1-column matrix containing the context counts or signature contributions
#' @export
extractSigsSv <- function(
   vcf.file=NULL, df=NULL, output='signatures', sample.name=NULL,
   sv.caller='gridss', half.tra.counts=F,
   sv.len.cutoffs=if(output=='signatures'){ c(10^c(3:7), Inf) } else { c(0, 10^c(3:7), Inf) },
   signature.profiles=SV_SIGNATURE_PROFILES,
   verbose=F, ...
){
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(
         vcf.file, vcf.fields=c('CHROM','POS','REF','ALT','FILTER','ID','INFO'),
         verbose=verbose, ...
      )
      df <- getContextsSv(df, sv.caller=sv.caller, verbose=verbose)

   } else if(!is.null(df)){
      #colnames(df) <- c('sv_type','sv_len')
      half.tra.counts <- F ## If providing dataframe as input default to 'manta'.
   } else {
      stop('Please specify either vcf.file or df as input')
   }

   if(verbose){ message('Creating SV type/length lookup table...') }
   sv_types <- c('DEL','DUP','INV') ## INS ignored. TRA/BND dealt with in a later step

   sv_contexts <- data.frame(
      sv_type = rep( sv_types, each = length(sv.len.cutoffs)-1 ),
      lower_cutoff = rep( sv.len.cutoffs[-length(sv.len.cutoffs)], length(sv_types) ),
      upper_cutoff = rep( sv.len.cutoffs[-1], length(sv_types) ),

      stringsAsFactors = F
   )

   sv_contexts$name <- with(sv_contexts,{
      v <- paste(
         sv_type,
         formatC(lower_cutoff, format = 'e', digits = 0),
         formatC(upper_cutoff, format = 'e', digits = 0),
         'bp',sep='_'
      )
      gsub('[+]','',v)
   })

   ## Deal with empty vcfs
   if(nrow(df)==0){
      context_counts <- rep(0, nrow(sv_contexts)+1)
   } else {
      if(verbose){ message('Counting DEL, DUP, and INV context occurrences...') }

      context_counts <- unlist(lapply(1:nrow(sv_contexts), function(i){
         row <- sv_contexts[i,]
         variants_ss <- df[
            df$sv_type == row$sv_type
            & df$sv_len >= row$lower_cutoff
            & df$sv_len < row$upper_cutoff
         ,]

         nrow(variants_ss)
      }))

      ## Count context occurrences for translocations
      if(verbose){ message('Counting TRA occurrences...') }
      translocation_counts <- nrow(df[df$sv_type == 'BND' | df$sv_type == 'TRA',])

      if(sv.caller=='manta' & half.tra.counts){ ## manta reports translocations twice (origin/destination)
         if(verbose){ message('Halving TRA counts...') }
         translocation_counts <- translocation_counts/2
      }

      context_counts <- c(context_counts,translocation_counts)
   }

   ## Assign context names
   names(context_counts) <- c(sv_contexts$name,'TRA')

   if(output == 'contexts'){
      if(verbose){ message('Returning SV contexts...') }
      out <- as.matrix(context_counts)

   } else if(output == 'signatures'){
      if(verbose){ message('Returning SV signatures...') }

      if(nrow(signature.profiles) != length(context_counts)){
         stop('The number of contexts in the signature profile matrix != the number of contexts in the context count vector.\n
              Check that the provided cutoffs sv.len.cutoffs also exists in signature.profiles')
      }

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
