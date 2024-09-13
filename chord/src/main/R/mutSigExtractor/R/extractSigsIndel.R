#' Extract indel signatures
#'
#' @param description Will return a 1-column matrix containing the absolute indel signature
#' contributions (i.e. the number of mutations contributing to each mutational signature).
#'
#' Two sets of indel contexts can be used: CHORD and PCAWG.
#'
#' For CHORD indel contexts, signatures used are insertions/deletions within repeat regions
#' (ins.rep, del.rep), insertions/deletions with flanking microhomology (ins.mh, del.mh), and
#' insertions/deletions which don't fall under the previous 2 categories (ins.none, del.none). Each
#' category is further stratified by the length of the indel.
#'
#' PCAWG indel contexts are described at: https://cancer.sanger.ac.uk/cosmic/signatures/ID/index.tt
#'
#' @param method Can be 'CHORD' or 'PCAWG'. Indicates the indel context type to extract.
#' @param vcf.file Path to the vcf file
#' @param df A dataframe containing the columns: chrom, pos, ref, alt. Alternative input option to vcf.file
#' @param output Output the absolute signature contributions (default, 'signatures'), indel
#' contexts ('contexts'), or an annotated bed-like dataframe ('df')
#' @param sample.name If a character is provided, the header for the output matrix will be named to
#'   this. If none is provided, the basename of the vcf file will be used.
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param indel.len.cap Specifies the max indel sequence length to consider when counting 'repeat'
#'   and 'none' contexts. Counts of longer indels will simply be binned to the counts of contexts at
#'   the max indel sequence length.
#' @param n.bases.mh.cap Specifies the max bases in microhomology to consider when counting repeat
#'   and microhomology contexts. Counts of longer indels will simply be binned to the counts of
#'   contexts at the max indel sequence length.
#' @param get.other.indel.allele Only applies when mode=='indel' For indels, some vcfs only report
#' the sequence of one allele (REF for deletions and ALT for insertions). If TRUE, the unreported
#' allele will be retrieved from the genome: a 5' base relative to the indel sequence. This base
#' will also be added to the indel sequence and the POS will be adjusted accordingly (POS=POS-1).
#' @param keep.indel.types A character vector of indel types to keep. Defaults to 'del' and 'ins' to
#' filter out MNVs (variants where REF and ALT length >= 2). MNV names are: 'mnv_neutral'
#' (REF lenth == ALT length), 'mnv_del' (REF length > ALT length), or 'mnv_ins' (REF length < ALT length).
#' @param verbose Print progress messages?
#' @param ... Other arguments that can be passed to variantsFromVcf()
#'
#' @return A 1-column matrix containing the context counts or signature contributions
#' @export
#'
extractSigsIndel <- function(..., method='CHORD'){
   method <- tolower(method)
   if(method=='chord'){
      extractSigsIndelChord(...)
   } else if(method=='pcawg'){
      extractSigsIndelPcawg(...)
   } else {
      stop("`method` must be 'CHORD' or 'PCAWG'")
   }
}

####################################################################################################
#' @rdname extractSigsIndel
extractSigsIndelPcawg <- function(
   vcf.file=NULL, df=NULL, output='contexts', sample.name=NULL, ref.genome=DEFAULT_GENOME,
   signature.profiles=INDEL_SIGNATURE_PROFILES, verbose=F, ...
){

   ## Init --------------------------------
   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(vcf.file, verbose=verbose, ...)
   }

   ## Filter variants
   df <- subsetSmnvs(df, type='indel', verbose=verbose)

   context_counts <- structure(
      rep(0, length(INDEL_CONTEXTS)),
      names=INDEL_CONTEXTS
   )

   ## Main ################################
   if(nrow(df)!=0){
      df$index <- 1:nrow(df)
      df$mut_len <- nchar(df$alt) - nchar(df$ref)

      ## --------------------------------
      if(verbose){ message('Identifying the deleted or inserted bases...') }
      df$mut_type <- ifelse(df$mut_len>=1, 'ins', 'del')
      df_split <- list( ## Don't use split to prevent errors when nrow(df)==0
         ins=df[df$mut_type=='ins',],
         del=df[df$mut_type=='del',]
      )
      df_split <- lapply(names(df_split), function(i){
         df_ss <- df_split[[i]]

         ## Deleted or inserted base(s)
         ## EXcludes the REF base for ins
         ## INcludes the REF base fo dels
         df_ss$mut_seq <-
            if(i=='del'){
               substring(df_ss$ref, first=2)
            } else {
               substring(df_ss$alt, first=2)
            }

         return(df_ss)
      })

      df <- do.call(rbind, df_split); rm(df_split); rownames(df) <- NULL
      df <- df[order(df$index),] ## Restore original row order

      ## --------------------------------
      if(verbose){ message("Getting the 3' and 5' flanking sequences...") }

      ## Convert string to variable name
      if(is.character(ref.genome)){ ref.genome <- eval(parse(text=ref.genome)) }

      ## EXcluding the REF base
      ## INcludes deleted sequence
      flank_size_r <- abs(df$mut_len) * 7

      df$flank_seq_r <- BSgenome::getSeq(
         x=ref.genome,
         names=df$chrom,
         start=df$pos+1,
         end  =df$pos+1+flank_size_r,
         as.character=T
      )

      ## INcluding the REF base; only relevant for determining del mh
      df$flank_seq_l <- BSgenome::getSeq(
         x=ref.genome,
         names=df$chrom,
         start=df$pos-abs(df$mut_len),
         end  =df$pos,
         as.character=T
      )

      ## --------------------------------
      if(verbose){ message('Determining repeat unit lengths...') }
      countRepeatUnits <- function(mut_seq, flank_seq){
         #mut_seq=c('AAG','TCTCTC','T','AG','A')
         #flank_seq=c('AAGAAGTAAG','TCTCTCTCTCTCTCTCTCAA','GTGGA','AGAGTTAG','AAAAAAAA')

         if(length(mut_seq) != length(flank_seq)){
            stop('length(mut_seq) must equal length(flank_seq)')
         }

         ## For each mut seq replace the inserted/deleted sequence in the flanking sequence with Zs
         flank_seq_z <- unlist(Map(gsub, pattern=mut_seq, replacement='Z', x=flank_seq, USE.NAMES=F))

         ## Remove all bases after the Zs and count how many bases are left
         flank_seq_z <- gsub("[^Z].*", "", flank_seq_z)
         nchar(flank_seq_z)
      }

      df$rep_len <- countRepeatUnits(df$mut_seq, df$flank_seq_r)

      ## --------------------------------
      if(verbose){ message('Determining homology lengths...') }
      countHomBases <- function(mut_seq, flank_seq){

         #mut_seq=c('A','A','TC')
         #flank_seq=c('AATTTTTT','AAAAAAAA','TAGCGGC')

         if(length(mut_seq) != length(flank_seq)){
            stop('length(mut_seq) must equal length(flank_seq)')
         }

         mut_seq_split <- strsplit(mut_seq, '')
         flank_seq_split <- strsplit(flank_seq, '')

         suppressWarnings({ ## Suppress warning when all chunks match
            out <- unlist(lapply(1:length(mut_seq), function(i){
               #i=1
               i_mut_seq <- mut_seq_split[[i]]

               i_flank_seq <- flank_seq_split[[i]]
               i_flank_seq <- i_flank_seq[1:length(i_mut_seq)]

               i_mut_seq==i_flank_seq

               first_no_match <- min(which(i_mut_seq!=i_flank_seq))
               first_no_match - 1
            }), use.names=T)
         })

         out[out==Inf] <- nchar(mut_seq)[out==Inf] ## When all bases match, replace with mut seq length

         return(out)
      }

      df$hom_len <- (function(){
         #countHomBases(mut_seq='A', flank_seq='AATTTTTT')
         #countHomBases(mut_seq=IRanges::reverse('TATC'), flank_seq=IRanges::reverse('ACCCGTC'))

         flank_seq_r_excl <- substring(df$flank_seq_r, nchar(df$mut_seq)+1)
         hom_len_r <- countHomBases(
            df$mut_seq,
            flank_seq_r_excl
         )
         hom_len_l <- countHomBases(
            IRanges::reverse(df$mut_seq),
            IRanges::reverse(df$flank_seq_l))

         pmax(hom_len_r, hom_len_l)
      })()

      ## --------------------------------
      if(verbose){ message('Determining indel subtype...') }
      #df$flank_seq_r <- NULL
      #df$flank_seq_l <- NULL

      ## Process 1bp and >=2bp indels separately
      df_split <- list(
         small = df[abs(df$mut_len)==1,],
         large = df[abs(df$mut_len)>=2,]
      )

      ## Get homopolymer type
      df_split$small <- within(df_split$small,{
         mut_subtype <- mut_seq
         mut_subtype[mut_subtype=='A'] <- 'T'
         mut_subtype[mut_subtype=='G'] <- 'C'
      })

      ## Initialize all large indels as of type repeat. Then override to type homology when applicable
      df_split$large$mut_subtype <- rep('rep',nrow(df_split$large))
      df_split$large <- within(df_split$large,{
         mut_subtype[
            mut_type=='del'
            & abs(mut_len)>=2
            & rep_len==1
            & hom_len>=1
         ] <- 'mh'
      })

      df <- do.call(rbind, df_split); rm(df_split); rownames(df) <- NULL
      df <- df[order(df$index),] ## Restore original row order

      ## --------------------------------
      if(verbose){ message('Calculating capped lengths...') }
      df <- within(df,{

         mut_len <- abs(mut_len)
         mut_len[mut_len>=5] <- '5+'

         rep_len[mut_type=='del' & rep_len>=6] <- '6+'
         rep_len[mut_type=='ins' & rep_len>=5] <- '5+'

         hom_len[hom_len>=5] <- '5+'

         mut_subtype_len <- rep_len
         mut_subtype_len[mut_subtype=='mh'] <- hom_len[mut_subtype=='mh']
      })

      if(verbose){ message('Determining indel context...') }
      if(nrow(df)==0){
         df$context <- character()
      } else {
         df$context <- with(df,{
            paste0(
               mut_type,'.',mut_len,'.',
               mut_subtype,'.',mut_subtype_len
            )
         })
      }

      if(verbose){ message('Counting context occurrences...') }
      df$context <- factor(df$context, names(context_counts))
      context_counts_new <- table(df$context)
      context_counts[names(context_counts_new)] <- context_counts_new
      context_counts <- context_counts[INDEL_CONTEXTS]
   }

   ## Output --------------------------------
   if(output == 'df'){
      if(verbose){ message('Returning annotated bed-like dataframe...') }
      return(df)
   }

   if(output == 'contexts'){
      if(verbose){ message('Returning context counts...') }
      out <- as.matrix(context_counts)

   } else if(output == 'signatures'){
      if(verbose){ message('Returning absolute signature contributions...') }
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

####################################################################################################
#' @rdname extractSigsIndel
extractSigsIndelChord <- function(
   vcf.file=NULL, df=NULL, sample.name=NULL, ref.genome=DEFAULT_GENOME, output='contexts',
   indel.len.cap=5, n.bases.mh.cap=5, get.other.indel.allele=F, keep.indel.types=c('del','ins'),
   verbose=F, ...
){

   ## Init --------------------------------
   if(verbose){ message('Loading variants...') }
   if(!is.null(vcf.file)){
      df <- variantsFromVcf(vcf.file, verbose=verbose, ...)
   }

   ## Filter variants
   df <- subsetSmnvs(df, type='indel', verbose=verbose)

   if(verbose){ message('Initializing indel signature output vector...') }
   indel_sig_names <- c(
      paste0('del.rep.len.', 1:indel.len.cap),
      paste0('ins.rep.len.', 1:indel.len.cap),
      paste0('del.mh.bimh.', 1:n.bases.mh.cap),
      paste0('ins.mh.bimh.', 1:n.bases.mh.cap),
      paste0('del.none.len.', 1:indel.len.cap),
      paste0('ins.none.len.', 1:indel.len.cap)
   )
   indel_sigs <- structure(rep(0,length(indel_sig_names)), names=indel_sig_names)

   ## Main ################################
   if(nrow(df)!=0){ ## Don't process empty vcfs

      df$ref_len <- nchar(df$ref)
      df$alt_len <- nchar(df$alt)

      ## Determine indel type
      df$indel_type <- with(df,{
         unlist(Map(function(ref_len, alt_len){
            if(ref_len >= 2 & alt_len >= 2){
               if(ref_len == alt_len){ 'mnv_neutral' }
               else if(ref_len > alt_len){ 'mnv_del' }
               else if(ref_len < alt_len){ 'mnv_ins' }
            } else if(ref_len > alt_len){
               'del'
            } else if(ref_len < alt_len){
               'ins'
            }
         },ref_len, alt_len, USE.NAMES=F))
      })

      ## Get REF or ALT allele for vcfs that only report one or the other --------------------------------
      ## Convert string to variable name
      if(is.character(ref.genome)){ ref.genome <- eval(parse(text=ref.genome)) }

      if(get.other.indel.allele==T){
         if(verbose){ message('Retrieving other indel allele...') }
         df_split <- lapply(
            list(del_type=c('del','mnv_del'),ins_type=c('ins','mnv_ins'),mnv_neutral='mnv_neutral'),
            function(i){ df[df$indel_type %in% i, ] }
         )

         if(nrow(df_split$del_type)!=0){
            ## Deletions
            ## ref:   'AGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
            ## alt:  'TAGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
            ## nchar('AGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC') = nchar(ref) = 46
            ## getSeq(x=ref.genome, names='chr4',start=84726292-1,84726292+46-1)
            ##       'TAGAACTACCATATGACCCAGCAGTCCCATTCTGGGTATATATCCAC'
            ## 5' base relative to ref ->
            ##    alt column
            ##    ref sequence
            df_split$del_type$alt <- with(df_split$del_type, {
               BSgenome::getSeq(
                  x=ref.genome,
                  names=chrom, start=pos-1,end=pos-1,
                  as.character=T
               )
            })
            df_split$del_type$ref <- with(df_split$del_type, { paste0(alt, ref) })
            df_split$del_type$pos <- df_split$del_type$pos-1
         }

         if(nrow(df_split$ins_type)!=0){
            ## Insertions
            ## ref:  ''
            ## alt:  'AGAGAGAGAGACAGAA'
            ## nchar('AGAGAGAGAGACAGAA') = nchar(alt) = 16
            ## getSeq(x=ref.genome, names='chr12',start=6902128-1,6902128+16-1)
            ##      'GAGAGAGAGAGACAGAA'
            ## 5' base relative to alt ->
            ##    ref column
            ##    alt sequence
            ## Substract 1 from pos
            df_split$ins_type$ref <- with(df_split$ins_type, {
               BSgenome::getSeq(
                  x=ref.genome,
                  names=chrom, start=pos-1,end=pos-1,
                  as.character=T
               )
            })
            df_split$ins_type$alt <- with(df_split$ins_type, { paste0(ref,alt) })
            df_split$ins_type$pos <- df_split$ins_type$pos-1
         }

         ## Unsplit df
         df <- do.call(rbind, df_split)
         rownames(df) <- NULL

         ## Recalculate ref/alt length
         df$ref_len <- nchar(df$ref)
         df$alt_len <- nchar(df$alt)
      }

      ## --------------------------------
      if(verbose){ message('Determining indel length and sequence...') }
      ## Determine indel length
      df$indel_len <- abs(df$alt_len-df$ref_len)

      ## Determine indel seq
      df$indel_seq <- with(df,{
         unlist(Map(function(ref,alt,indel_type,indel_len){
            indel_start_pos <- 2
            if(indel_type %in% c('del','mnv_del')){ ## dels
               substring(ref, indel_start_pos, indel_start_pos+indel_len-1)
            } else if(indel_type %in% c('ins','mnv_ins')){ ## ins
               substring(alt, indel_start_pos, indel_start_pos+indel_len-1)
            } else {
               NA
            }
         },ref,alt,indel_type,indel_len, USE.NAMES=F))
      })

      df <- df[df$indel_type %in% keep.indel.types,]
      if(nrow(df)==0){ warning('No variants remain after filtering for: ',paste(keep.indel.types, collapse=', ')) }
   }

   if(nrow(df)!=0){
      ## Subset for relevant columns --------------------------------

      df <- df[,c('chrom','pos','ref','alt','indel_len','indel_type','indel_seq')]

      ## Pre-calculations for repeat and microhomology contexts --------------------------------
      if(verbose){ message('Determining the start/end positions for the left/right flanks of each indel...') }
      ## Getting flank start/end positions first, then retrieving the sequences using getSeq with a
      ## vectorized input improves speed significantly.
      ## Cap n.indel.lengths.r to 3 (length of 3' (right-hand side) sequence to retrieve),
      ## since the max n_copies_along_flank condition used below caps at >=2.
      ## This improves speed significantly
      flanks_start_end <- with(df,{
         do.call(rbind, Map(
            f = indelSeqFlanksStartEnd,
            chrom, pos, indel_len, indel_type,
            n.indel.lengths.r=3
         ))
      })

      if(verbose){ message('Retrieving flanking sequences...') }
      l_flank <- BSgenome::getSeq(
         x = ref.genome,
         names = df$chrom,
         start = flanks_start_end[,'l_start'],
         end = flanks_start_end[,'l_end'],
         as.character = T
      )

      r_flank <- BSgenome::getSeq(
         x = ref.genome,
         names = df$chrom,
         start = flanks_start_end[,'r_start'],
         end = flanks_start_end[,'r_end'],
         as.character = T
      )

      ## Repeat contexts --------------------------------
      if(verbose){ message("Calculating the number of copies of the indel sequence are present in the 3' flanking sequence...") }
      df$n_copies_along_flank <- unlist(Map(nCopiesAlongFlank, df$indel_seq, r_flank, USE.NAMES=F))

      ## Microhomology contexts --------------------------------
      if(verbose){ message("Calculating the (max) number of bases that are homologous to the 5'/3' flanking sequence...") }
      df$n_bases_mh <- unlist(Map(function(indel_seq, l_flank, r_flank){
         mh_l <- nBasesMH(IRanges::reverse(indel_seq), IRanges::reverse(l_flank))
         mh_r <- nBasesMH(indel_seq, r_flank)

         max(mh_l,mh_r)

      }, df$indel_seq, l_flank, r_flank, USE.NAMES=F))

      ## Assign repeat, microhomology, or no context --------------------------------
      if(verbose){ message('Determining indel contexts...') }
      df$context <- unlist(Map(function(n_copies_along_flank, n_bases_mh, indel_len){
         if (n_copies_along_flank >= 2){
            if(indel_len < 50){ context <-'rep' }
            else { context <- 'mh' }

         } else if(n_copies_along_flank >= 1 && n_bases_mh >= 2) {
            context <- 'mh'
         } else if(n_copies_along_flank >= 1 && n_bases_mh >= 1 && indel_len > 3 ) {
            context <- 'mh'

         } else {
            context <- 'none'
         }

         return(context)

      }, df$n_copies_along_flank, df$n_bases_mh, df$indel_len))

      if(output=='df'){ return(df) }

      ## Gather components for counting final contexts/signatures --------------------------------
      if(verbose){ message('Counting indel context occurrences...') }
      ## Slightly redundant (could have just assigned components to a dataframe). But easier to debug
      sig_parts <- data.frame(
         indel_type = df$indel_type,
         context=df$context,
         indel_len = df$indel_len,
         n_copies_along_flank=df$n_copies_along_flank,
         n_bases_mh=df$n_bases_mh
      )

      ## Bin values larger than cap into one bin for indel_len and n_bases_mh
      sig_parts <- within(sig_parts,{
         indel_len[indel_len >= indel.len.cap] <- indel.len.cap
         n_bases_mh[n_bases_mh >= n.bases.mh.cap] <- n.bases.mh.cap
      })

      ## Count occurrences of each signature
      sig_occurrences <- table(with(sig_parts,{
         unlist(Map(function(indel_type,context,indel_len,n_copies_along_flank,n_bases_mh){
            if(context == 'mh'){
               paste(indel_type, context, 'bimh', n_bases_mh, sep = '.')
            } else {
               paste(indel_type, context, 'len', indel_len, sep = '.')
            }
         },
         indel_type, context, indel_len, n_copies_along_flank, n_bases_mh))
      }))

      ## Fill in indel signature matrix that was initiated at the start of the function
      indel_sigs[names(sig_occurrences)] <- sig_occurrences
   }

   ## Output --------------------------------
   if(verbose){ message('Returning indel context counts...') }
   if(verbose & !is.data.frame(df)){ warning("Input contained no variants. Returning dataframe of 0's") }
   out <- matrix(indel_sigs, ncol = 1)
   rownames(out) <- names(indel_sigs)

   colnames(out) <-
      if(is.null(sample.name)){
         if(!is.null(vcf.file)){ basename(vcf.file) }
         else { 'unknown_sample' }
      } else {
         sample.name
      }

   return(out)
}

