#' Extract relevant variant info for extracting SNV, indel, and SV signatures.
#'
#' @param vcf.file Path to the vcf file
#' @param vcf.filter A character or character vector to specifying which variants to keep,
#' corresponding to the values in the vcf FILTER column
#' @param vcf.fields A character vector specifying the vcf columns to retrieve
#' @param keep.chroms A character vector specifying which chromosomes to keep (chromosome names
#' should be in the style of the vcf). To keep autosomal and sex chromosomes for example use:
#' keep.chroms=c(1:22,'X','Y')
#' @param chrom.name.style A string indicating the chromosome naming style. This can be for example
#' 'UCSC' ('chrX') or 'NCBI' ('X'). See the documentation for `GenomeInfoDb::seqlevelsStyle()` for
#' more details.
#' @param merge.consecutive Some vcfs report MNVs as consecutive variants. For these vcfs, such rows
#' need to be merged into one row for proper function of downstream mutSigExtractor functions.
#' @param verbose Print progress messages?
#'
#' @return A data frame with each vcf field as a column
#' @export
variantsFromVcf <- function(
   vcf.file, vcf.filter=NA, vcf.fields=c('CHROM','POS','REF','ALT','FILTER'),
   keep.chroms=NULL, chrom.name.style='UCSC',
   merge.consecutive=F, verbose=F
){
   # vcf.file='/Users/lnguyen/hpc/cuppen/shared_resources/PCAWG/final_consensus_12oct/PASS_vcfs/SNV/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv/0009b464-b376-4fbc-8a56-da538269a02f.consensus.20160830.somatic.snv_mnv_PASS.vcf.gz'

   if(verbose){ message('Reading in vcf file...') }
   vcf <- readVcfFields(vcf.file, vcf.fields)
   colnames(vcf) <- tolower(colnames(vcf))

   if(nrow(vcf)==0){
      if(verbose){ warning('VCF contains no rows. Returning empty dataframe') }
      return(data.frame())
   }

   ## Set chromosome names to the same used in the supplied ref genome
   vcf$chrom <- as.character(vcf$chrom)
   if(!is.null(chrom.name.style)){
      if(verbose){ message('Converting chrom name style to style in ref.genome...') }
      GenomeInfoDb::seqlevelsStyle(vcf$chrom)<- chrom.name.style
   }

   ## Keep certain chromosome types
   if(!is.null(keep.chroms)){
      if(verbose){ message('Keeping chromosomes as indicated in keep.chroms...') }
      ## Force chromosome name style to that in ref genome
      GenomeInfoDb::seqlevelsStyle(keep.chroms)<- chrom.name.style
      vcf <- vcf[vcf$chrom %in% keep.chroms,]
   }

   ## Filter vcf
   if(!is.na(vcf.filter)){
      if(verbose){ message('Keeping variants where FILTER is ', paste(vcf.filter,collapse=', '),' ...') }
      vcf <- vcf[vcf$filter %in% vcf.filter,]
      vcf$filter <- NULL
   }

   if(nrow(vcf)==0){
      if(verbose){ warning('After filtering, VCF contains no rows. Returning empty dataframe') }
      return(data.frame())
   }

   ## Flatten MNVs that are reported on consecutive rows
   if(merge.consecutive){
      if(verbose){ message('Merging consecutive rows...') }

      ## Find consecutive variants per chromosome
      vcf$chrom <- factor(vcf$chrom, unique(vcf$chrom))
      vcf_split <- lapply(split(vcf, vcf$chrom), function(i){
         #i=split(vcf, vcf$chrom)[[1]]
         pos_diff <- c(0,diff(i$pos))

         is_consecutive <- pos_diff==1
         is_consecutive[which(is_consecutive)-1] <- TRUE

         i$is_consecutive <- is_consecutive

         return(i)
      })

      vcf <- do.call(rbind, vcf_split)

      ## Groups rows into consecutive/non-consecutive
      rle_out <- rle(vcf$is_consecutive)
      vcf$group <- rep(
         1:length(rle_out$lengths),
         times=rle_out$lengths
      )
      vcf$group <- factor(vcf$group, unique(vcf$group))

      ## Merge consecutive variants
      vcf_split <- split(vcf, vcf$group)
      vcf <- do.call(rbind, lapply(vcf_split, function(i){
         if(!i$is_consecutive[1]){ return(i) }

         out <- i[1,]
         out$ref <- paste(i$ref,collapse='')
         out$alt <- paste(i$alt,collapse='')

         return(out)
      }))

      rownames(vcf) <- NULL
      vcf$is_consecutive <- vcf$group <- NULL
      vcf$chrom <- as.character(vcf$chrom)
   }

   return(vcf)
}
