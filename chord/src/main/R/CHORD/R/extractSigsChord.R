#' Extract mutation contexts in the format compatible with CHORD
#' 
#' @description This function is a wrapper for the 3 functions from mutSigExtractor:
#' extractSigsSnv(), extractSigsIndel(), extractSigsSv(). Some post-processing is done to produce
#' compatible input for CHORD
#' 
#' @param vcf.snv Path to the vcf file containing SNVs
#' @param vcf.indel Path to the vcf file containing indels. By default vcf.indel=vcf.snv
#' @param vcf.sv Path to the vcf file containing SVs
#' @param df.snv A dataframe containing the columns: chrom, pos, ref, alt.
#' @param df.indel A dataframe containing the columns: chrom, pos, ref, alt.
#' @param df.sv A dataframe with the columns: sv_type, sv_len. sv_type can be DEL, DUP, INV, TRA, BND.
#' @param sample.name The name of the sample as a character. Defaults to 'sample' if none is 
#' provided.
#' @param sv.caller SV vcfs are not standardized and therefore need to be parsed differently
#' depending on the caller. Currently supports 'manta' or 'gridss'.
#' @param vcf.filters A list of in the form list(snv=character(),indel=character(),sv=character()) 
#' indicated which variants to keep, corresponding to the values in the vcf FILTER column. NA can be
#' specified for each list item to ignore filtering of a vcf.
#' @param output.path If a path is specified, the output is written to this path.
#' @param ref.genome A BSgenome reference genome. Default is BSgenome.Hsapiens.UCSC.hg19. If another
#' reference genome is indicated, it will also need to be installed.
#' @param verbose Whether to print progress messages
#'
#' @return A 1-row data frame containing the mutational signature contributions
#' @export
#'
extractSigsChord <- function(
  vcf.snv=NULL, vcf.indel=vcf.snv, vcf.sv=NULL,
  df.snv=NULL, df.indel=df.snv, df.sv=NULL,
  sample.name='sample',
  vcf.filters=list(snv=NA,indel=NA,sv=NA),
  sv.caller='gridss', output.path=NULL, ref.genome=mutSigExtractor:::DEFAULT_GENOME, verbose=F
){

  if(is.null(vcf.snv) & is.null(df.snv)){
    stop('Either vcf.snv or df.snv inputs are required')
  }
  
  if(is.null(vcf.indel) & is.null(df.indel)){
    stop('Either vcf.indel or df.indel inputs are required')
  }
  
  if(is.null(vcf.sv) & is.null(df.sv)){
    stop('Either vcf.sv or df.sv inputs are required')
  }
  
  ######### Loading vcfs and/or dataframe inputs #########
  if(verbose){ message('\n#====== Loading variants from vcfs ======#') }
  
  variants <- list()
  
  if(verbose){ message('\n## SNVs') }
  if(!is.null(vcf.snv)){
    variants$snv <- mutSigExtractor::variantsFromVcf(
      vcf.snv,
      vcf.filter=vcf.filters$snv,
      verbose=verbose
    )
  } else {
    variants$snv <- df.snv
  }

  if(verbose){ message('\n## Indels') }
  if(!is.null(vcf.indel)){
    if(vcf.indel==vcf.snv){
      if(verbose){ message('vcf file is the same for both SNVs and indels. Skipping reading vcf for indels') }
      variants$indel <- variants$snv
    } else {
      variants$indel <- mutSigExtractor::variantsFromVcf(
        vcf.indel,
        vcf.filter=vcf.filters$indel, 
        verbose=verbose
      )
    }
  } else {
    if(identical(df.indel,df.snv)){
      variants$indel <- df.snv
    } else {
      variants$indel <- df.indel
    }
    
  }
  
  if(verbose){ message('\n## SVs') }
  if(!is.null(vcf.sv)){
    variants$sv <- mutSigExtractor::variantsFromVcf(vcf.sv, vcf.filter=vcf.filters$sv, vcf.fields=c('CHROM','POS','REF','ALT','FILTER','ID','INFO'), verbose=verbose)
    variants$sv <- mutSigExtractor::getContextsSv(variants$sv, sv.caller=sv.caller)
  } else {
    variants$sv <- df.sv
  }
  
  ######### Count contexts #########
  if(verbose){ message('\n#====== Counting mutation contexts ======#') }
  sigs <- list()
  
  if(verbose){ message('\n## Single base substitutions') }
  sigs$snv <- mutSigExtractor::extractSigsSnv(df=variants$snv, output='contexts', ref.genome=ref.genome, verbose=verbose)
  
  if(verbose){ message('\n## Indel contexts (types x lengths)') }
  sigs$indel <- mutSigExtractor::extractSigsIndel(df=variants$indel, method='chord', ref.genome=ref.genome, verbose=verbose)
  
  if(verbose){ message('\n## SV contexts (types x lengths)') }
  sigs$sv <- mutSigExtractor::extractSigsSv(
    df=variants$sv, sv.caller=sv.caller, output='contexts',
    sv.len.cutoffs = c(0, 10^3, 10^4, 10^5, 10^6, 10^7,Inf),
    verbose=verbose
  )
  
  ######### Export #########
  if(verbose){ message('\n#====== Exporting output =========#') }
  out <- do.call(cbind,lapply(sigs,t))
  rownames(out) <- sample.name
  
  if(is.null(output.path)){
    if(verbose){ message('output.path not specified. Directly returning output') }
    return(out)
  } else {
    if(verbose){ message('Writing tsv file') }
    write.table(out, output.path, sep='\t', quote=F)
  }
}
