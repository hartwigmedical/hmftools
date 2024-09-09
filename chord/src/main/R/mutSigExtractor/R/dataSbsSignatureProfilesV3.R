#' SBS signature profiles (new; v3, May 2019)
#'
#' A matrix
#' rows: 96-trinucleotide context
#' cols: single base substitution (SBS) signatures
#'
#' Source: https://cancer.sanger.ac.uk/sigs-assets-20/COSMIC_Mutational_Signatures_v3.1.xlsx
#'
#' @docType data
#'
#' @usage data(SBS_SIGNATURE_PROFILES_V3)
'SBS_SIGNATURE_PROFILES_V3'

if(F){
   ## Code to create RData --------
   #df <- read.table('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/raw/COSMIC_v3.1_SBS_GRCh37.txt',sep='\t', header=1)
   df <- read.table('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/raw/COSMIC_v3.2_SBS_GRCh37.txt',sep='\t', header=1)
   df <- df[order(stringr::str_extract(df$Type,'\\w>\\w')),]


   m <- df[,grep('^SBS',colnames(df))]
   rownames(m) <- df$Type

   sig_metadata <- read.delim(
      '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/sigs_v3.2_metadata.txt'
   )
   excl_sigs <- subset(sig_metadata, nchar(exclude_reason)!=0, sig_name, drop=T)
   m <- m[,!(colnames(m) %in% excl_sigs)]

   #m[,'SBS1'] - SBS_SIGNATURE_PROFILES_V3[,'SBS1']
   #cbind(rownames(m), rownames(SBS_SIGNATURE_PROFILES_V3))

   SBS_SIGNATURE_PROFILES_V3 <- m
   save(
      SBS_SIGNATURE_PROFILES_V3,
      file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/SBS_SIGNATURE_PROFILES_V3.RData'
   )
}

