#' Indel signature profiles
#'
#' A matrix
#' rows: indel context
#' cols: indel signature
#'
#' Source: https://cancer.sanger.ac.uk/sigs-assets-20/COSMIC_Mutational_Signatures_v3.1.xlsx
#'
#' @docType data
#'
#' @usage data(INDEL_SIGNATURE_PROFILES)
'INDEL_SIGNATURE_PROFILES'

# ## Code to create RData --------
# df <- openxlsx::read.xlsx(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/COSMIC_Mutational_Signatures_v3.1.xlsx',
#    sheet='ID_GRCh37'
# )
#
# m <- df[!is.na(df$Type2),]
# rownames(m) <- m$Type2
# m <- m[,grep('^ID',colnames(df))]
# m <- as.matrix(m)
#
# excl_sigs <- read.delim(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/sigs_v3.1_exclusion.txt'
# )
# excl_sigs <- subset(excl_sigs, nchar(exclude_reason)!=0, sig_name, drop=T)
# m <- m[,!(colnames(m) %in% excl_sigs)]
#
# INDEL_SIGNATURE_PROFILES <- m
# save(
#    INDEL_SIGNATURE_PROFILES,
#    file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/INDEL_SIGNATURE_PROFILES.RData'
# )


