#' DBS signature profiles
#'
#' A matrix
#' rows: dinucleotide base substitution (DBS) context
#' cols: DBS signatures
#'
#' Source: https://cancer.sanger.ac.uk/sigs-assets-20/COSMIC_Mutational_Signatures_v3.1.xlsx
#'
#' @docType data
#'
#' @usage data(DBS_SIGNATURE_PROFILES)
'DBS_SIGNATURE_PROFILES'

# ## Code to create RData --------
# df <- openxlsx::read.xlsx(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/COSMIC_Mutational_Signatures_v3.1.xlsx',
#    sheet='DBS_GRCh37'
# )
#
# # identical(
# #    df$Type,
# #    rownames(DBS_SIGNATURE_PROFILES)
# # )
#
# m <- df[,grep('^DBS',colnames(df))]
# rownames(m) <- df$Type
#
# excl_sigs <- read.delim(
#    '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/sigs_v3.1_exclusion.txt'
# )
# excl_sigs <- subset(excl_sigs, nchar(exclude_reason)!=0, sig_name, drop=T)
# m <- m[,!(colnames(m) %in% excl_sigs)]
#
# DBS_SIGNATURE_PROFILES <- m
# save(
#    DBS_SIGNATURE_PROFILES,
#    file='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/data/DBS_SIGNATURE_PROFILES.RData'
# )
