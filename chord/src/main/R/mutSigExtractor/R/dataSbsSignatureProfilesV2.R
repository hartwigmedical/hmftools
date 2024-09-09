#' SBS signature profiles (original 30)
#'
#' A matrix
#' rows: 96-trinucleotide context
#' cols: 30 single base substitution (SBS) signatures
#'
#' Source: https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#'
#' @docType data
#'
#' @usage data(SBS_SIGNATURE_PROFILES_V2)
'SBS_SIGNATURE_PROFILES_V2'

## Code to create RData file
# sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
# cancer_signatures_raw <- read.table(sp_url, sep = "\t", header = TRUE)
# cancer_signatures_raw <- cancer_signatures_raw[order(cancer_signatures_raw[,1]),]
#
# cancer_signatures <- as.matrix(cancer_signatures_raw[,4:33])
# colnames(cancer_signatures) <- paste0('SBS', formatC(1:30, width = 2, format = "d", flag = "0"))
# rownames(cancer_signatures) <- cancer_signatures_raw$Somatic.Mutation.Type
#
# SBS_SIGNATURE_PROFILES_V2 <- cancer_signatures
# save(SBS_SIGNATURE_PROFILES_V2, file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/SBS_SIGNATURE_PROFILES_V2.RData')


