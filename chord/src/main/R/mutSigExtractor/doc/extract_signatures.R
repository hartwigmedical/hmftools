#' Look up table for doublet base substitution types and reverse complement sequences
#'
#' Excel sheet source: https://cancer.sanger.ac.uk/signatures/DBS_vignettes/DBS-doublet-base-substitution-classification.xlsx
#'
#' @docType data
#'
#' @usage data(DBS_TYPES)
'DBS_TYPES'

# ## Code to create RData file
# library(openxlsx)
# DBS_TYPES <- read.xlsx(
#    '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/raw/DBS-doublet-base-substitution-classification.xlsx',
#    sheet='proc'
# )
#
# DBS_TYPES <- DBS_TYPES[,c('ref','context','context_rev_comp')]
#
# save(DBS_TYPES, file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/data/DBS_TYPES.RData')
