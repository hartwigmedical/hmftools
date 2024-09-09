#' SV signature profiles
#'
#' A matrix.
#' rows: For deletions, duplications and inversions, SV type/SV length context. Translocations take up one row as they do not
#' have length.
#' cols: 6 SV signatures
#'
#' Source: https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip
#' (Supplementary.Table.21.Signatures.v3.xlsx)
#'
#' Note that the probabilities of the clustered and non-clustered rearrangements have been combined. In other words, whether the
#' rearrangements were clustered/non-clustered were not considered.
#'
#' @docType data
#'
#' @usage data(SV_SIGNATURE_PROFILES)
'SV_SIGNATURE_PROFILES'
