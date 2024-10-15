#' Various transformation of contexts (e.g. into mutational signatures)
#'
#' @description Performs the following transformations of each mut type (snv, dbs, indel, sv) in the
#' following order: simplify, least-squares fitting, calculate relative counts. Simplify flattens
#' the 96 snv contexts into 6 contexts, indel contexts into their types (mh, rep, none),
#' and SV type/length contexts into SV type. Least-squares fitting converts contexts into
#' signatures.
#'
#' @param contexts A list containing data frames of snv, indel, and/or sv contexts. Alternatively,
#' these can be specified with the arguments: snv, indel, or sv. This is handy if one only wants to
#' transform a matrix of a certain mut type
#' @param snv See argument 'contexts'.
#' @param indel See argument 'contexts'.
#' @param sv See argument 'contexts'.
#' @param simplify.types Which types to flatten. Accepts: 'snv','indel','sv','all'
#' @param lsqnonneg.types Which types to fit to signatures. Accepts: 'snv', 'sv','all'
#' @param rel.types Which types to convert to relative contribution. Accepts: 'snv','indel','sv','all'
#' @param sig.profiles.snv SNV sig profiles used for lsqnonneg. Defaults to 30 COSMIC signatures
#' @param sig.profiles.dbs DBS sig profiles used for lsqnonneg. Defaults to PCAWG DBS signatures
#' @param sig.profiles.sv SV sig profiles used for lsqnonneg. Defaults to 6 SV signatures from 560
#' breast cancer paper
#' @param export.list Output a list with the split mutation types rather than a matrix?
#'
#' @return A matrix or data frame
#' @export
#'
#' @examples
#' base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/Breast_Organoids/matrices/'
#' contexts <- list(
#'    snv = readSigsAsDf(paste0(base_dir,'/snv_contexts')),
#'    indel = readSigsAsDf(paste0(base_dir,'/indel')),
#'    sv = readSigsAsDf(paste0(base_dir,'/sv_contexts'))
#' )
#' transformContexts(contexts, simplify.types=c('snv','indel','sv'), rel.types=c('snv','indel','sv'))

transformContexts <- function(
   contexts=NULL, snv=NULL, dbs=NULL, indel=NULL, sv=NULL,
   simplify.types=NULL, lsqnonneg.types=NULL, rel.types=NULL,

   sig.profiles.snv=SBS_SIGNATURE_PROFILES_V2,
   sig.profiles.dbs=DBS_SIGNATURE_PROFILES,
   sig.profiles.sv=SV_SIGNATURE_PROFILES,

   export.list=F
){

   if(!is.null(contexts)){
      snv <- if(!is.null(contexts$snv)){ contexts$snv }
      dbs <- if(!is.null(contexts$dbs)){ contexts$dbs }
      indel <- if(!is.null(contexts$indel)){ contexts$indel }
      sv <- if(!is.null(contexts$sv)){ contexts$sv }
   }

   out <- list()

   if(!is.null(snv)){
      if(any(c('snv','all') %in% simplify.types)){
         snv_split <- splitDfRegex(SUBSTITUTIONS, df = snv)
         snv <- do.call(cbind, lapply(snv_split, rowSums))
         colnames(snv) <- gsub('>','.',SUBSTITUTIONS)
      }

      if(any(c('snv','all') %in% lsqnonneg.types)){
         snv <- fitToSignatures(mut.context.counts=contexts$snv, signature.profiles=sig.profiles.snv)
      }

      if(any(c('snv','all') %in% rel.types)){
         snv <- snv/rowSums(snv)
         snv[is.na(snv)] <- 0
      }

      out$snv <- snv
   }

   if(!is.null(dbs)){

      if(any(c('dbs','all') %in% lsqnonneg.types)){
         dbs <- fitToSignatures(mut.context.counts=contexts$dbs, signature.profiles=sig.profiles.dbs)
      }

      if(any(c('dbs','all') %in% rel.types)){
         dbs <- dbs/rowSums(dbs)
         dbs[is.na(dbs)] <- 0
      }

      out$dbs <- dbs
   }

   if(!is.null(indel)){
      if(any(c('indel','all') %in% simplify.types)){
         indel_types <- c('del.rep', 'ins.rep', 'del.mh', 'ins.mh', 'del.none', 'ins.none')
         indel_split <- splitDfRegex(indel_types, df = indel)
         indel <- do.call(cbind, lapply(indel_split, rowSums))
         colnames(indel) <- indel_types
      }

      if('indel' %in% lsqnonneg.types){
         stop('Indel least-squares fit signatures have not been implemented in this version')
      }

      if(any(c('indel','all') %in% rel.types)){
         indel <- indel/rowSums(indel)
         indel[is.na(indel)] <- 0
      }

      out$indel <- indel
   }

   if(!is.null(sv)){
      if(any(c('sv','all') %in% simplify.types)){
         sv_types <- c('DEL', 'DUP', 'INV', 'TRA')
         sv_split <- splitDfRegex(sv_types, df = sv)
         sv <- do.call(cbind, lapply(sv_split, rowSums))
         colnames(sv) <- paste0('SV.',sv_types)
      }

      if(any(c('sv','all') %in% lsqnonneg.types)){
         sv <- fitToSignatures(mut.context.counts=contexts$sv, signature.profiles=sig.profiles.sv)
      }

      if(any(c('sv','all') %in% rel.types)){
         sv <- sv/rowSums(sv)
         sv[is.na(sv)] <- 0
      }

      out$sv <- sv
   }
   if(export.list){ return(out) }
   return( do.call(cbind, unname(out)) )

}


