####################################################################################################
#' Get the start/end positions for the left/right flanks of an indel
#'
#' @description Helper function for extractSigsIndel(), getting flank start/end positions first, then retrieving the sequences
#' using getSeq with a vectorized input improves speed significantly.
#'
#' @param chrom A character stating the chromosome
#' @param pos The REF position as an integer
#' @param indel.len The length of the indel as an integer
#' @param indel.type A character stating whether the indel is an insertion ('ins') or deletion ('del')
#' @param n.indel.lengths.l The length of the flanking sequence to return (measured in indel lengths)
#' @param n.indel.lengths.r See n.indel.lengths.l
#'
#' @return A vector of the left flank start/end position and right flank start/end position
indelSeqFlanksStartEnd <- function(chrom, pos, indel.len, indel.type, n.indel.lengths.l=1, n.indel.lengths.r=1){
   if(indel.type == 'del'){
      r_flank <- c(
         r_start = pos + 1 + indel.len,
         r_end = pos + indel.len + indel.len*n.indel.lengths.r
      )
   } else if(indel.type == 'ins'){
      r_flank <- c(
         r_start = pos + 1,
         r_end = pos + indel.len*n.indel.lengths.r
      )
   }

   l_flank <- c(
      l_start = pos - indel.len*n.indel.lengths.l + 1,
      l_end = pos
   )

   out <- l_flank
   out <- c(out, r_flank)
   return(out)
}

####################################################################################################
#' Calculate the number of copies of the indel sequence are present in the flank sequence
#'
#' @description Helper function for extractSigsIndel(). Scans along the flanking sequence for
#' copies of the indel sequence using a sliding window with a length equal to the length of
#' the indel sequence. The sliding window jumps in increments equal to the indel length.
#' Scanning stops prematurely when the indel sequence does not match the one in the sliding
#' window. The indel sequence itself counts as one copy.
#'
#' Scanning needs only to be done from the 5' -> 3' direction as the reported POS of an indel
#' in a repeat region is always the first position of the repeat region.
#'
#' @param indel.seq The indel sequence as a character
#' @param flank.seq The flanking sequence as a character
#'
#' @return An integer stating the number of copies found in the flanking sequence
nCopiesAlongFlank <- function(indel.seq, flank.seq){
   # indel.seq = "T"
   # flank.seq = "TTTTGCG"

   indel_length <- nchar(indel.seq)
   rail_seq <- paste0(indel.seq, flank.seq)

   count <- 0
   seq_window <- substr(rail_seq, 1, indel_length)
   while(indel.seq == seq_window){
      count <- count + 1
      seq_window <- substr(
         rail_seq,
         count*indel_length + 1,
         count*indel_length + indel_length
      )
   }

   return(count)
}

####################################################################################################
#' Calculate the number of bases that are homologous to the 3' flanking sequence
#'
#' @description Helper function for extractSigsIndel(). Scans a maximum of 1 indel length from
#' 5' to 3' in the flanking sequence to detect identical bases. Stops prematurely if a
#' non-identical base is found.
#'
#' DSBs can be repaired using 3' microhomology, which can be present relative to either the +
#' or - strand. Therefore, both left and right flanking sequences of the indel need to be
#' considered. When calculating left flanking homology (i.e. homology in the 3' direction for
#' the - strand), the reverse complement of the indel sequence and flanking sequence need to
#' be taken. However, the reverse can be taken for the calculation to save computation.
#'
#' @param indel.seq The indel sequence as a character
#' @param flank.seq The flanking sequence as a character
#'
#' @return An integer stating the number of bases in microhomology
nBasesMH <- function(indel.seq, flank.seq, indel.len){
   #indel.sequence = "CTA"
   #flank.sequence = "C"

   indel_len <- nchar(indel.seq)
   indel.seq <- unlist(strsplit(indel.seq, ''))
   flank.seq <- unlist(strsplit(flank.seq, ''))[1:indel_len]

   n_bases <- 0
   for(i in 1:length(indel.seq)){
      if(indel.seq[i] != flank.seq[i]){
         break
      } else {
         n_bases <- n_bases + 1
      }
   }
   return(n_bases)
}



