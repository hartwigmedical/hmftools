package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion

object CiderConstants
{
    // for partially rearranged VDJs, we designate that the "interesting"
    // sequence is 60 bases before the J anchor or 60 bases after the V anchor
    const val PARTIAL_VDJ_UNANCHORED_LENGTH_BASES: Int = 60

    // following are constants for BLOSUM search
    const val MAX_BLOSUM_DIFF_PER_AA: Int = 3
    const val BLOSUM_SIMILARITY_SCORE_CONSTANT: Int = 6
    const val BLOSUM_UNKNOWN_BASE_PENALTY: Int = 2

    const val CANDIDATE_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 10
    const val VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 3

    const val MAX_READ_DISTANCE_FROM_ANCHOR: Int = 50

    // minimum amount of bases a read must overlap with a layout for it to be added
    const val LAYOUT_MIN_READ_OVERLAP_BASES: Int = 20
    const val MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES: Int = 20

    // word size of 8 is experimentally decided to provide good speed up in finding VJ overlaps
    const val VJ_JOIN_HASH_WORD_SIZE: Int = 8
    const val MIN_ANCHOR_LENGTH_BASES: Int = 3
    const val LAYOUT_MIN_SUPPORT_TO_SEAL_NODE: Int = 5

    // minimum fraction of reads that support a certain base in a layout
    // this is very important such that we do not let one or two reads stop
    // layouts from being joined together
    const val MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION: Double = 0.02
    const val MIN_POLY_G_TRIM_COUNT: Int = 2
    const val POLY_G_TRIM_EXTRA_BASE_COUNT: Int = 5

    const val MIN_NON_SPLIT_READ_STRADDLE_LENGTH: Int = 30

    // filter out matches that have too low identity
    // reason for doing this is that we use match/mismatch of 1/-4, in worst case we can
    // get 1 mismatch for every 4 matches, and could find alignments with 80% identity.
    // those are probably too different to use. We use 90% for V / J identities, and 95%
    // cut off for full match
    const val ALIGNMENT_MATCH_MIN_VJ_IDENTITY = 90
    const val ALIGNMENT_MATCH_FULL_MATCH_IDENTITY = 95

    // blast uses v38
    val BLAST_REF_GENOME_VERSION = RefGenomeVersion.V38
}