package com.hartwig.hmftools.cider

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
    const val VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 1

    const val MAX_READ_DISTANCE_FROM_ANCHOR: Int = 50

    // minimum amount of bases a read must overlap with a layout for it to be added
    const val LAYOUT_MIN_READ_OVERLAP_BASES: Int = 20
    const val MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES: Int = 20
    const val MIN_ANCHOR_LENGTH_BASES: Int = 3

    // minimum fraction of reads that support a certain base in a layout
    const val MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION: Double = 0.02
    const val MIN_POLY_G_TRIM_COUNT: Int = 6
    const val POLY_G_TRIM_EXTRA_BASE_COUNT: Int = 5
    const val MIN_CDR3_LENGTH_BASES: Int = 9

    const val MIN_NON_SPLIT_READ_STRADDLE_LENGTH: Int = 30

    enum class VjAnchorTemplateTsvColumn
    {
        gene,
        allele,
        chr,
        posStart,
        posEnd,
        strand,
        anchorStart,
        anchorEnd,
        anchorSequence,
        anchorAA,
        sequence
    }
}