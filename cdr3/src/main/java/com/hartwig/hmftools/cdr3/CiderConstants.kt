package com.hartwig.hmftools.cdr3

object CiderConstants
{
    // following are constants for BLOSUM search
    const val MAX_BLOSUM_DIFF_PER_AA: Int = 2
    const val BLOSUM_SIMILARITY_SCORE_CONSTANT: Int = 5
    const val CANDIDATE_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 10
    const val VDJ_MIN_PARTIAL_ANCHOR_AA_LENGTH: Int = 1

    const val MIN_CANDIDATE_READ_ANCHOR_OVERLAP: Int = 6
    const val MIN_VJ_LAYOUT_JOIN_OVERLAP_BASES: Int = 20
    const val MIN_POLY_G_TRIM_COUNT: Int = 6
    const val POLY_G_TRIM_EXTRA_BASE_COUNT: Int = 5
    const val MIN_CDR3_LENGTH_BASES: Int = 9
}