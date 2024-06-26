package com.hartwig.hmftools.teal

import htsjdk.samtools.util.SequenceUtil

object TealConstants
{
    const val DEFAULT_PARTITION_SIZE = 10_000_000
    const val DEFAULT_MIN_TELE_SEQ_COUNT = 2
    const val POLY_G_THRESHOLD = 0.9
    const val MIN_TELOMERE_MATCH_BASES = 12
    const val CANONICAL_TELOMERE_SEQ = "TTAGGG"

    // NOTE: cannot use TealUtils.reverseComplementSeq, cause we would run into static initialisation
    // ordering issue, since TealUtils uses TealConstants
    val CANONICAL_TELOMERE_SEQ_REV = SequenceUtil.reverseComplement("TTAGGG")

    /*
    val TELOMERE_HEXAMERS: Array<String> = arrayOf(
        CANONICAL_TELOMERE_SEQ,
        "TCAGGG",
        "TTCGGG",
        "GTAGGG",
        "TGAGGG",
        "TTGGGG",
        "TAAGGG",
        "ATAGGG",
        "CTAGGG",
        "TTTGGG"
    )
    val TELOMERE_HEXAMERS_REV: Array<String> = TELOMERE_HEXAMERS.map({ obj: String -> SequenceUtil.reverseComplement(obj) }).toTypedArray()
     */

    val CANONICAL_TELOMERE_SEQUENCES = arrayOf(
        CANONICAL_TELOMERE_SEQ.repeat(DEFAULT_MIN_TELE_SEQ_COUNT),
        CANONICAL_TELOMERE_SEQ_REV.repeat(DEFAULT_MIN_TELE_SEQ_COUNT))
}