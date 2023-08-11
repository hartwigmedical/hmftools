package com.hartwig.hmftools.cider.blastn

import com.hartwig.hmftools.common.genome.region.Strand

data class BlastnMatch(
    val querySeqLen: Int,
    val subjectTitle: String,
    val percentageIdent: Double,
    val queryCoverage: Double,
    val alignmentLength: Int,
    val numMismatch: Int,
    val numGapOpenings: Int,
    val queryAlignStart: Int,
    val queryAlignEnd: Int,
    val subjectAlignStart: Int,
    val subjectAlignEnd: Int,
    val queryFrame: Int,
    val subjectFrame: Strand,
    val expectedValue: Double,
    val bitScore: Double,
    val alignedPartOfQuerySeq: String,
    val alignedPartOfSubjectSeq: String
)
