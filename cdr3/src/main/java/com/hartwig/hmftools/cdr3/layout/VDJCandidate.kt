package com.hartwig.hmftools.cdr3.layout

data class VDJCandidate(
    val vLayout: ReadLayout,
    val jLayout: ReadLayout,
    val vOverlapStart: Int, // where does the overlap start in V
    val jOverlapEnd: Int, // where does the overlap ends in J
    val overlapSeq: String
)
{
    val vdjSequence: String get()
    {
        val vSeq = vLayout.consensusSequence().substring(0, vOverlapStart)
        val jSeq = jLayout.consensusSequence().substring(jOverlapEnd)
        return vSeq + overlapSeq + jSeq
    }
}