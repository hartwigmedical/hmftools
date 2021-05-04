package com.hartwig.hmftools.lilackt.candidates

import com.hartwig.hmftools.lilackt.SequenceCount
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci

class AminoAcidFiltering(private val aminoAcidBoundaries: Set<Int>) {

    fun aminoAcidCandidates(candidates: List<HlaSequenceLoci>, aminoAcidCount: SequenceCount): List<HlaSequenceLoci> {
        var result = candidates
        val locations = (0 until aminoAcidCount.length).toSet() subtract aminoAcidBoundaries
        for (loci: Int in locations) {
            //TODO: charles?
//            val depth = aminoAcidCount.depth(loci)
//            if (depth >= config.someParameter) {
                val expectedSequences = aminoAcidCount.sequenceAt(loci)
                result = result.filter { it.consistentWithAny(expectedSequences, loci) }
//            }
        }
        return result
    }



}