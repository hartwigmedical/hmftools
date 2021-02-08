package com.hartwig.hmftools.lilac.candidates

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci

class AminoAcidFiltering(private val aminoAcidBoundaries: Set<Int>) {

    fun aminoAcidCandidates(candidates: List<HlaSequenceLoci>, aminoAcidCount: SequenceCount): List<HlaSequenceLoci> {
        var result = candidates
        val locations = (0 until aminoAcidCount.length).toSet() subtract aminoAcidBoundaries
        for (loci in locations) {
            val expectedSequences = aminoAcidCount.sequenceAt(loci)

            result = result.filter { it.consistentWithAny(expectedSequences, loci) }
        }
        return result
    }



}