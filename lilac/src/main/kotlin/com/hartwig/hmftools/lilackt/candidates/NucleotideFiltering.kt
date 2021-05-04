package com.hartwig.hmftools.lilackt.candidates

import com.hartwig.hmftools.lilackt.nuc.NucleotideFragment
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci

class NucleotideFiltering(private val minNucleotideCount: Int, private val aminoAcidBoundaries: Set<Int>) {

    fun filterCandidatesOnAminoAcidBoundaries(candidates: Collection<HlaSequenceLoci>, fragments: List<NucleotideFragment>): List<HlaSequenceLoci> {
        var result = candidates.toList()
        for (aminoAcid in aminoAcidBoundaries) {
            val nucleotideStart = aminoAcid * 3
            val startSequences = nucleotideSequence(fragments, nucleotideStart)
            val endSequences = nucleotideSequence(fragments, nucleotideStart + 1, nucleotideStart + 2)

            result = result.filter { it.consistentWithAny(nucleotideStart, startSequences, endSequences) }
        }

        return result
    }

    private fun HlaSequenceLoci.consistentWithAny(startLoci: Int, startSequences: List<String>, endSequences: List<String>): Boolean {
        return consistentWithAny(startSequences, startLoci) && consistentWithAny(endSequences, startLoci + 1, startLoci + 2)
    }

    private fun nucleotideSequence(fragments: List<NucleotideFragment>, vararg nucleotideIndices: Int): List<String> {
        return fragments
                .filter { it.containsAllNucleotides(*nucleotideIndices) }
                .map { it.nucleotides(*nucleotideIndices) }
                .groupingBy { it }
                .eachCount()
                .filter { it.value >= minNucleotideCount }
                .map { it.key }

    }


}