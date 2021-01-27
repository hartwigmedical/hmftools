package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.lilac.seq.HlaSequence

class NucleotideFiltering(private val minBaseCount: Int, private val aminoAcidBoundaries: Set<Int>) {

    fun filterCandidatesOnAminoAcidBoundaries(candidates: Collection<HlaSequence>, fragments: List<NucleotideFragment>): List<HlaSequence> {
        var result = candidates.toList()
        for (aminoAcid in aminoAcidBoundaries) {
            val nucleotideStart = aminoAcid * 3
            val startSequences = nucleotideSequence(fragments, nucleotideStart)
            val endSequences = nucleotideSequence(fragments, nucleotideStart + 1, nucleotideStart  + 2)

            result = result.filter { it.consistentWith(nucleotideStart, startSequences, endSequences) }
        }

        return result
    }

    private fun HlaSequence.consistentWith(startLoci: Int, startSequences: List<CharArray>, endSequences: List<CharArray>): Boolean {
        val startIsConsistent = startSequences.isEmpty() || this.consistentWith(listOf(startLoci).toIntArray(), startSequences)
        val endIsConsistent = endSequences.isEmpty() || this.consistentWith(listOf(startLoci + 1, startLoci + 2).toIntArray(), endSequences)
        return startIsConsistent && endIsConsistent
    }



    private fun nucleotideSequence(fragments: List<NucleotideFragment>, vararg nucleotideIndices: Int): List<CharArray> {
        return fragments
                .filter { it.containsAllNucleotides(*nucleotideIndices) }
                .map { it.nucleotides(*nucleotideIndices) }
                .groupingBy { it }
                .eachCount()
                .filter { it.value >= minBaseCount }
                .keys
                .map { it.toCharArray() }

    }


}