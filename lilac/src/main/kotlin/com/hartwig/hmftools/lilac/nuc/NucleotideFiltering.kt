package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.lilac.seq.HlaSequence

class NucleotideFiltering(private val minBaseCount: Int, private val aminoAcidBoundaries: Set<Int>) {

    fun filterCandidatesOnAminoAcidBoundaries(candidates: Collection<HlaSequence>, fragments: List<NucleotideFragment>): List<HlaSequence> {
        var result = candidates.toList()
        for (aminoAcid in aminoAcidBoundaries) {
            result = filterCandidateOnAminoAcid(result, fragments, aminoAcid)
        }

        return result
    }

    private fun filterCandidateOnAminoAcid(candidates: Collection<HlaSequence>, fragments: List<NucleotideFragment>, aminoAcid: Int): List<HlaSequence> {
        val aminoAcidPredicate = aminoAcidPredicate(fragments, aminoAcid)
        return candidates.filter(aminoAcidPredicate)
    }

    private fun aminoAcidPredicate(fragments: List<NucleotideFragment>, aminoAcidIndex: Int): (HlaSequence) -> Boolean {
        val firstPredicate = nucleotidePredicate(fragments, aminoAcidIndex * 3)
        val remainingPredicate = nucleotidePredicate(fragments, aminoAcidIndex * 3 + 1, aminoAcidIndex * 3 + 2)
        return { x: HlaSequence -> firstPredicate(x) && remainingPredicate(x) }
    }

    private fun nucleotidePredicate(fragments: List<NucleotideFragment>, vararg nucleotideIndices: Int): (HlaSequence) -> Boolean {
        val fragmentsContainingSpecifiedNucleotides = fragments
                .filter { it.containsAllNucleotides(*nucleotideIndices) }
                .map { it.nucleotides(*nucleotideIndices) }
                .groupingBy { it }
                .eachCount()
                .filter { it.value >= minBaseCount }
                .keys
                .map { it.toCharArray() }

        return { x: HlaSequence -> fragmentsContainingSpecifiedNucleotides.isEmpty() || x.consistentWith(nucleotideIndices, fragmentsContainingSpecifiedNucleotides) }
    }

}