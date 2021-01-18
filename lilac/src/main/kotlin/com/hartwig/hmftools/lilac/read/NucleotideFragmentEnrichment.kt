package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.nuc.SequenceCount

class NucleotideFragmentEnrichment(nucleotideExonBoundaryStarts: Collection<Int>, private var nucleotideCounts: SequenceCount) {

    private val homLoci = nucleotideCounts.homozygousIndices().toSet()
    private val homStarts = nucleotideExonBoundaryStarts.filter { homLoci.contains(it) }
    private val homEnds = nucleotideExonBoundaryStarts.filter { homLoci.contains(it + 1) && homLoci.contains(it + 2) }

    fun enrichHomSpliceJunctions(fragments: List<NucleotideFragment>): List<NucleotideFragment> {
        val result = mutableListOf<NucleotideFragment>()

        for (fragment in fragments) {
            var enriched = fragment

            for (homStart in homStarts) {

                if (missingStart(homStart, enriched)) {
                    enriched = enriched.addStart(homStart)
                }
            }

            for (homEnd in homEnds) {
                if (missingEnd(homEnd, enriched)) {
                    enriched = enriched.addEnd(homEnd)
                }
            }

            result.add(enriched)
        }

        return result
    }

    private fun missingStart(index: Int, fragment: NucleotideFragment): Boolean {
        return !fragment.containsNucleotide(index) && fragment.containsAllNucleotides(index + 1, index + 2)
    }

    private fun missingEnd(index: Int, fragment: NucleotideFragment): Boolean {
        return fragment.containsNucleotide(index) && !fragment.containsNucleotide(index + 1) && !fragment.containsNucleotide(index + 2)
    }

    private fun NucleotideFragment.addStart(index: Int): NucleotideFragment {
        return this.enrich(index, nucleotideCounts.sequenceAt(index).first())
    }

    private fun NucleotideFragment.addEnd(index: Int): NucleotideFragment {
        return this
                .enrich(index + 1, nucleotideCounts.sequenceAt(index + 1).first())
                .enrich(index + 2, nucleotideCounts.sequenceAt(index + 2).first())
    }
}