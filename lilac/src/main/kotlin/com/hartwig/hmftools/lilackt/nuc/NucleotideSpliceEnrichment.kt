package com.hartwig.hmftools.lilackt.nuc

import com.hartwig.hmftools.lilackt.SequenceCount

/**
 * Add missing homozygous exon boundaries to fragments
 */
class NucleotideSpliceEnrichment(private val minBaseQuality: Int, private val minBaseCount: Int, private val aminoAcidBoundary: Set<Int>) {

    fun enrich(fragments: List<NucleotideFragment>): List<NucleotideFragment> {
        val filteredNucleotides = fragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }

        val nucleotideCounts = SequenceCount.nucleotides(minBaseCount, filteredNucleotides)
        val nucleotideExonBoundaryStarts = aminoAcidBoundary.map { 3 * it }
        val homLoci = nucleotideCounts.homozygousIndices().toSet()
        val homStarts = nucleotideExonBoundaryStarts.filter { homLoci.contains(it) }
        val homEnds = nucleotideExonBoundaryStarts.filter { homLoci.contains(it + 1) && homLoci.contains(it + 2) }

        val result = mutableListOf<NucleotideFragment>()

        for (fragment in fragments) {
            var enriched = fragment
            for (homStart in homStarts) {
                if (missingStart(homStart, enriched)) {
                    enriched = enriched.addStart(homStart, nucleotideCounts)
                }
            }

            for (homEnd in homEnds) {
                if (missingEnd(homEnd, enriched)) {
                    enriched = enriched.addEnd(homEnd, nucleotideCounts)
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

    private fun NucleotideFragment.addStart(index: Int, nucleotideCounts: SequenceCount): NucleotideFragment {
        return this.enrich(index, nucleotideCounts.sequenceAt(index).first(), minBaseQuality)
    }

    private fun NucleotideFragment.addEnd(index: Int, nucleotideCounts: SequenceCount): NucleotideFragment {
        return this
                .enrich(index + 1, nucleotideCounts.sequenceAt(index + 1).first(), minBaseQuality)
                .enrich(index + 2, nucleotideCounts.sequenceAt(index + 2).first(), minBaseQuality)
    }
}