package com.hartwig.hmftools.lilackt.nuc

import com.hartwig.hmftools.lilackt.SequenceCount

/**
 * Only permit high quality nucleotide values, ie, nucleotides that have at least 1 high quality (> [minBaseQuality]) and at least [minEvidence] instances in aggregate
 */
class NucleotideQualEnrichment(private val minBaseQuality: Int, private val minEvidence: Int) {

    fun enrich(nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        val highQualityFragments = nucleotideFragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        val highQualityCounts = SequenceCount.nucleotides(1, highQualityFragments)
        val rawCounts = SequenceCount.nucleotides(minEvidence, nucleotideFragments)
        val result = nucleotideFragments.map { enrich(it, highQualityCounts, rawCounts) }
        return result
    }

    private fun enrich(fragment: NucleotideFragment, highQualityCount: SequenceCount, rawCount: SequenceCount): NucleotideFragment {
        val nucleotideLoci = mutableListOf<Int>()
        val nucleotideQuality = mutableListOf<Int>()
        val nucleotides = mutableListOf<String>()

        for (i in fragment.nucleotideLoci().indices) {
            val loci = fragment.nucleotideLoci()[i]

            val currentQuality = fragment.nucleotideQuality()[i]
            val fragmentNucleotide = fragment.nucleotides()[i]
            val highQualitySequences = highQualityCount.sequenceAt(loci)
            val rawSequences = rawCount.sequenceAt(loci)
            val allowedSequences = rawSequences intersect highQualitySequences

            if (fragmentNucleotide in allowedSequences) {
                nucleotideLoci.add(loci)
                nucleotideQuality.add(currentQuality)
                nucleotides.add(fragmentNucleotide)
            }
        }

        return NucleotideFragment(fragment.id, fragment.genes, nucleotideLoci, nucleotideQuality, nucleotides)
    }

}