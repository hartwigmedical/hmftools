package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.lilac.SequenceCount

/**
 * Only permit high quality nucleotide values, ie, nucleotides that have at least [minEvidence] instances with quality > [minBaseQuality] in aggregate
 */
class NucleotideQualEnrichment(private val minBaseQuality: Int, private val minEvidence: Int) {

    fun enrich(nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        val qualityFilteredFragments = nucleotideFragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        val counts = SequenceCount.nucleotides(minEvidence, qualityFilteredFragments)
        val result = nucleotideFragments.map { enrich(it, counts) }
        return result
    }

    private fun enrich(fragment: NucleotideFragment, count: SequenceCount): NucleotideFragment {
        val nucleotideLoci = mutableListOf<Int>()
        val nucleotideQuality = mutableListOf<Int>()
        val nucleotides = mutableListOf<Char>()

        for (i in fragment.nucleotideLoci().indices) {
            val loci = fragment.nucleotideLoci()[i]

            val currentQuality = fragment.nucleotideQuality()[i]
            val fragmentNucleotide = fragment.nucleotides()[i]
            val countNucleotides = count.sequenceAt(loci)

            if (fragmentNucleotide in countNucleotides) {
                val updatedQuality = currentQuality.coerceAtLeast(minBaseQuality)
                nucleotideLoci.add(loci)
                nucleotideQuality.add(updatedQuality)
                nucleotides.add(fragmentNucleotide)
            }
        }

        return NucleotideFragment(fragment.id, fragment.genes, nucleotideLoci, nucleotideQuality, nucleotides)
    }

}