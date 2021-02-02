package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment

/**
 * Only permit high quality amino acids, ie, amino acids that have at least [minEvidence] instances with quality > [minBaseQuality] in aggregate
 */
class AminoAcidQualEnrichment(private val minBaseQuality: Int, private val minEvidence: Int) {

    fun enrich(nucleotideFragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val qualityFilteredNucleotideFragments = nucleotideFragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        val qualityFilteredAminoAcidFragments = qualityFilteredNucleotideFragments.map { it.toAminoAcidFragment() }
        val highQualityAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, qualityFilteredAminoAcidFragments)

        val unfilteredAminoAcidFragments = nucleotideFragments
                .filter { it.isNotEmpty() }
                .map { it.toAminoAcidFragment() }
        val result = unfilteredAminoAcidFragments.map { enrich(it, highQualityAminoAcidCounts) }

        return result
    }

    private fun enrich(fragment: AminoAcidFragment, count: SequenceCount): AminoAcidFragment {
        val initialIntersect = fragment.aminoAcidIndices()
        val filteredIntersect = initialIntersect.filter { loci ->
            val allowed = count.sequenceAt(loci)
            val actual = fragment.aminoAcid(loci)

            return@filter allowed.contains(actual)
        }

        return fragment.intersectAminoAcidLoci(filteredIntersect)
    }

}