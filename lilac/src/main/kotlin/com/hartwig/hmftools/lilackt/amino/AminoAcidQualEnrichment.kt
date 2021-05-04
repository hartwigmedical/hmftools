package com.hartwig.hmftools.lilackt.amino

import com.hartwig.hmftools.lilackt.SequenceCount
import com.hartwig.hmftools.lilackt.nuc.NucleotideFragment

/**
 * Only permit high quality amino acids, ie, amino acids that have at least [minEvidence]
 */
class AminoAcidQualEnrichment(private val minEvidence: Int) {

    fun enrich(nucleotideFragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val qualityFilteredAminoAcidFragments = nucleotideFragments.map { it.toAminoAcidFragment() }
        val highQualityAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, qualityFilteredAminoAcidFragments)

        val unfilteredAminoAcidFragments = nucleotideFragments
                .filter { it.isNotEmpty() }
                .map { it.toAminoAcidFragment() }
        val result = unfilteredAminoAcidFragments.map { enrich(it, highQualityAminoAcidCounts) }

        return result
    }

    private fun enrich(fragment: AminoAcidFragment, count: SequenceCount): AminoAcidFragment {
        val initialIntersect = fragment.aminoAcidLoci()
        val filteredIntersect = initialIntersect.filter { loci ->
            val allowed = count.sequenceAt(loci)
            val actual = fragment.aminoAcid(loci)

            return@filter allowed.contains(actual)
        }

        return fragment.intersectAminoAcidLoci(filteredIntersect)
    }

}