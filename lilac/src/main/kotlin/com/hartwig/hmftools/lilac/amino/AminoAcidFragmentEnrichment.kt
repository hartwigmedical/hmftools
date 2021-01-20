package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment

class AminoAcidFragmentEnrichment(private val minBaseQuality: Int, private val minAminoAcidCount: Int) {

    fun heterozygousIntersect(nucleotideFragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val qualityFilteredNucleotideFragments = nucleotideFragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        val qualityFilteredAminoAcidFragments = qualityFilteredNucleotideFragments.map { it.toAminoAcidFragment() }

        val highQualityAminoAcidCounts = SequenceCount.aminoAcids(minAminoAcidCount, qualityFilteredAminoAcidFragments)
        highQualityAminoAcidCounts.writeVertically("/Users/jon/hmf/analysis/hla/aminoacids.count.hla-a.before.txt")

        val unfilteredAminoAcidFragments = nucleotideFragments.map { it.toAminoAcidFragment() }
        val result =  unfilteredAminoAcidFragments.map { it.intersect(highQualityAminoAcidCounts) }
        val updatedAminoAcidCounts = SequenceCount.aminoAcids(minAminoAcidCount, result)
        updatedAminoAcidCounts.writeVertically("/Users/jon/hmf/analysis/hla/aminoacids.count.hla-a.after.txt")

        return result
    }

    private fun AminoAcidFragment.intersect(count: SequenceCount): AminoAcidFragment {
        val initialIntersect = this.aminoAcidIndices()
        val filteredIntersect = initialIntersect.filter { loci ->
            val allowed = count.sequenceAt(loci)
            val actual = this.aminoAcid(loci)

            return@filter allowed.contains(actual)
        }

        return this.intersectAminoAcidLoci(filteredIntersect)
    }

    private fun AminoAcidFragment.intersect(aminoAcidHeterozygousLoci: List<Int>, aminoAcidHeterozygousChar: List<Collection<Char>>): AminoAcidFragment {
        val initialIntersect = this.aminoAcidIndices()
        val filteredIntersect = initialIntersect.filter { loci ->
            if (loci == 10) {
                println("SDf")
            }

            val allowed = aminoAcidHeterozygousChar[aminoAcidHeterozygousLoci.indexOf(loci)]
            val actual = this.aminoAcid(loci)

            return@filter allowed.contains(actual)
        }

        return this.intersectAminoAcidLoci(filteredIntersect)
    }

}