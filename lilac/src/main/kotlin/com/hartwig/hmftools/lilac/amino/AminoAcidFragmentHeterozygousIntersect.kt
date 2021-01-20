package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.read.Fragment

class AminoAcidFragmentHeterozygousIntersect(private val minAminoAcidCount: Int, private val minBaseQuality: Int) {

    fun heterozygousIntersect(nucleotideFragments: List<NucleotideFragment>): List<Fragment> {
        val qualityFilteredNucleotideFragments = nucleotideFragments.map { it.qualityFilter(minBaseQuality) }
        val qualityFilteredAminoAcidFragments = qualityFilteredNucleotideFragments.map { it.toAminoAcidFragment() }

        val aminoAcidCounts = SequenceCount.nucleotides(minAminoAcidCount, qualityFilteredAminoAcidFragments)
        val aminoAcidHeterozygousLoci = aminoAcidCounts.heterozygousLoci()
        val aminoAcidHeterozygousChar = aminoAcidHeterozygousLoci.map { aminoAcidCounts.sequenceAt(it) }


        val unfilteredAminoAcidFragments = nucleotideFragments.map { it.toAminoAcidFragment() }
        return unfilteredAminoAcidFragments.map { it.intersect(aminoAcidHeterozygousLoci, aminoAcidHeterozygousChar) }
    }

    private fun Fragment.intersect(aminoAcidHeterozygousLoci: List<Int>, aminoAcidHeterozygousChar: List<Collection<Char>>): Fragment {
        val initialIntersect = this.aminoAcidIndices() intersect aminoAcidHeterozygousLoci
        val filteredIntersect = initialIntersect.filter { loci ->
            val index = aminoAcidHeterozygousLoci.indexOf(loci)
            return@filter aminoAcidHeterozygousChar[index].contains(this.aminoAcid(loci))
        }

        return this.intersectAminoAcidLoci(filteredIntersect)
    }

}