package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.lilac.LilacApplication
import com.hartwig.hmftools.lilac.SequenceCount

class NucleotideFragmentFactory(private val minBaseCount: Int, private val rawNucleotideFragments: List<NucleotideFragment>, private val aBoundaries: Set<Int>,
                                private val bBoundaries: Set<Int>, private val cBoundaries: Set<Int>) {

    fun allNucleotides(): List<NucleotideFragment> {
        val allProteinExonBoundaries = (aBoundaries + bBoundaries + cBoundaries)
        return allEvidence(allProteinExonBoundaries, rawNucleotideFragments, setOf(), setOf(), setOf())
    }

    fun typeANucleotides(): List<NucleotideFragment> {
        val inANotInB = aBoundaries subtract bBoundaries
        val inANotInC = aBoundaries subtract cBoundaries

        val inBNotInA = bBoundaries subtract aBoundaries
        val inCNotInA = cBoundaries subtract aBoundaries

        val bExcluded = inANotInB union inBNotInA
        val cExcluded = inANotInC union inCNotInA

        return allEvidence(aBoundaries, rawNucleotideFragments, setOf(), bExcluded, cExcluded)
    }

    fun typeBNucleotides(): List<NucleotideFragment> {
        val inBNotInA = bBoundaries subtract aBoundaries
        val inBNotInC = bBoundaries subtract cBoundaries

        val inANotInB = aBoundaries subtract bBoundaries
        val inCNotInB = cBoundaries subtract bBoundaries

        val aExcluded = inBNotInA union inANotInB
        val cExcluded = inBNotInC union inCNotInB

        return allEvidence(bBoundaries, rawNucleotideFragments, aExcluded, setOf(), cExcluded)
    }

    fun typeCNucleotides(): List<NucleotideFragment> {
        val inCNotInA = cBoundaries subtract aBoundaries
        val inCNotInB = cBoundaries subtract bBoundaries

        val inANotInC = aBoundaries subtract cBoundaries
        val inBNotInC = bBoundaries subtract cBoundaries

        val aExcluded = inCNotInA union inANotInC
        val bExcluded = inCNotInB union inBNotInC

        return allEvidence(cBoundaries, rawNucleotideFragments, aExcluded, bExcluded, setOf())
    }


    private fun allEvidence(
            aminoAcidBoundaries: Set<Int>, rawNucleotideFragments: List<NucleotideFragment>,
            aExcludedAminoAcidReads: Set<Int>, bExcludedAminoAcidReads: Set<Int>, cExcludedAminoAcidReads: Set<Int>): List<NucleotideFragment> {

        val aExcludedNucleotides = aExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
        val bExcludedNucleotides = bExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
        val cExcludedNucleotides = cExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }

        val filteredNucleotideFragments = rawNucleotideFragments
                .filter { !exclude(it, LilacApplication.HLA_A, aExcludedNucleotides) && !exclude(it, LilacApplication.HLA_B, bExcludedNucleotides) && !exclude(it, LilacApplication.HLA_C, cExcludedNucleotides) }
        val filteredNucleotideFragmentCounts = SequenceCount.nucleotides(minBaseCount, filteredNucleotideFragments)

        return NucleotideFragmentEnrichment(aminoAcidBoundaries.map { it * 3 }, filteredNucleotideFragmentCounts).enrichHomSpliceJunctions(filteredNucleotideFragments)

    }


    private fun exclude(fragment: NucleotideFragment, gene: String, excludedNucleotides: Collection<Int>): Boolean {
        return fragment.alignedGene == gene && excludedNucleotides.any { fragment.containsNucleotide(it) }
    }

}