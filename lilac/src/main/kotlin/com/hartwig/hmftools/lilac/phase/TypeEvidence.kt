package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.LilacApplication2
import com.hartwig.hmftools.lilac.nuc.SequenceCount
import com.hartwig.hmftools.lilac.read.NucleotideFragment
import com.hartwig.hmftools.lilac.read.NucleotideFragmentEnrichment

class TypeEvidence(private val minBaseCount: Int, private val minFragmentCount: Int, private val rawNucleotideFragments: List<NucleotideFragment>,
                   private val aBoundaries: Set<Int>, private val bBoundaries: Set<Int>, private val cBoundaries: Set<Int>) {


    fun typeAEvidence(): List<PhasedEvidence> {
        val inANotInB = aBoundaries subtract bBoundaries
        val inANotInC = aBoundaries subtract cBoundaries

        val inBNotInA = bBoundaries subtract aBoundaries
        val inCNotInA = cBoundaries subtract aBoundaries

        val bExcluded = inANotInB union inBNotInA
        val cExcluded = inANotInC union inCNotInA

        return allEvidence(aBoundaries, rawNucleotideFragments, setOf(), bExcluded, cExcluded)
    }

    fun typeBEvidence(): List<PhasedEvidence> {
        val inBNotInA = bBoundaries subtract aBoundaries
        val inBNotInC = bBoundaries subtract cBoundaries

        val inANotInB = aBoundaries subtract bBoundaries
        val inCNotInB = cBoundaries subtract bBoundaries

        val aExcluded = inBNotInA union inANotInB
        val cExcluded = inBNotInC union inCNotInB

        return allEvidence(aBoundaries, rawNucleotideFragments, aExcluded, setOf(), cExcluded)
    }

    fun typeCEvidence(): List<PhasedEvidence> {
        val inCNotInA = cBoundaries subtract aBoundaries
        val inCNotInB = cBoundaries subtract bBoundaries

        val inANotInC = aBoundaries subtract cBoundaries
        val inBNotInC = bBoundaries subtract cBoundaries

        val aExcluded = inCNotInA union inANotInC
        val bExcluded = inCNotInB union inBNotInC

        return allEvidence(aBoundaries, rawNucleotideFragments, aExcluded, bExcluded, setOf())
    }


    private fun allEvidence(
            aminoAcidBoundaries: Set<Int>, rawNucleotideFragments: List<NucleotideFragment>,
            aExcludedAminoAcidReads: Set<Int>, bExcludedAminoAcidReads: Set<Int>, cExcludedAminoAcidReads: Set<Int>): List<PhasedEvidence> {

        val aExcludedNucleotides = aExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
        val bExcludedNucleotides = bExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
        val cExcludedNucleotides = cExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }

        val filteredNucleotideFragments = rawNucleotideFragments
                .filter { !exclude(it, LilacApplication2.HLA_A, aExcludedNucleotides) && !exclude(it, LilacApplication2.HLA_B, bExcludedNucleotides) && !exclude(it, LilacApplication2.HLA_C, cExcludedNucleotides) }
        val filteredNucleotideFragmentCounts = SequenceCount.nucleotides(minBaseCount, filteredNucleotideFragments)

        val enrichedNucleotideFragments = NucleotideFragmentEnrichment(aminoAcidBoundaries.map { it * 3 }, filteredNucleotideFragmentCounts).enrichHomSpliceJunctions(filteredNucleotideFragments)
        val aminoAcidFragments = enrichedNucleotideFragments.map { it.toAminoAcidFragment() }
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)

        val heterozygousIndices = aminoAcidCounts.heterozygousIndices()
        val heterozygousEvidence = ExtendedEvidence(minBaseCount, minFragmentCount, heterozygousIndices, aminoAcidFragments)

        val allEvidence = mutableSetOf<PhasedEvidence>()
        val initialEvidence = heterozygousEvidence.initialEvidence()
        var unprocessedEvidence = initialEvidence

        allEvidence.addAll(initialEvidence)

        while (unprocessedEvidence.isNotEmpty()) {
            val topEvidence = unprocessedEvidence[0]
            allEvidence.add(topEvidence)

            val newEvidence = heterozygousEvidence.extendConsecutive(topEvidence, allEvidence)
            allEvidence.addAll(newEvidence)

            val updatedEvidence = mutableSetOf<PhasedEvidence>()
            updatedEvidence.addAll(unprocessedEvidence.drop(1))
            updatedEvidence.addAll(newEvidence)

            unprocessedEvidence = updatedEvidence.sorted()
        }

        return longestEvidence(allEvidence)

    }

    private fun longestEvidence(evidence: Collection<PhasedEvidence>): List<PhasedEvidence> {
        fun Collection<PhasedEvidence>.otherContains(victim: PhasedEvidence): Boolean {
            return this.any { it != victim && it.contains(victim) }
        }
        return evidence
                .filter { !evidence.otherContains(it) }
                .sortedBy { it.aminoAcidIndices[0] }
    }


    private fun exclude(fragment: NucleotideFragment, gene: String, excludedNucleotides: Collection<Int>): Boolean {
        return fragment.alignedGene == gene && excludedNucleotides.any { fragment.containsNucleotide(it) }
    }

}


