package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles

class PhasedEvidenceFactory(private val minFragmentsPerAllele: Int, private val minFragmentsToRemoveSingles: Int, private val minEvidence: Int) {

    fun evidence(expectedAlleles: ExpectedAlleles, aminoAcidAminoAcidFragments: List<AminoAcidFragment>): List<PhasedEvidence> {
        val aminoAcidCounts = SequenceCount.aminoAcids(minEvidence, aminoAcidAminoAcidFragments)
        val heterozygousIndices = aminoAcidCounts.heterozygousLoci()
        println(heterozygousIndices)

        val heterozygousEvidence = ExtendEvidence(minFragmentsPerAllele, minFragmentsToRemoveSingles, heterozygousIndices, aminoAcidAminoAcidFragments, expectedAlleles)

        val finalisedEvidence = mutableSetOf<PhasedEvidence>()
        val unprocessedEvidence = mutableListOf<PhasedEvidence>()
        unprocessedEvidence.addAll(heterozygousEvidence.pairedEvidence())


        while (unprocessedEvidence.isNotEmpty()) {
            val top = unprocessedEvidence.removeAt(0)
            println("Processing: $top")

            val (parent, children) = heterozygousEvidence.merge(top, finalisedEvidence + unprocessedEvidence)

            if (children.isNotEmpty()) {
                println("Produced:   $parent")
                finalisedEvidence.removeAll(children)
                unprocessedEvidence.removeAll(children)
                unprocessedEvidence.add(parent)
            } else {
                finalisedEvidence.add(parent)
            }

            unprocessedEvidence.sort()
        }

//        return longestEvidence(finalisedEvidence)
        return finalisedEvidence.sortedBy { it.aminoAcidIndices[0] }

    }

    private fun longestEvidence(evidence: Collection<PhasedEvidence>): List<PhasedEvidence> {
        fun Collection<PhasedEvidence>.otherContains(victim: PhasedEvidence): Boolean {
            return this.any { it != victim && it.contains(victim) }
        }
        return evidence
                .filter { !evidence.otherContains(it) }
                .sortedBy { it.aminoAcidIndices[0] }
    }

}


