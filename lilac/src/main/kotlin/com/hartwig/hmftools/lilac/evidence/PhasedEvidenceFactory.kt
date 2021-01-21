package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment

class PhasedEvidenceFactory(private val minEvidence: Int, private val minFragments: Int) {

    fun evidence(aminoAcidAminoAcidFragments: List<AminoAcidFragment>): List<PhasedEvidence> {
        val aminoAcidCounts = SequenceCount.aminoAcids(minEvidence, aminoAcidAminoAcidFragments)

        val heterozygousIndices = aminoAcidCounts.heterozygousLoci()
        val heterozygousEvidence = ExtendEvidence(minFragments, heterozygousIndices, aminoAcidAminoAcidFragments)

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

}


