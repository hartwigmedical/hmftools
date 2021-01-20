package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import java.util.*

class ExtendEvidence(private val minTotalEvidence: Int, private val heterozygousLoci: List<Int>, private val aminoAcidFragments: List<AminoAcidFragment>) {

    fun initialEvidence(): List<PhasedEvidence> {
        val result = mutableListOf<PhasedEvidence>()

        for (i in 0..(heterozygousLoci.size - 2)) {
            val evidence = PhasedEvidence.evidence(aminoAcidFragments, heterozygousLoci[i], heterozygousLoci[i + 1]
            )
            if (evidence.evidence.isNotEmpty()) {
                result.add(evidence)
            }
        }

        return result.sorted().filter { it.totalEvidence() > minTotalEvidence }
    }

    fun extendConsecutive(current: PhasedEvidence, others: Set<PhasedEvidence>): List<PhasedEvidence> {

//        val expected = setOf(1, 3)
//        if (expected.all { current.aminoAcidIndices.contains(it) }) {
//            println("HERE")
//        }

        val existingIndices = current.aminoAcidIndices
        val remainingIndices = heterozygousLoci.filter { it !in existingIndices }

        val minExisting = existingIndices.min()!!
        val maxExisting = existingIndices.max()!!

        val remainingIndicesAbove = remainingIndices.filter { it > maxExisting }.sorted()
        val remainingIndicesBelow = remainingIndices.filter { it < minExisting }.sorted().reversed()

        val result = mutableListOf<PhasedEvidence>()
        if (remainingIndicesAbove.isNotEmpty()) {
            val unambiguousIndices = current.unambiguousTailIndices() + remainingIndicesAbove[0]
            val allNewIndices = current.aminoAcidIndices + remainingIndicesAbove[0]
            val next = next(current, unambiguousIndices, allNewIndices, others)
            if (next != null) {
                result.add(next)
            }
        }

        if (remainingIndicesBelow.isNotEmpty()) {
            val unambiguousIndices = (current.unambiguousHeadIndices() + remainingIndicesBelow[0]).sortedArray()
            val allNewIndices = (current.aminoAcidIndices + remainingIndicesBelow[0]).sortedArray()
            val next = next(current, unambiguousIndices, allNewIndices, others)
            if (next != null) {
                result.add(next)
            }
        }

        return result.sorted().filter { it.totalEvidence() >= minTotalEvidence }
    }

    private fun next(current: PhasedEvidence, unambiguousIndices: IntArray, allIndices: IntArray, others: Set<PhasedEvidence>): PhasedEvidence? {
        val fake = PhasedEvidence(allIndices, Collections.emptyMap())
        if (!others.contains(fake)) {
            val newEvidence = PhasedEvidence.evidence(aminoAcidFragments, *unambiguousIndices)
            if (newEvidence.evidence.isNotEmpty()) {
                if (newEvidence.aminoAcidIndices.size == current.aminoAcidIndices.size + 1) {
                    return newEvidence
                } else {
                    val combinedEvidence = CombineEvidence.combineOverlapping(current, newEvidence)
                    return combinedEvidence
                }
            }
        }
        return null
    }


}