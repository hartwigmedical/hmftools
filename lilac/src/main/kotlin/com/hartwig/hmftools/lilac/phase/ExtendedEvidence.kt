package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.Fragment
import java.util.*

class ExtendedEvidence(private val minMinEvidence: Int, private val minTotalEvidence: Int, private val heterozygousIndices: List<Int>, private val fragments: List<Fragment>) {

    fun initialEvidence(): List<PhasedEvidence> {
        val result = mutableListOf<PhasedEvidence>()

        for (i in 0..(heterozygousIndices.size - 2)) {
            val evidence = PhasedEvidence.evidence(fragments, heterozygousIndices[i], heterozygousIndices[i + 1]
            )
            if (evidence.evidence.isNotEmpty()) {
                result.add(evidence)
            }
        }

        return result.sorted().filter { it.totalEvidence() > minTotalEvidence }
    }

    fun extendConsecutive(current: PhasedEvidence, others: Set<PhasedEvidence>): List<PhasedEvidence> {

//        val expected = setOf(117, 118, 120)
//        if (expected.all { current.aminoAcidIndices.contains(it) }) {
//            println("HERE")
//        }

        val existingIndices = current.aminoAcidIndices
        val remainingIndices = heterozygousIndices.filter { it !in existingIndices }

        val minExisting = existingIndices.min()!!
        val maxExisting = existingIndices.max()!!

        val remainingIndicesAbove = remainingIndices.filter { it > maxExisting }.sorted()
        val remainingIndicesBelow = remainingIndices.filter { it < minExisting }.sorted().reversed()

        val result = mutableListOf<PhasedEvidence>()
        if (remainingIndicesAbove.isNotEmpty()) {
            val unambiguousIndices = current.unambiguousTailIndices() + remainingIndicesAbove[0]
            val allNewIndices = current.aminoAcidIndices + remainingIndicesAbove[0]
            val fake = PhasedEvidence(allNewIndices, Collections.emptyMap())
            if (!others.contains(fake)) {
                val newEvidence = PhasedEvidence.evidence(fragments, *unambiguousIndices)
                if (newEvidence.evidence.isNotEmpty() && newEvidence.minEvidence() >= minMinEvidence) {
                    if (newEvidence.aminoAcidIndices.size == current.aminoAcidIndices.size + 1) {
                        result.add(newEvidence)
                    } else {
                        val combinedEvidence = PhasedEvidence.combineOverlapping(current, newEvidence)
                        result.add(combinedEvidence)
                    }
                }
            }
        }

        if (remainingIndicesBelow.isNotEmpty()) {
            val unambiguousIndices = (current.unambiguousHeadIndices() + remainingIndicesBelow[0]).sortedArray()
            val allNewIndices = (current.aminoAcidIndices + remainingIndicesBelow[0]).sortedArray()
            val fake = PhasedEvidence(allNewIndices, Collections.emptyMap())
            if (!others.contains(fake)) {
                val newEvidence = PhasedEvidence.evidence(fragments, *unambiguousIndices)
                if (newEvidence.evidence.isNotEmpty() && newEvidence.minEvidence() >= minMinEvidence) {
                    if (newEvidence.aminoAcidIndices.size == current.aminoAcidIndices.size + 1) {
                        result.add(newEvidence)
                    } else {
                        val combinedEvidence = PhasedEvidence.combineOverlapping(newEvidence, current)
                        result.add(combinedEvidence)
                    }
                }
            }
        }

        return result.sorted().filter { it.totalEvidence() >= minTotalEvidence }
    }

}