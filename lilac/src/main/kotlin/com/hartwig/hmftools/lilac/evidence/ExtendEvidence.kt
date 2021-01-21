package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import java.util.*

class ExtendEvidence(private val minEvidenceTotal: Int, private val heterozygousLoci: List<Int>, private val aminoAcidFragments: List<AminoAcidFragment>) {

    companion object {
        const val MIN_EVIDENCE_SEQUENCE = 1
    }

    fun initialEvidence(): List<PhasedEvidence> {
        val result = mutableListOf<PhasedEvidence>()

        for (i in 0..(heterozygousLoci.size - 2)) {
            val evidence = PhasedEvidence.evidence(MIN_EVIDENCE_SEQUENCE, aminoAcidFragments, heterozygousLoci[i], heterozygousLoci[i + 1]
            )
            if (evidence.evidence.isNotEmpty()) {
                result.add(evidence)
            }
        }

        return result.sorted().filter { it.totalEvidence() > minEvidenceTotal }
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
            val next = next(true, current, unambiguousIndices, allNewIndices, others)
            if (next != null) {
                result.add(next)
            }
        }

        if (remainingIndicesBelow.isNotEmpty()) {
            val unambiguousIndices = (current.unambiguousHeadIndices() + remainingIndicesBelow[0]).sortedArray()
            val allNewIndices = (current.aminoAcidIndices + remainingIndicesBelow[0]).sortedArray()
            val next = next(false, current, unambiguousIndices, allNewIndices, others)
            if (next != null) {
                result.add(next)
            }
        }

        return result.sorted().filter { it.totalEvidence() >= minEvidenceTotal }
    }

    private fun next(currentIsLeft: Boolean, current: PhasedEvidence, unambiguousIndices: IntArray, allIndices: IntArray, others: Set<PhasedEvidence>): PhasedEvidence? {
        val fake = PhasedEvidence(allIndices, Collections.emptyMap())
        if (!others.contains(fake)) {
            val newEvidence = PhasedEvidence.evidence(MIN_EVIDENCE_SEQUENCE, aminoAcidFragments, *unambiguousIndices)
            if (newEvidence.totalEvidence() < minEvidenceTotal) {
                return null
            }

            if (newEvidence.evidence.isNotEmpty()) {
                val allIndicesInNewEvidence = newEvidence.aminoAcidIndices.size == allIndices.size
                val left = if (currentIsLeft) current else newEvidence
                val right = if (currentIsLeft) newEvidence else current
                if (!CombineEvidence.canCombine(left, right)) {
//                    println("BAD MERGE: $allIndicesInNewEvidence")
//                    println(left)
//                    println(right)
                    return null
                }

                return if (allIndicesInNewEvidence) {
                    newEvidence
                } else {
                    CombineEvidence.combineOverlapping(left, right)
                }
            }
        }
        return null
    }


}