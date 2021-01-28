package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles
import java.util.*

class ExtendEvidence(
        private val minFragmentsPerAllele: Int,
        private val minFragmentsToRemoveSingles: Int,
        private val heterozygousLoci: List<Int>,
        private val aminoAcidFragments: List<AminoAcidFragment>,
        private val expectedAlleles: ExpectedAlleles) {

    fun pairedEvidence(): List<PhasedEvidence> {
        val result = mutableListOf<PhasedEvidence>()

        for (i in 0..(heterozygousLoci.size - 2)) {
            val indices = listOf(heterozygousLoci[i], heterozygousLoci[i + 1])
            val filteredFragments = aminoAcidFragments.filter { it.containsAll(indices) }
            if (filteredFragments.isNotEmpty()) {
                val minTotalFragments = minTotalFragments(indices)

                val left = PhasedEvidence.evidence(aminoAcidFragments, indices[0])
                val right = PhasedEvidence.evidence(aminoAcidFragments, indices[1])
                val combinedEvidence = PhasedEvidence.evidence(aminoAcidFragments, *indices.toIntArray()).removeSingles(minFragmentsToRemoveSingles)
                if (combinedEvidence.totalEvidence() >= minTotalFragments && CombineEvidence.canCombine(left, combinedEvidence, right)) {
                    result.add(combinedEvidence)
                } else {
                    println("FAIL:" + combinedEvidence)
                }
            }
        }

        return result.sorted()
    }


    fun merge(current: PhasedEvidence, others: Set<PhasedEvidence>): Pair<PhasedEvidence, Set<PhasedEvidence>> {

//        val expected = setOf(206, 207, 212, 217, 242)
//        if (current.aminoAcidIndices.size == expected.size && expected.all { current.aminoAcidIndices.contains(it) }) {
//            println("HERE")
//        }

        val existingIndices = current.aminoAcidIndices
        val minExisting = existingIndices.min()!!
        val maxExisting = existingIndices.max()!!

        val othersContainingMax = others.filter { it != current && it.aminoAcidIndices.contains(maxExisting) }
        val othersContainingMin = others.filter { it != current && it.aminoAcidIndices.contains(minExisting) }

        if (othersContainingMin.isNotEmpty() && othersContainingMax.isNotEmpty()) {
            val left = othersContainingMin[0]
            val right = othersContainingMax[0]
            return if (left.totalEvidence() > right.totalEvidence()) {
                val result = merge(current, left, current)
                if (result.second.isEmpty()) {
                    merge(current, current, right)
                } else {
                    result
                }
            } else {
                val result = merge(current, current, right)
                if (result.second.isEmpty()) {
                    merge(current, left, current)
                } else {
                    result
                }
            }
        }

        if (othersContainingMin.isNotEmpty()) {
            val left = othersContainingMin[0]
            return merge(current, left, current)
        }

        if (othersContainingMax.isNotEmpty()) {
            val right = othersContainingMax[0]
            return merge(current, current, right)
        }

        return Pair(current, setOf())
    }


    private fun minTotalFragments(indices: List<Int>): Int {
        return expectedAlleles.expectedAlleles(indices) * minFragmentsPerAllele
    }

    private fun merge(current: PhasedEvidence, left: PhasedEvidence, right: PhasedEvidence): Pair<PhasedEvidence, Set<PhasedEvidence>> {

        val leftTail = left.unambiguousTailIndices()
        val rightHead = right.unambiguousHeadIndices()
        val mergeIndices = (leftTail + rightHead).distinct().sorted()

        val filteredFragments = aminoAcidFragments.filter { it.containsAll(mergeIndices) }
        if (filteredFragments.isNotEmpty()) {
            val minTotalFragments = minTotalFragments(mergeIndices)
            val mergeEvidence = PhasedEvidence.evidence(filteredFragments, *mergeIndices.toIntArray()).removeSingles(minFragmentsToRemoveSingles)
            if (CombineEvidence.canCombine(left, mergeEvidence, right)) {
                val combined = CombineEvidence.combine(left, mergeEvidence, right)
                if (combined.totalEvidence() >= minTotalFragments) {
                    return Pair(combined, setOf(left, right))
                }
            }
        }

        return Pair(current, setOf())
    }

    fun extendConsecutive(current: PhasedEvidence, others: Set<PhasedEvidence>): List<PhasedEvidence> {


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

        return result.sorted().filter { it.totalEvidence() >= 30 }
    }

    private fun next(currentIsLeft: Boolean, current: PhasedEvidence, unambiguousIndices: IntArray, allIndices: IntArray, others: Set<PhasedEvidence>): PhasedEvidence? {
        val fake = PhasedEvidence(allIndices, Collections.emptyMap())
        if (!others.contains(fake)) {

            val newEvidence = PhasedEvidence.evidence(aminoAcidFragments, *unambiguousIndices)
            if (newEvidence.totalEvidence() < 30) {
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