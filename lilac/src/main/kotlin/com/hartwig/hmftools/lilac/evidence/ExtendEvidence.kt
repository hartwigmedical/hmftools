package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles

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
//                    println("FAIL:$combinedEvidence")
                }
            }
        }

        return result.sorted()
    }


    fun merge(current: PhasedEvidence, others: Set<PhasedEvidence>): Pair<PhasedEvidence, Set<PhasedEvidence>> {

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


}