package com.hartwig.hmftools.lilac.evidence

import kotlin.math.min

object CombineEvidence {

    fun canCombine(left: PhasedEvidence, right: PhasedEvidence): Boolean {

        val indexIntersection = (left.aminoAcidIndices.toSet() intersect right.aminoAcidIndices.toSet()).toList().sorted()
        if (indexIntersection.isEmpty()) {
            return false
        }

        val uniqueToLeft = left.aminoAcidIndices.filter { it !in indexIntersection }
        if (uniqueToLeft.isEmpty()) {
            return false
        }

        val uniqueToRight = right.aminoAcidIndices.filter { it !in indexIntersection }
        if (uniqueToLeft.isEmpty()) {
            return false
        }

        val maxUniqueToLeft = uniqueToLeft.max()!!
        val minUniqueToRight = uniqueToRight.min()!!
        if (maxUniqueToLeft > minUniqueToRight) {
            return false
        }

        val leftCommon = left.evidence.keys.map { it.substring(uniqueToLeft.size) }.toSet()
        val rightCommon = right.evidence.keys.map { it.substring(0, indexIntersection.size) }.toSet()

        return leftCommon == rightCommon
    }


    fun combineOverlapping(left: PhasedEvidence, right: PhasedEvidence): PhasedEvidence {
        assert(canCombine(left, right))

        val indexUnion = (left.aminoAcidIndices.toSet() union right.aminoAcidIndices.toSet()).toList().sorted()
        val indexIntersection = (left.aminoAcidIndices.toSet() intersect right.aminoAcidIndices.toSet()).toList().sorted()
        val uniqueToLeft = left.aminoAcidIndices.filter { it !in indexIntersection }

        val result = mutableMapOf<String, Int>()
        for ((leftSequence, leftCount) in left.evidence) {
            val leftUnique = leftSequence.substring(0, uniqueToLeft.size)
            val leftCommon = leftSequence.substring(uniqueToLeft.size)
            for ((rightSequence, rightCount) in right.evidence) {
                val rightCommon = rightSequence.substring(0, indexIntersection.size)
                if (leftCommon == rightCommon) {
                    val combined = leftUnique + rightSequence
                    result[combined] = min(leftCount, rightCount)
                }
            }
        }
        return PhasedEvidence(indexUnion.toIntArray(), result)
    }

}