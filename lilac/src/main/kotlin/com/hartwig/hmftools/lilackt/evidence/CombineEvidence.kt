package com.hartwig.hmftools.lilackt.evidence

import kotlin.math.min

object CombineEvidence {

    fun combine(left: PhasedEvidence, common: PhasedEvidence, right: PhasedEvidence): PhasedEvidence {
        val leftWithCommon = combineOverlapping(left, common)
        val result =  combineOverlapping(leftWithCommon, right)
        return result
    }

    fun canCombine(left: PhasedEvidence, common: PhasedEvidence, right: PhasedEvidence): Boolean {
        return canCombine(left, common) && canCombine(common, right)
    }

    fun canCombine(left: PhasedEvidence, right: PhasedEvidence): Boolean {
        val indexIntersection = (left.aminoAcidIndices.toSet() intersect right.aminoAcidIndices.toSet()).toList().sorted()
        if (indexIntersection.isEmpty()) {
            return false
        }

        val uniqueToLeft = left.aminoAcidIndices.filter { it !in indexIntersection }
        val uniqueToRight = right.aminoAcidIndices.filter { it !in indexIntersection }

        if (uniqueToRight.isNotEmpty() && uniqueToLeft.any { it > uniqueToRight.min()!! }) {
            return false
        }

        if (uniqueToLeft.isNotEmpty() && uniqueToRight.any { it < uniqueToLeft.max()!! }) {
            return false
        }

        val leftCommon = left.evidence.keys.map { it.substring(left.aminoAcidIndices.size - indexIntersection.size) }.toSet()
        val rightCommon = right.evidence.keys.map { it.substring(0, indexIntersection.size) }.toSet()

        return leftCommon == rightCommon
    }


    fun combineOverlapping(left: PhasedEvidence, right: PhasedEvidence): PhasedEvidence {

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