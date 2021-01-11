package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.Fragment
import kotlin.math.min


data class PhasedEvidence(val aminoAcidIndices: IntArray, val evidence: Map<String, Int>) : Comparable<PhasedEvidence> {

    private fun evidenceIndex(index: Int): Int {
        for (i in aminoAcidIndices.indices) {
            if (aminoAcidIndices[i] == index) {
                return i
            }
        }

        return -1
    }

    fun contains(other: PhasedEvidence): Boolean {
        val overlap = (other.aminoAcidIndices.toSet() intersect aminoAcidIndices.toSet()).size
        return overlap == other.aminoAcidIndices.size
    }

    fun overlaps(other: PhasedEvidence): Boolean {
        val overlap = (other.aminoAcidIndices.toSet() intersect aminoAcidIndices.toSet()).size
        return overlap > 0 && overlap < aminoAcidIndices.size && overlap < other.aminoAcidIndices.size
    }


    fun combine(other: PhasedEvidence): PhasedEvidence {
        if (this.contains(other)) {
            return this
        }

        if (other.contains(this)) {
            return other
        }



        val indexUnion = (other.aminoAcidIndices.toSet() union aminoAcidIndices.toSet()).toList().sorted()
        val indexIntersection = (other.aminoAcidIndices.toSet() intersect aminoAcidIndices.toSet()).toList().sorted()
        val evidenceIndices = indexIntersection.map { Pair(this.evidenceIndex(it), other.evidenceIndex(it)) }

        return this
    }



    fun minOverlapEnd(): Int {
        for (i in 1..aminoAcidIndices.size) {
            evidence.map { it.key }.map { it.substring(it.length - 1 - i, it.length) }.groupingBy { }.eachCount()


        }

        return aminoAcidIndices.size
    }


    fun minEvidence(): Int {
        return evidence.values.min() ?: 0
    }

    fun totalEvidence(): Int {
        return evidence.values.sum()
    }


    companion object {

        fun evidence(minQual: Int, fragments: List<Fragment>, vararg indices: Int): PhasedEvidence {
            val filteredFragments = fragments.filter { it.containsAll(minQual, indices) }
            val aminoAcidEvidence = filteredFragments.map { it.toAminoAcids(minQual, indices) }.groupingBy { it }.eachCount()
            return PhasedEvidence(indices, aminoAcidEvidence)
        }

        private fun Fragment.toAminoAcids(minQual: Int, indices: IntArray): String {
            return indices.map { this.aminoAcid(it, minQual) }.joinToString("")
        }

        private fun Fragment.containsAll(minQual: Int, indices: IntArray): Boolean {
            return indices.all { this.containsAminoAcid(it) && this.aminoAcid(it, minQual) != '.' }
        }

        fun combineOverlapping(left: PhasedEvidence, right: PhasedEvidence): PhasedEvidence {
            assert(left.overlaps(right))
            assert(right.overlaps(left))
            assert(left.aminoAcidIndices[0] < right.aminoAcidIndices[0])

            val indexUnion = (left.aminoAcidIndices.toSet() union right.aminoAcidIndices.toSet()).toList().sorted()
            val indexIntersection = (left.aminoAcidIndices.toSet() intersect right.aminoAcidIndices.toSet()).toList().sorted()

            val uniqueToLeft = left.aminoAcidIndices.filter { it !in indexIntersection}
            val uniqueToRight = right.aminoAcidIndices.filter { it !in indexIntersection}

            assert(uniqueToLeft.all { it < indexIntersection.min() ?: 0 })
            assert(uniqueToRight.all { it > indexIntersection.max() ?: 0 })

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


    override fun toString(): String {
        return "PhasedEvidence(bases=${ 3 * (aminoAcidIndices[aminoAcidIndices.size - 1] - aminoAcidIndices[0])} loci=${aminoAcidIndices.size} types=${evidence.size} indices=${aminoAcidIndices.contentToString()}, evidence=$evidence, total=${totalEvidence()})"
    }

    override fun compareTo(other: PhasedEvidence): Int {
        val totalCompare = -totalEvidence().compareTo(other.totalEvidence())
        if (totalCompare != 0) {
            return totalCompare
        }
        return minEvidence().compareTo(other.minEvidence())
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is PhasedEvidence) return false

        if (!aminoAcidIndices.contentEquals(other.aminoAcidIndices)) return false

        return true
    }

    override fun hashCode(): Int {
        return aminoAcidIndices.contentHashCode()
    }


}
