package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.Fragment
import java.util.*
import kotlin.math.min


data class PhasedEvidence(val aminoAcidIndices: IntArray, val evidence: Map<String, Int>) : Comparable<PhasedEvidence> {

    fun unambiguousHeadIndices(): IntArray {
        return aminoAcidIndices.filterIndexed { index, _ -> index < unambiguousHeadLength() }.toIntArray()
    }

    fun unambiguousTailIndices(): IntArray {
        val minIndex = aminoAcidIndices.size - unambiguousTailLength()
        return aminoAcidIndices.filterIndexed { index, _ -> index >= minIndex }.toIntArray()
    }

    private fun unambiguousTailLength(): Int {
        for (i in aminoAcidIndices.indices) {
            val endIndex = aminoAcidIndices.size
            val startIndex = endIndex - i - 1
            val evidenceTails = evidence.keys.map { it.substring(startIndex, endIndex) }.toSet()
            if (evidenceTails.size == evidence.size) {
                return i + 1
            }

        }

        return aminoAcidIndices.size
    }

    private fun unambiguousHeadLength(): Int {
        for (i in aminoAcidIndices.indices) {
            val length = i + 1
            val evidenceTails = evidence.keys.map { it.substring(0, length) }.toSet()
            if (evidenceTails.size == evidence.size) {
                return length
            }

        }

        return aminoAcidIndices.size
    }

    fun contains(other: PhasedEvidence): Boolean {
        val overlap = (other.aminoAcidIndices.toSet() intersect aminoAcidIndices.toSet()).size
        return overlap == other.aminoAcidIndices.size
    }

    fun overlaps(other: PhasedEvidence): Boolean {
        val overlap = (other.aminoAcidIndices.toSet() intersect aminoAcidIndices.toSet()).size
        return overlap > 0 && overlap < aminoAcidIndices.size && overlap < other.aminoAcidIndices.size
    }

    fun minEvidence(): Int {
        return evidence.values.min() ?: 0
    }

    fun totalEvidence(): Int {
        return evidence.values.sum()
    }


    companion object {

        fun evidence(fragments: List<Fragment>, vararg indices: Int): PhasedEvidence {
            val filteredFragments = fragments.filter { it.containsAll(indices) }
            val aminoAcidEvidence = filteredFragments.map { it.toAminoAcids(indices) }.groupingBy { it }.eachCount()
            return PhasedEvidence(indices, aminoAcidEvidence)
        }

        private fun Fragment.toAminoAcids(indices: IntArray): String {
            return indices.map { this.aminoAcid(it) }.joinToString("")
        }

        private fun Fragment.containsAll(indices: IntArray): Boolean {
            return indices.all { this.containsAminoAcid(it) && this.aminoAcid(it) != '.' }
        }

        fun combineOverlapping(left: PhasedEvidence, right: PhasedEvidence): PhasedEvidence {
            assert(left.overlaps(right))
            assert(right.overlaps(left))
            assert(left.aminoAcidIndices[0] < right.aminoAcidIndices[0])

            val indexUnion = (left.aminoAcidIndices.toSet() union right.aminoAcidIndices.toSet()).toList().sorted()
            val indexIntersection = (left.aminoAcidIndices.toSet() intersect right.aminoAcidIndices.toSet()).toList().sorted()

            val uniqueToLeft = left.aminoAcidIndices.filter { it !in indexIntersection }
            val uniqueToRight = right.aminoAcidIndices.filter { it !in indexIntersection }

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

    fun evidenceString(): String {
        val resultBuilder = StringJoiner(" ")
        evidence.forEach { (evidence, count) -> resultBuilder.add(toEvidenceString(evidence, count)) }
        return "{$resultBuilder}"
    }

    fun toEvidenceString(evidence: String, count: Int): String {
        val resultBuilder = StringJoiner("")
        var evidenceIndex = 0
        for (i in aminoAcidIndices[0]..aminoAcidIndices[aminoAcidIndices.lastIndex]) {
            if (i in aminoAcidIndices) {
                resultBuilder.add(evidence[evidenceIndex].toString())
                evidenceIndex++
            } else {
                resultBuilder.add("-")
            }
        }

        return resultBuilder.add("=").add(count.toString()).toString()
    }

    override fun toString(): String {
        val uniqueTail = unambiguousTailLength()
        val uniqueHead = unambiguousHeadLength()
        val bases = 3 * (aminoAcidIndices[aminoAcidIndices.size - 1] - aminoAcidIndices[0])

//        return "PhasedEvidence(head=${uniqueHead} tail=${uniqueTail} loci=${aminoAcidIndices.size} types=${evidence.size} indices=${aminoAcidIndices.contentToString()}, evidence=${evidenceString()}, total=${totalEvidence()})"
        return "PhasedEvidence(head=${uniqueHead} tail=${uniqueTail} loci=${aminoAcidIndices.size} types=${evidence.size} indices=${aminoAcidIndices.contentToString()}, evidence=${evidence}, total=${totalEvidence()})"
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
