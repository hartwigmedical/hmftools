package com.hartwig.hmftools.lilackt.evidence

import com.hartwig.hmftools.lilackt.amino.AminoAcidFragment
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci
import java.util.*


data class PhasedEvidence(val aminoAcidIndices: IntArray, val evidence: Map<String, Int>) : Comparable<PhasedEvidence> {

    fun inconsistentEvidence(candidates: Collection<HlaSequenceLoci>): PhasedEvidence {
        fun consistentWithAny(sequence: String): Boolean {
            return candidates.any { it.consistentWith(sequence, *aminoAcidIndices) }
        }
        return PhasedEvidence(aminoAcidIndices, evidence.filter { !consistentWithAny(it.key) })
    }

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


    fun minEvidence(): Int {
        return evidence.values.min() ?: 0
    }

    fun totalEvidence(): Int {
        return evidence.values.sum()
    }

    fun removeSingles(minTotalEvidence: Int): PhasedEvidence {
        if (totalEvidence() >= minTotalEvidence) {
            return PhasedEvidence(aminoAcidIndices, evidence.filter { it.value > 1 })
        }

        return this
    }

    companion object {

        fun evidence(aminoAcidFragments: List<AminoAcidFragment>, vararg indices: Int): PhasedEvidence {
            val filteredFragments = aminoAcidFragments.filter { it.containsAll(indices) }
            val aminoAcidEvidence = filteredFragments.map { it.toAminoAcids(indices) }.groupingBy { it }.eachCount()
            return PhasedEvidence(indices, aminoAcidEvidence)
        }

        private fun AminoAcidFragment.toAminoAcids(indices: IntArray): String {
            return indices.map { this.aminoAcid(it) }.joinToString("")
        }

        private fun AminoAcidFragment.containsAll(indices: IntArray): Boolean {
            return indices.all { this.containsAminoAcid(it) }
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
