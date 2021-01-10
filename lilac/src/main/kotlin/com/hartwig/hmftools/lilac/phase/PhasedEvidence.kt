package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.Fragment


data class PhasedEvidence(val aminoAcidIndices: IntArray, val evidence: Map<String, Int>, val fragments: List<Fragment>) : Comparable<PhasedEvidence> {


    fun minOverlapEnd(): Int {
        for (i in 1..aminoAcidIndices.size) {
            evidence.map { it.key }.map { it.substring(it.length - 1 - i, it.length ) }.groupingBy {  }.eachCount()


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
            return PhasedEvidence(indices, aminoAcidEvidence, filteredFragments)
        }

        private fun Fragment.toAminoAcids(minQual: Int, indices: IntArray): String {
            return indices.map { this.aminoAcid(it, minQual) }.joinToString("")
        }

        private fun Fragment.containsAll(minQual: Int, indices: IntArray): Boolean {
            return indices.all { this.containsAminoAcid(it) && this.aminoAcid(it, minQual) != '.' }
        }
    }


    override fun toString(): String {
        return "PhasedEvidence(indices=${aminoAcidIndices.contentToString()}, evidence=$evidence, total=${totalEvidence()})"
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
