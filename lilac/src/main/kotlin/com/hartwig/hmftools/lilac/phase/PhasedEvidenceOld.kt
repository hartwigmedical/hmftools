package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.Fragment


data class PhasedEvidenceOld(val indices: IntArray, val aminoAcids: CharArray, val evidence: Int) {

    companion object {

        fun evidence(minQual: Int, fragments: List<Fragment>, vararg indices: Int): List<PhasedEvidenceOld> {
            val filteredFragments = fragments.filter { it.containsAll(minQual, indices) }
            val aminoAcidEvidence = filteredFragments.map { it.toAminoAcids(minQual, indices) }.groupingBy { it }.eachCount()
            return aminoAcidEvidence.entries.map { (aminoAcids, count) -> PhasedEvidenceOld(indices, aminoAcids.toCharArray(), count) }
        }

        private fun Fragment.toAminoAcids(minQual: Int, indices: IntArray): String {
            return indices.map { this.aminoAcid(it, minQual) }.joinToString("")
        }

        private fun Fragment.containsAll(minQual: Int, indices: IntArray): Boolean {
            return indices.all { this.containsAminoAcid(it) && this.aminoAcid(it, minQual) != '.' }
        }
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is PhasedEvidenceOld) return false

        if (!indices.contentEquals(other.indices)) return false
        if (!aminoAcids.contentEquals(other.aminoAcids)) return false
        if (evidence != other.evidence) return false

        return true
    }

    override fun hashCode(): Int {
        var result = indices.contentHashCode()
        result = 31 * result + aminoAcids.contentHashCode()
        result = 31 * result + evidence
        return result
    }

    override fun toString(): String {
        return "PhasedEvidence(indices=${indices.contentToString()}, aminoAcids=${aminoAcids.contentToString()}, evidence=$evidence)"
    }


}
