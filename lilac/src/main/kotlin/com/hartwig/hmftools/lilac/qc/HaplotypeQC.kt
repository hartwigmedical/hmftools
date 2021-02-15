package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.evidence.PhasedEvidence
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager
import kotlin.math.max

data class HaplotypeQC(val unusedHaplotypes: Int, val unusedHaplotypeMaxSupport: Int, val unusedHaplotypeMaxLength: Int) {

    fun header(): List<String> {
        return listOf("unusedHaplotypes", "unusedHaplotypeMaxSupport", "unusedHaplotypeMaxLength")
    }

    fun body(): List<String> {
        return listOf(unusedHaplotypes.toString(), unusedHaplotypeMaxSupport.toString(), unusedHaplotypeMaxLength.toString())
    }

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        fun create(minEvidence: Int, winners: List<HlaSequenceLoci>, evidence: List<PhasedEvidence>): HaplotypeQC {
            val unused = evidence
                    .map { it.unmatchedHaplotype(minEvidence, winners) }
                    .filter { it.haplotypes.isNotEmpty() }
                    .distinct()
                    .sortedBy { it.loci.min() }

            var unusedCount = 0
            var maxSupport = 0
            var maxLength = 0
            for (unusedEvidence in unused) {
                for ((haplotype, count) in unusedEvidence.haplotypes) {
                    maxSupport = max(maxSupport, count)
                    maxLength = max(maxLength, haplotype.length)
                    unusedCount++
                }

                logger.warn("UNMATCHED_HAPLTOYPE - $unusedEvidence")

            }

            return HaplotypeQC(unusedCount, maxSupport, maxLength)
        }

        fun PhasedEvidence.unmatchedHaplotype(minEvidence: Int, winners: Collection<HlaSequenceLoci>): UnmatchedHaplotype {
            fun consistentWithAny(sequence: String): Boolean {
                return winners.any { it.consistentWith(sequence, *aminoAcidIndices) }
            }

            val unmatched = evidence
                    .filter { !consistentWithAny(it.key) }
                    .filter { it.value >= minEvidence }

            return UnmatchedHaplotype(aminoAcidIndices, unmatched)
        }
    }
}


data class UnmatchedHaplotype(val loci: IntArray, val haplotypes: Map<String, Int>) {

    fun totalEvidence(): Int {
        return haplotypes.values.sum()
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is UnmatchedHaplotype) return false

        if (!loci.contentEquals(other.loci)) return false
        if (haplotypes != other.haplotypes) return false

        return true
    }

    override fun hashCode(): Int {
        var result = loci.contentHashCode()
        result = 31 * result + haplotypes.hashCode()
        return result
    }

    override fun toString(): String {
        return "loci=${loci.contentToString()}, haplotypes=$haplotypes)"
    }
}