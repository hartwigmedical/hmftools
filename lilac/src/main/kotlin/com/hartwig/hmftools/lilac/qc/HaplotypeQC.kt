package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.SequenceCount
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

        fun create(minEvidence: Int, winners: List<HlaSequenceLoci>, evidence: List<PhasedEvidence>, aminoAcidCount: SequenceCount): HaplotypeQC {
            val unused = evidence
                    .flatMap { it.unmatchedHaplotype(minEvidence, winners, aminoAcidCount) }
                    .distinct()
                    .sortedBy { it.startLoci}

            var unusedCount = 0
            var maxSupport = 0
            var maxLength = 0
            for (unusedEvidence in unused) {
                maxSupport = max(maxSupport, unusedEvidence.supportingFragments)
                maxLength = max(maxLength, unusedEvidence.haplotype.length)
                unusedCount++

                logger.warn("UNMATCHED_HAPLTOYPE - $unusedEvidence")

            }

            return HaplotypeQC(unusedCount, maxSupport, maxLength)
        }

        fun PhasedEvidence.unmatchedHaplotype(minEvidence: Int, winners: Collection<HlaSequenceLoci>, aminoAcidCount: SequenceCount): List<UnmatchedHaplotype> {
            fun consistentWithAny(sequence: String): Boolean {
                return winners.any { it.consistentWith(sequence, *aminoAcidIndices) }
            }

            val unmatched = evidence
                    .filter { !consistentWithAny(it.key) }
                    .filter { it.value >= minEvidence }

            if (unmatched.isEmpty()) {
                return listOf()
            }

            return unmatched.map { Pair(it.key, it.value) }.map { UnmatchedHaplotype.create(this.aminoAcidIndices, it, aminoAcidCount) }
        }
    }
}


data class UnmatchedHaplotype(val startLoci: Int, val endLoci: Int, val supportingFragments: Int, val haplotype: String) {
    companion object {
        fun create(aminoAcidIndices: IntArray, evidence: Pair<String, Int>, aminoAcidCount: SequenceCount): UnmatchedHaplotype {
            require(aminoAcidIndices.isNotEmpty())
            val startLoci = aminoAcidIndices.min()!!
            val endLoci = aminoAcidIndices.max()!!
            val sparseHaplotype = evidence.first
            val completeHaplotype = (startLoci..endLoci).joinToString("") { if (aminoAcidIndices.contains(it)) sparseHaplotype[aminoAcidIndices.indexOf(it)].toString() else aminoAcidCount.sequenceAt(it).first() }
            return UnmatchedHaplotype(startLoci, endLoci, evidence.second, completeHaplotype)
        }
    }

    override fun toString(): String {
        return "startLoci=$startLoci, endLoci=$endLoci, supportingFragments=$supportingFragments, haplotype='$haplotype'"
    }


}