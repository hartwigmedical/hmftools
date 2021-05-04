package com.hartwig.hmftools.lilackt.qc

import com.hartwig.hmftools.lilackt.SequenceCount
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager
import kotlin.math.max

data class AminoAcidQC(val unusedAminoAcids: Int, val unusedAminoAcidMaxSupport: Int) {

    fun header(): List<String> {
        return listOf("unusedAminoAcids", "unusedAminoAcidMaxSupport")
    }

    fun body(): List<String> {
        return listOf(unusedAminoAcids.toString(), unusedAminoAcidMaxSupport.toString())
    }

    companion object {
        private const val MIN_COUNT = 3
        private val logger = LogManager.getLogger(this::class.java)

        fun create(winners: Set<HlaSequenceLoci>, aminoAcidCount: SequenceCount): AminoAcidQC {

            var unused = 0
            var largest = 0

            for (locus in aminoAcidCount.heterozygousLoci()) {
                val expected = aminoAcidCount[locus]
                val actual = winners.filter { locus < it.sequences.size }.map { it.sequence(locus) }.toSet()
                for ((sequence, count) in expected) {
                    if (count >= MIN_COUNT && sequence !in actual) {
                        unused++
                        largest = max(largest, count)
                        logger.warn("    UNMATCHED_AMINO_ACID - amino acid '$sequence' with $count support at locus $locus not in winning solution")
                    }
                }
            }

            return AminoAcidQC(unused, largest)
        }
    }

}