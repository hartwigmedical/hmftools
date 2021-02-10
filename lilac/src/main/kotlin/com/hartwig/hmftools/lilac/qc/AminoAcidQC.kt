package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager
import kotlin.math.max

data class AminoAcidQC(val unusedAminoAcids: Int, val unusedAminoAcidMaxSupport: Int) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        fun create(winners: List<HlaSequenceLoci>, aminoAcidCount: SequenceCount): AminoAcidQC {

            var unused = 0
            var largest = 0

            for (locus in aminoAcidCount.heterozygousLoci()) {
                val expected = aminoAcidCount[locus]
                val actual = winners.filter { locus < it.sequences.size }.map { it.sequence(locus) }.toSet()
                for ((sequence, count) in expected) {
                    if (sequence !in actual) {
                        unused++
                        largest = max(largest, count)
                        logger.warn("UNMATCHED_AMINO_ACID - amino acid '$sequence' at locus $locus not in winning solution")
                    }
                }
            }

            return AminoAcidQC(unused, largest)
        }
    }

}