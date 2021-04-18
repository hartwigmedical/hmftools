package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.read.SAMRecordReader
import org.apache.logging.log4j.LogManager

data class BamQC(val discardedAlignmentFragments: Int, val discardedIndelFragments: Int, val discardedIndelMaxCount: Int, val discardedPonIndelFragments: Int, val discardedPonIndelMaxCount: Int) {

    companion object {
        private const val MIN_SUPPORT = 3
        private val logger = LogManager.getLogger(this::class.java)

        fun create(reader: SAMRecordReader): BamQC {
            val fragmentsWithUnmatchedPonIndel = reader.unmatchedPonIndels(MIN_SUPPORT)
            val fragmentsWithUnmatchedIndel = reader.unmatchedIndels(MIN_SUPPORT)
            for ((indel, count) in fragmentsWithUnmatchedIndel) {
                logger.warn("    UNMATCHED_INDEL - $count fragments excluded with unmatched indel $indel")
            }

            for ((indel, count) in fragmentsWithUnmatchedPonIndel) {
                logger.info("    UNMATCHED_PON_INDEL - $count fragments excluded with unmatched PON indel $indel")
            }

            return BamQC(reader.alignmentFiltered(),
                    fragmentsWithUnmatchedIndel.size, fragmentsWithUnmatchedIndel.values.max() ?: 0,
                    fragmentsWithUnmatchedPonIndel.size, fragmentsWithUnmatchedPonIndel.values.max() ?: 0)
        }
    }

    fun header(): List<String> {
        return listOf("discardedIndelFragments", "discardedIndelMaxCount", "discardedAlignmentFragments")
    }

    fun body(): List<String> {
        return listOf(discardedIndelFragments.toString(), discardedIndelMaxCount.toString(), discardedAlignmentFragments.toString())
    }
}

