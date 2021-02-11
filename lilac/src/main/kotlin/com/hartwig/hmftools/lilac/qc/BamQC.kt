package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.read.SAMRecordReader

data class BamQC(val discardedIndelFragments: Int, val discardedIndelMaxCount: Int) {

    companion object {
        fun create(reader: SAMRecordReader): BamQC {
            val fragmentsWithUnmatchedIndel = reader.unmatchedIndels(2)
            for ((indel, count) in fragmentsWithUnmatchedIndel) {
                LilacQC.logger.warn("UNMATCHED_INDEL - $count fragments excluded with unmatched indel $indel")
            }
            return BamQC(fragmentsWithUnmatchedIndel.size, fragmentsWithUnmatchedIndel.values.max() ?: 0)
        }
    }

}

