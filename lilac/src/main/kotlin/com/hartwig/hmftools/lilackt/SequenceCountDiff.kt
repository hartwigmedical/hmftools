package com.hartwig.hmftools.lilackt

import kotlin.math.min

data class SequenceCountDiff(val loci: Int, val sequence: String, val referenceCount: Int, val referenceDepth: Int, val tumorCount: Int, val tumorDepth: Int) {

    companion object {
        fun create(referenceCount: SequenceCount, tumorCount: SequenceCount): List<SequenceCountDiff> {
            val result = mutableListOf<SequenceCountDiff>()
            val minLength = min(referenceCount.length, tumorCount.length)

            for (loci in 0 until minLength) {
                val refDepth = referenceCount.depth(loci)
                val tumorDepth = tumorCount.depth(loci)

                val ref = referenceCount[loci]
                val tumor = tumorCount[loci]
                val sequences = ref.keys + tumor.keys
                for (sequence in sequences) {
                    if (sequence !in tumor || sequence !in ref) {
                        result.add(SequenceCountDiff(loci, sequence, ref[sequence] ?: 0, refDepth, tumor[sequence] ?: 0, tumorDepth))
                    }
                }
            }

            for (loci in minLength until tumorCount.length) {
                val tumorDepth = tumorCount.depth(loci)
                for ((sequence, count) in tumorCount[loci]) {
                    result.add(SequenceCountDiff(loci, sequence, 0, 0, count, tumorDepth))
                }
            }

            for (loci in minLength until referenceCount.length) {
                val refDepth = referenceCount.depth(loci)
                for ((sequence, count) in referenceCount[loci]) {
                    result.add(SequenceCountDiff(loci, sequence, count, refDepth, 0, 0))
                }
            }


            return result
        }
    }

}