package com.hartwig.hmftools.lilackt.nuc

import kotlin.math.max
import kotlin.math.min

class ExpectedAlleles(private val alleleCount: IntArray) {

    fun expectedAlleles(loci: Collection<Int>): Int {
        return loci.map { expectedAlleles(it) }.min()!!
    }

    fun expectedAlleles(loci: Int): Int {
        return if (loci >= alleleCount.size) 2 else alleleCount[loci]
    }

    companion object {
        fun expectedAlleles(otherMin1: Int, otherMin2: Int): ExpectedAlleles {
            val min = min(otherMin1, otherMin2)
            val max = max(otherMin1, otherMin2)

            val three = (1 until min).map { 6 }
            val two = (min until max).map { 4 }
            return ExpectedAlleles((three + two).toIntArray())

        }
    }
}


