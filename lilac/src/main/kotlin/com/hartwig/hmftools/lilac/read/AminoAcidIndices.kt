package com.hartwig.hmftools.lilac.read

object AminoAcidIndices {

    fun indices(nucStartIndex: Int, nucEndIndex: Int): IntRange {
        val start = nucStartIndex / 3 + (if (nucStartIndex % 3 == 0) 0 else 1)
        val end = (nucEndIndex + 1) / 3 - 1

        return IntRange(start, end)
    }

}