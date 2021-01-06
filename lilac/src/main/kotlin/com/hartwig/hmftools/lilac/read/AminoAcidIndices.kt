package com.hartwig.hmftools.lilac.read

object AminoAcidIndices {

    fun indices(startIndex: Int, endIndex: Int): IntRange {
        val start = startIndex / 3 + (if (startIndex % 3 == 0) 0 else 1)
        val end = (endIndex + 1) / 3 - 1

        return IntRange(start, end)
    }

}