package com.hartwig.hmftools.lilac.prot

object ProteinSequenceBoundaries {

    fun aBoundaries(): List<Int> {
        val boundaries = listOf(24, 114, 206, 298, 337, 348, 365)
        val offsets = listOf(2, 2, 20, 20, 20, 20, 20)
        return boundaries.zip(offsets) { x, y -> x + y }
    }

    fun bBoundaries(): List<Int> {
        val boundaries = listOf(24, 114, 206, 298, 337, 348)
        val offsets = listOf(0, 13, 13, 13, 14, 14)
        return boundaries.zip(offsets) { x, y -> x + y }
    }

    fun cBoundaries(): List<Int> {
        val boundaries = listOf(24, 114, 206, 298, 338, 349, 365)
        val offsets = listOf(0, 5, 24, 24, 30, 36, 36)
        return boundaries.zip(offsets) { x, y -> x + y }
    }

}