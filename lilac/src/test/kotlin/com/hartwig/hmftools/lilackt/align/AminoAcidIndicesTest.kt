package com.hartwig.hmftools.lilackt.align

import com.hartwig.hmftools.lilackt.read.AminoAcidIndices
import junit.framework.Assert.assertEquals
import org.junit.Test

class AminoAcidIndicesTest {

    @Test
    fun testIndices() {
        assertRange(0, -1, AminoAcidIndices.indices(0, 0))
        assertRange(0, -1, AminoAcidIndices.indices(0, 1))
        assertRange(0, 0, AminoAcidIndices.indices(0, 2))
        assertRange(1, 0, AminoAcidIndices.indices(1, 2))
        assertRange(1, 0, AminoAcidIndices.indices(2, 2))

        assertRange(1, 1, AminoAcidIndices.indices(3, 5))
        assertRange(1, 1, AminoAcidIndices.indices(3, 6))
        assertRange(1, 1, AminoAcidIndices.indices(3, 7))
        assertRange(1, 2, AminoAcidIndices.indices(3, 8))
        assertRange(1, 2, AminoAcidIndices.indices(3, 9))
        assertRange(1, 2, AminoAcidIndices.indices(3, 10))

    }

    private fun assertRange(expectedStart: Int, expectedEnd: Int, victim: IntRange) {
        assertEquals(expectedStart, victim.first)
        assertEquals(expectedEnd, victim.last)
    }

}