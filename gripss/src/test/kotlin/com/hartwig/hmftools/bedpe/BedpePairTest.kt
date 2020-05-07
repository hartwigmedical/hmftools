package com.hartwig.hmftools.bedpe

import org.junit.Assert.assertEquals
import org.junit.Test

class BedpePairTest {

    @Test
    fun testDecode() {
        val entry = "11\t118297204\t118397539\t10\t27035521\t27160016\tKMT2A-ABI1\t+\t-\t."
        val start = BedpeEntry("11", 118297205, 118397539, 1)
        val end = BedpeEntry("10", 27035522, 27160016, -1)

        assertEquals(BedpePair(start, end), BedpePair.fromString(entry))
    }

}