package com.hartwig.hmftools.lilac.hla

import junit.framework.Assert.assertEquals
import org.junit.Test

class HlaAlleleTest {

    @Test
    fun testDecode() {
        assertContig("A*01:01:01:01", "A*01:01:01:01")
        assertContig("A*26:22", "A*26:22")
    }


    fun assertContig(expected: String, contig: String) {
        assertEquals(expected, HlaAllele(contig).toString())
    }

}