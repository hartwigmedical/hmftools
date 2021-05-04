package com.hartwig.hmftools.lilackt.hla

import junit.framework.Assert.assertEquals
import org.junit.Test

class HlaAlleleTest {

    @Test
    fun testDecode() {
        assertContig("A*01:01:01:01", "A*01:01:01:01")
        assertContig("A*01:01:01", "A*01:01:01")
        assertContig("A*26:22", "A*26:22")
        assertContig("A*26", "A*26")
    }

    @Test
    fun testReduce() {
        val eightDigit = HlaAllele("A*01:01:01:01")
        assertEquals(HlaAllele("A*01:01:01"), eightDigit.asSixDigit())
        assertEquals(HlaAllele("A*01:01"), eightDigit.asFourDigit())
        assertEquals(HlaAllele("A*01"), eightDigit.asAlleleGroup())

    }


    fun assertContig(expected: String, contig: String) {
        assertEquals(expected, HlaAllele(contig).toString())
    }

}