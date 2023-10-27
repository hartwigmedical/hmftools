package com.hartwig.hmftools.teal.util

import org.junit.Before
import kotlin.test.*

internal class TelomereMatcherTest
{
    @Before
    fun setUp()
    {
        // org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE)
    }

    @Test
    fun testCalcTelomereMatch()
    {
        // try a simple one
        val seq = "ACTACCATTAGGGTTAGGGACTA"

        val telomereMatch = TelomereMatcher.findGTelomereSegment(seq, 0.9)
        assertNotNull(telomereMatch)
        assertEquals(telomereMatch.numMatchedBases, 12)
        assertEquals(telomereMatch.matchStart, 7)
        assertEquals(telomereMatch.matchEnd, 19)
        assertEquals(telomereMatch.matchedSequence, "TTAGGGTTAGGG")

        assertTrue(TelomereMatcher.matchesGTelomere(seq, 0.9, 12))
    }

    @Test
    fun testCalcTelomereNoMatch()
    {
        // no match
        val seq = "TGCAGCTTAACTGAGAGCCGCTCCTCTTCTCT"

        val telomereMatch = TelomereMatcher.findGTelomereSegment(seq, 0.9)
        assertNull(telomereMatch)
        assertFalse(TelomereMatcher.matchesGTelomere(seq, 0.9, 12))
    }

    /*
    @Test
    fun testCalcTelomereMatchShort()
    {
        // try a simple one
        val seq = "ACTATTAGGGCACTA"

        val telomereMatch = TelomereMatcher.findGTelomereSegment(seq, 0.9)
        assertNotNull(telomereMatch)
        assertEquals(telomereMatch.numMatchedBases, 6)
        assertEquals(telomereMatch.matchStart, 4)
        assertEquals(telomereMatch.matchEnd, 10)
        assertEquals(telomereMatch.matchedSequence, "TTAGGG")

        assertFalse(TelomereMatcher.matchesGTelomere(seq, 0.9, 12))
    }
     */

    @Test
    fun testCalcTelomereMatch2()
    {
        // try a simple one
        val seq = "ACTATTAGGGCTTAGGGACTA"

        val telomereMatch = TelomereMatcher.findGTelomereSegment(seq, 0.9)
        assertNotNull(telomereMatch)
        assertEquals(telomereMatch.numMatchedBases, 12)
        assertEquals(telomereMatch.matchStart, 4)
        assertEquals(telomereMatch.matchEnd, 17)
        assertEquals(telomereMatch.matchedSequence, "TTAGGGCTTAGGG")

        assertTrue(TelomereMatcher.matchesGTelomere(seq, 0.9, 12))
    }

    @Test
    fun testTelomereMatch3()
    {
        val seq = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGGTAGGGTGAGGGTTAGGGTGAGGGTTAGGGTGAGGGTGAGGGTTAGGGTGAGGGTTAGGGTGAGGGAGAGGGTTAGGG"

        val telomereMatch = TelomereMatcher.findGTelomereSegment(seq, 0.9)
        assertNotNull(telomereMatch)
        assertEquals(telomereMatch.numMatchedBases, 110)
        assertEquals(telomereMatch.matchStart, 0)
        assertEquals(telomereMatch.matchEnd, 120)
        assertEquals(telomereMatch.matchedSequence, seq)

        assertTrue(TelomereMatcher.matchesGTelomere(seq, 0.9, 12))
    }


    @Test
    fun testCalcTelomereMatch4()
    {
        val seq = "GTACTATTAGGGCTTAGGGACTA"

        val telomereMatch = TelomereMatcher.findGTelomereSegment(seq, 0.9)
        assertNotNull(telomereMatch)
        assertEquals(telomereMatch.numMatchedBases, 12)
        assertEquals(telomereMatch.matchStart, 6)
        assertEquals(telomereMatch.matchEnd, 19)
        assertEquals(telomereMatch.matchedSequence, "TTAGGGCTTAGGG")

        assertTrue(TelomereMatcher.matchesGTelomere(seq, 0.9, 12))
    }
}
