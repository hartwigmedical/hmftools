package com.hartwig.hmftools.teal.util

import kotlin.test.*

internal class TelomereMatcherTest
{
    fun setUp()
    {
        //Configurator.setRootLevel(Level.TRACE)
    }

    @Test
    fun testCalcTelomereMatch()
    {
        setUp()

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
    fun testCalcTelomereMatchShort()
    {
        setUp()

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

    @Test
    fun testCalcTelomereMatch2()
    {
        setUp()

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
        setUp()

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
        setUp()

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
