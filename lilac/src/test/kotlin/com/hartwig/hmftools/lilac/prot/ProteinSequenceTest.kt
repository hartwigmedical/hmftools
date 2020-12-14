package com.hartwig.hmftools.lilac.prot

import junit.framework.Assert.assertEquals
import junit.framework.Assert.assertTrue
import org.junit.Test

class ProteinSequenceTest {

    @Test
    fun testExonicBoundaryInMiddle() {
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(3, listOf(3))
        assertEquals(2, exonicKmers.size)
        assertTrue("ABC" in exonicKmers)
        assertTrue("EFG" in exonicKmers)
    }

    @Test
    fun testExonicBoundaryAtEnd() {
        val proteins = "ABCDEFG"
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(5, listOf(proteins.length - 1))
        assertEquals(2, exonicKmers.size)
        assertTrue("ABCDE" in exonicKmers)
        assertTrue("BCDEF" in exonicKmers)
    }

    @Test
    fun testExonicBoundaryPastEnd() {
        val proteins = "ABCDEFG"
        val victim = ProteinSequence("A*Test", proteins)
        val exonicKmers = victim.uniqueExonicKmers(6, listOf(proteins.length))
        assertEquals(2, exonicKmers.size)
        assertTrue("ABCDEF" in exonicKmers)
        assertTrue("BCDEFG" in exonicKmers)
    }

    @Test
    fun testExonicBoundaryAtStart() {
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(6, listOf(0))
        assertEquals(1, exonicKmers.size)
        assertTrue("BCDEFG" in exonicKmers)
    }

    @Test
    fun testRollingKmers4() {
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(4, listOf())
        assertEquals(4, exonicKmers.size)
        assertTrue("ABCD" in exonicKmers)
        assertTrue("BCDE" in exonicKmers)
        assertTrue("CDEF" in exonicKmers)
        assertTrue("DEFG" in exonicKmers)
    }

    @Test
    fun testRollingKmers6() {
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(6, listOf())
        assertEquals(2, exonicKmers.size)
        assertTrue("ABCDEF" in exonicKmers)
        assertTrue("BCDEFG" in exonicKmers)
    }


    @Test
    fun testRollingKmersMinLength() {
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(7, 20, listOf())
        assertEquals(1, exonicKmers.size)
        assertTrue("ABCDEFG" in exonicKmers)
    }

    @Test
    fun testRollingKmersInsufficientMinLength() {
        val victim = ProteinSequence("A*Test", "ABCDEFG")
        val exonicKmers = victim.uniqueExonicKmers(8, 20, listOf())
        assertEquals(0, exonicKmers.size)
    }

}