package com.hartwig.hmftools.gripss

import org.junit.Assert.assertEquals
import org.junit.Test

class VariantTypeTest {

    @Test
    fun testBnd() {
        val start = VariantType.create("1", 100, "A", "AAT[2:110[")
        assertEquals(Translocation("A", "AT","2", 110, 1, -1), start)

        val end = VariantType.create("2", 110, "T", "]1:100]AAT")
        assertEquals(Translocation("T", "AA","1", 100, -1, 1), end)
    }

    @Test
    fun testDup() {
        val start = VariantType.create("1", 110, "A", "AAT[1:100[")
        assertEquals(Duplication("A", "AT","1", 100, 1, -1, 10), start)

        val end = VariantType.create("1", 100, "T", "]1:110]AAT")
        assertEquals(Duplication("T", "AA","1", 110, -1, 1, 10), end)
    }

    @Test
    fun testDel() {
        val start = VariantType.create("1", 100, "A", "AAT[1:110[")
        assertEquals(Deletion("A", "AT","1", 110, 1, -1, 10), start)

        val end = VariantType.create("1", 110, "T", "]1:100]AAT")
        assertEquals(Deletion("T", "AA","1", 100, -1, 1, 10), end)
    }

    @Test
    fun testIns() {
        val start = VariantType.create("1", 100, "A", "AAT[1:101[")
        assertEquals(Insertion("A", "AT","1", 101, 1, -1, 1), start)

        val end = VariantType.create("1", 101, "T", "]1:100]AAT")
        assertEquals(Insertion("T", "AA","1", 100, -1, 1, 1), end)
    }

    @Test
    fun testToString() {
        assertToString(".GATAC")
        assertToString("GATAC.")
        assertToString("CAT[1:3123[")
        assertToString("CAT]1:3123]")
        assertToString("[1:3123[CAT")
        assertToString("]1:3123]CAT")
    }

    private fun assertToString(expected: String) {
        val variant = VariantType.create("1", 100, "A", expected)
        assertEquals(expected, variant.toString())
    }
}