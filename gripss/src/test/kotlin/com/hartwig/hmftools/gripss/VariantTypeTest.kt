package com.hartwig.hmftools.gripss

import org.junit.Assert.assertEquals
import org.junit.Test

class VariantTypeTest {

    @Test
    fun testBnd() {
        val startAlt = "AAT[2:110["
        val start = VariantType.create("1", 100, "A", "AAT[2:110[")
        assertEquals(Translocation("AT","2", 110, 1, -1), start)

        val end = VariantType.create("2", 110, "T", "]1:100]AAT")
        assertEquals(Translocation("AA","1", 100, -1, 1), end)
    }

    @Test
    fun testDup() {
        val start = VariantType.create("1", 110, "A", "AAT[1:100[")
        assertEquals(Duplication("AT","1", 100, 1, -1, 10), start)

        val end = VariantType.create("1", 100, "T", "]1:110]AAT")
        assertEquals(Duplication("AA","1", 110, -1, 1, 10), end)
    }

    @Test
    fun testDel() {
        val start = VariantType.create("1", 100, "A", "AAT[1:110[")
        assertEquals(Deletion("AT","1", 110, 1, -1, 10), start)

        val end = VariantType.create("1", 110, "T", "]1:100]AAT")
        assertEquals(Deletion("AA","1", 100, -1, 1, 10), end)
    }

    @Test
    fun testIns() {
        val start = VariantType.create("1", 100, "A", "AAT[1:101[")
        assertEquals(Insertion("AT","1", 101, 1, -1, 1), start)

        val end = VariantType.create("1", 101, "T", "]1:100]AAT")
        assertEquals(Insertion("AA","1", 100, -1, 1, 1), end)
    }

}