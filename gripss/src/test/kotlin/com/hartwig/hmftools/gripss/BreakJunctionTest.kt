package com.hartwig.hmftools.gripss

import org.junit.Assert.assertEquals
import org.junit.Test

class BreakJunctionTest {

    @Test
    fun testBreakPoint() {
        var victim = BreakJunction.create("ACTACCCCAACCTCCCCCAT[1:811432[")
        assertEquals(BreakPoint(1, "1", 811432, -1), victim)

        victim = BreakJunction.create("ACTACCCCAACCTCCCCCAT]1:811432]")
        assertEquals(BreakPoint(1, "1", 811432, 1), victim)

        victim = BreakJunction.create("[1:811432[ACTACCCCAACCTCCCCCAT")
        assertEquals(BreakPoint(-1, "1", 811432, -1), victim)

        victim = BreakJunction.create("]1:811432]ACTACCCCAACCTCCCCCAT")
        assertEquals(BreakPoint(-1, "1", 811432, 1), victim)
    }


    @Test
    fun testBreakEnd() {
        var victim = BreakJunction.create("ACTACCCCAACCTCCCCCAT.")
        assertEquals(BreakEnd(1), victim)

        victim = BreakJunction.create(".ACTACCCCAACCTCCCCCAT")
        assertEquals(BreakEnd(-1), victim)
    }
}