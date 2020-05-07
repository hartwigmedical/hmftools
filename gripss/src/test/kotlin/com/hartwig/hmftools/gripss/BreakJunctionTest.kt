package com.hartwig.hmftools.gripss

import org.junit.Assert.assertEquals
import org.junit.Test

class BreakJunctionTest {

    @Test
    fun testBreakPoint() {
        var victim = BreakJunction.create("AC", "ACTACCCCAACCTCCCCCAT[1:811432[")
        assertEquals(BreakPoint("TACCCCAACCTCCCCCAT",1, "1", 811432, -1), victim)

        victim = BreakJunction.create("A", "ACTACCCCAACCTCCCCCAT]1:811432]")
        assertEquals(BreakPoint("CTACCCCAACCTCCCCCAT", 1, "1", 811432, 1), victim)

        victim = BreakJunction.create("T", "[1:811432[ACTACCCCAACCTCCCCCAT")
        assertEquals(BreakPoint("ACTACCCCAACCTCCCCCA", -1, "1", 811432, -1), victim)

        victim = BreakJunction.create("AT", "]1:811432]ACTACCCCAACCTCCCCCAT")
        assertEquals(BreakPoint("ACTACCCCAACCTCCCCC", -1, "1", 811432, 1), victim)
    }


    @Test
    fun testBreakEnd() {
        var victim = BreakJunction.create("A", "ACTACCCCAACCTCCCCCAT.")
        assertEquals(BreakEnd("CTACCCCAACCTCCCCCAT", 1), victim)

        victim = BreakJunction.create("T", ".ACTACCCCAACCTCCCCCAT")
        assertEquals(BreakEnd("ACTACCCCAACCTCCCCCA", -1), victim)
    }
}