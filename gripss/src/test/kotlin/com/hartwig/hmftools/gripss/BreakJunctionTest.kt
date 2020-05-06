package com.hartwig.hmftools.gripss

import org.junit.Test

class BreakJunctionTest {


    @Test
    fun testAltPartner() {
//        val alt = PartnerAlt.create("[1:10178[GGTAGGG")
//        println(alt)

        val alt2 = BreakJunction.create("ACTACCCCAACCTCCCCCAT]1:811432]")
        println(alt2)

        val alt3 = BreakJunction.create("AAA.")
        println(alt3)
    }
}