package com.hartwig.hmftools.lilackt.variant

import com.hartwig.hmftools.common.variant.CodingEffect
import com.hartwig.hmftools.lilackt.hla.HlaAllele
import com.hartwig.hmftools.lilackt.variant.SomaticCodingCount.Companion.addVariant
import org.junit.Assert.assertEquals
import org.junit.Test

class SomaticCodingCountTest {

    @Test
    fun testRepeatAllele() {
        val winningAlleles = listOf(
                HlaAllele("A*01:01"), HlaAllele("A*01:01"),
                HlaAllele("B*01:01"), HlaAllele("B*01:02"),
                HlaAllele("C*01:01"), HlaAllele("C*01:03"))

        var count = SomaticCodingCount.create(winningAlleles)
        assertEquals(6, count.size)

        count = count.addVariant(true, CodingEffect.MISSENSE, setOf(HlaAllele("A*01:01"), HlaAllele("B*01:01")))
        count = count.addVariant(false, CodingEffect.MISSENSE, setOf(HlaAllele("B*01:01")))

        assertEquals(0.5, count[0].total, 0.001)
        assertEquals(0.5, count[0].inframeIndel, 0.001)

        assertEquals(0.0, count[1].total, 0.001)
        assertEquals(1.5, count[2].total, 0.001)
        assertEquals(0.5, count[2].inframeIndel, 0.001)
        assertEquals(1.0, count[2].missense, 0.001)
    }

}