package com.hartwig.hmftools.paddle.dnds

import org.junit.Assert.assertEquals
import org.junit.Test

class DndsCvTest {

    @Test
    fun testExpectedDrivers() {
        val victim = DndsCv(115, 3.4091126)
        assertEquals(72.78686, victim.expectedDrivers(103), 0.00001)
    }
}
