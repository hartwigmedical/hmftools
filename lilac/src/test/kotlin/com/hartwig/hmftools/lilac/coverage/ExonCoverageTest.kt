package com.hartwig.hmftools.lilac.coverage

import org.junit.Assert.assertEquals
import org.junit.Assert.assertTrue
import org.junit.Test

class ExonCoverageTest {

    @Test
    fun testCoverage() {
        val victim = ExonCoverage("POPSICLE")
        victim.addCoverageFromEnd(3)
        assertEquals(intArrayOf(0, 0, 0, 0, 0, 1, 1, 1).joinToString { it.toString() }, victim.coverage.joinToString { it.toString() })

        victim.addCoverageFromEnd(4)
        assertEquals(intArrayOf(0, 0, 0, 0, 1, 2, 2, 2).joinToString { it.toString() }, victim.coverage.joinToString { it.toString() })

        victim.addCoverageFromStart(1)
        assertEquals(intArrayOf(1, 0, 0, 0, 1, 2, 2, 2).joinToString { it.toString() }, victim.coverage.joinToString { it.toString() })
        assertEquals(0, victim.minCoverage())
        assertEquals(2, victim.maxCoverage())

        victim.addCompleteCoverage()
        assertEquals(intArrayOf(2, 1, 1, 1, 2, 3, 3, 3).joinToString { it.toString() }, victim.coverage.joinToString { it.toString() })
        assertEquals(1, victim.minCoverage())
        assertEquals(3, victim.maxCoverage())
        assertTrue(victim.isCovered())
    }

}