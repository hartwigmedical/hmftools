package com.hartwig.hmftools.lilackt.evidence

import org.junit.Assert.assertFalse
import org.junit.Assert.assertTrue
import org.junit.Test

class CombineEvidenceTest {

    @Test
    fun testCheckSingles() {
        val left = PhasedEvidence(listOf(344).toIntArray(), mapOf(Pair("T", 13), Pair("S", 71)))
        val right = PhasedEvidence(listOf(348).toIntArray(), mapOf(Pair("S", 20), Pair("C", 20)))
        val combined = PhasedEvidence(listOf(344, 348).toIntArray(), mapOf(Pair("SS", 9), Pair("SC", 4)))
        assertFalse(CombineEvidence.canCombine(left, combined, right))
    }

    @Test
    fun testRealLife() {
        val left = PhasedEvidence(listOf(1, 3).toIntArray(), mapOf(Pair("AM", 36), Pair("LM", 16),  Pair("RT", 11),  Pair("RM", 35)))
        val right = PhasedEvidence(listOf(1, 3, 7).toIntArray(), mapOf(Pair("AMT", 31), Pair("AMA", 1), Pair("LMT", 15),  Pair("RTT", 9),  Pair("RTA", 1),  Pair("RMA", 16)))
        assertTrue(CombineEvidence.canCombine(left, right))

    }

    @Test
    fun canCombineExtraBaseToRight() {
        val left = PhasedEvidence(listOf(3, 7).toIntArray(), mapOf(Pair("AR", 2), Pair("TT", 3)))
        val right = PhasedEvidence(listOf(3, 7, 10).toIntArray(), mapOf(Pair("ARS", 2), Pair("TTH", 3)))
        assertTrue(CombineEvidence.canCombine(left, right))
    }

    @Test
    fun cannotCombineExtraBaseToRight() {
        val left = PhasedEvidence(listOf(3, 7).toIntArray(), mapOf(Pair("AR", 2), Pair("TT", 3)))
        val right = PhasedEvidence(listOf(3, 7, 10).toIntArray(), mapOf(Pair("ARS", 2)))
        assertFalse(CombineEvidence.canCombine(left, right))
    }

    @Test
    fun canCombineExtraBaseToLeft() {
        val left = PhasedEvidence(listOf(2, 3, 7).toIntArray(), mapOf(Pair("HAR", 2), Pair("ATT", 3)))
        val right = PhasedEvidence(listOf(3, 7).toIntArray(), mapOf(Pair("AR", 2), Pair("TT", 3)))
        assertTrue(CombineEvidence.canCombine(left, right))
    }

    @Test
    fun canCombineOverlappingAndMatchingButNotWhenReversed() {
        val left = PhasedEvidence(listOf(3, 7).toIntArray(), mapOf(Pair("AR", 2), Pair("TT", 3)))
        val right = PhasedEvidence(listOf(7, 10).toIntArray(), mapOf(Pair("RS", 2), Pair("TS", 3)))
        assertTrue(CombineEvidence.canCombine(left, right))
        assertFalse(CombineEvidence.canCombine(right, left))
    }

    @Test
    fun cannotCombineWhenSequencesAreNotAccountedFor() {
        val left = PhasedEvidence(listOf(3, 7).toIntArray(), mapOf(Pair("AR", 2), Pair("TT", 3)))
        val right = PhasedEvidence(listOf(7, 10).toIntArray(), mapOf(Pair("SS", 2), Pair("TS", 3)))
        assertFalse(CombineEvidence.canCombine(left, right))
    }

    @Test
    fun cannotCombineNonOverlapping() {
        val left = PhasedEvidence(listOf(3, 7).toIntArray(), mapOf(Pair("AR", 2), Pair("TR", 3)))
        val right = PhasedEvidence(listOf(9, 10).toIntArray(), mapOf(Pair("TS", 2), Pair("QR", 3)))
        assertFalse(CombineEvidence.canCombine(left, right))
    }

}