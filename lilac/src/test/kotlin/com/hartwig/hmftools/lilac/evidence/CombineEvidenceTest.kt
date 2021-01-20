package com.hartwig.hmftools.lilac.evidence

import org.junit.Assert.assertFalse
import org.junit.Assert.assertTrue
import org.junit.Test

class CombineEvidenceTest {

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