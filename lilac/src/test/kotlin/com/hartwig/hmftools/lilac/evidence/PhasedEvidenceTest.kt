package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.seq.HlaSequence
import org.junit.Assert.assertEquals
import org.junit.Assert.assertTrue
import org.junit.Test

class PhasedEvidenceTest {

    @Test
    fun testInconsistentEvidence() {
        val victim = PhasedEvidence(listOf(0,1,3).toIntArray(), mapOf(Pair("CAT", 4), Pair("ATC", 5)))
        val catCandidate = HlaSequence("A*01:01", "CART")
        val atcCandidate = HlaSequence("A*01:02", "ATRC")
        val wildAtcCandidate = HlaSequence("A*01:03", "*TRC")
        val wildCandidate = HlaSequence("A*01:04", "****")

        val noMissing = victim.inconsistentEvidence(listOf(catCandidate, atcCandidate))
        assertEquals(0, noMissing.totalEvidence())

        val catMissing = victim.inconsistentEvidence(listOf(atcCandidate))
        assertEquals(1, catMissing.evidence.size)
        assertTrue(catMissing.evidence.containsKey("CAT"))

        val wildAtcMatch = victim.inconsistentEvidence(listOf(wildAtcCandidate))
        assertEquals(1, wildAtcMatch.evidence.size)
        assertTrue(wildAtcMatch.evidence.containsKey("CAT"))

        val wildMatch = victim.inconsistentEvidence(listOf(wildCandidate))
        assertEquals(0, wildMatch.totalEvidence())
    }

}