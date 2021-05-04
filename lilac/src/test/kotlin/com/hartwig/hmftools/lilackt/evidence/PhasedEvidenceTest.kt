package com.hartwig.hmftools.lilackt.evidence

import com.hartwig.hmftools.lilackt.seq.HlaSequence
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci
import org.junit.Assert.assertEquals
import org.junit.Assert.assertTrue
import org.junit.Test

class PhasedEvidenceTest {

    @Test
    fun testInconsistentEvidence() {
        val victim = PhasedEvidence(listOf(0,1,3).toIntArray(), mapOf(Pair("CAT", 4), Pair("ATC", 5)))
        val catCandidate = create(HlaSequence("A*01:01", "CART"))
        val atcCandidate = create(HlaSequence("A*01:02", "ATRC"))
        val wildAtcCandidate = create(HlaSequence("A*01:03", "*TRC"))
        val wildCandidate = create(HlaSequence("A*01:04", "****"))

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

    fun create(sequences: HlaSequence): HlaSequenceLoci {
        return HlaSequenceLoci.create(sequences.allele, sequences.rawSequence, sequences.rawSequence)
    }

}