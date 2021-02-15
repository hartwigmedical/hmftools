package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.evidence.PhasedEvidence
import com.hartwig.hmftools.lilac.qc.HaplotypeQC.Companion.unmatchedHaplotype
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import org.junit.Assert
import org.junit.Test

class HaplotypeQCTest {

    @Test
    fun testUnmatchedHaplotype() {
        val victim = PhasedEvidence(listOf(0,1,3).toIntArray(), mapOf(Pair("CAT", 4), Pair("ATC", 5)))
        val catCandidate = create(HlaSequence("A*01:01", "CART"))
        val atcCandidate = create(HlaSequence("A*01:02", "ATRC"))
        val wildAtcCandidate = create(HlaSequence("A*01:03", "*TRC"))
        val wildCandidate = create(HlaSequence("A*01:04", "****"))

        val noMissing = victim.unmatchedHaplotype(0,listOf(catCandidate, atcCandidate))
        Assert.assertEquals(0, noMissing.haplotypes.size)

        val catMissing = victim.unmatchedHaplotype(4,listOf(atcCandidate))
        Assert.assertEquals(1, catMissing.haplotypes.size)
        Assert.assertTrue(catMissing.haplotypes.containsKey("CAT"))

        val catNotMissingBecauseOfMinEvidence = victim.unmatchedHaplotype(5,listOf(atcCandidate))
        Assert.assertEquals(0, catNotMissingBecauseOfMinEvidence.haplotypes.size)

        val wildAtcMatch = victim.unmatchedHaplotype(0,listOf(wildAtcCandidate))
        Assert.assertEquals(1, wildAtcMatch.haplotypes.size)
        Assert.assertTrue(wildAtcMatch.haplotypes.containsKey("CAT"))

        val wildMatch = victim.unmatchedHaplotype(0,listOf(wildCandidate))
        Assert.assertEquals(0, wildMatch.totalEvidence())
    }

    private fun create(sequences: HlaSequence): HlaSequenceLoci {
        return HlaSequenceLoci.create(sequences.allele, sequences.rawSequence, sequences.rawSequence)
    }

}