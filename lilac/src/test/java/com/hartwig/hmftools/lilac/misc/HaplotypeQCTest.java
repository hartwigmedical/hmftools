package com.hartwig.hmftools.lilac.misc;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class HaplotypeQCTest
{
    @Test
    public void testUnmatchedHaplotype()
    {
        val index0Map = mapOf(Pair("C", 4), Pair("A", 4))
        val index1Map = mapOf(Pair("A", 4), Pair("T", 4))
        val index2Map = mapOf(Pair("R", 8))
        val index3Map = mapOf(Pair("T", 4), Pair("C", 4))
        val array = Lists.newArrayList(index0Map, index1Map, index2Map, index3Map)

        val aminoAcidCount = SequenceCount(1, array)
        val victim = PhasedEvidence(listOf(0, 1, 3).toIntArray(), mapOf(Pair("CAT", 4), Pair("ATC", 5)))
        val catCandidate = create(HlaSequence("A*01:01", "CART"))
        val atcCandidate = create(HlaSequence("A*01:02", "ATRC"))
        val wildAtcCandidate = create(HlaSequence("A*01:03", "*TRC"))
        val wildCandidate = create(HlaSequence("A*01:04", "****"))

        val noMissing = victim.unmatchedHaplotype(0, listOf(catCandidate, atcCandidate), aminoAcidCount)
        Assert.assertEquals(0, noMissing.size)

        val catMissing = victim.unmatchedHaplotype(4, listOf(atcCandidate), aminoAcidCount)
        Assert.assertEquals(1, catMissing.size)
        Assert.assertTrue(catMissing[0].haplotype == "CART")

        val catNotMissingBecauseOfMinEvidence = victim.unmatchedHaplotype(5, listOf(atcCandidate), aminoAcidCount)
        Assert.assertEquals(0, catNotMissingBecauseOfMinEvidence.size)

        val wildAtcMatch = victim.unmatchedHaplotype(0, listOf(wildAtcCandidate), aminoAcidCount)
        Assert.assertEquals(1, wildAtcMatch.size)
        Assert.assertTrue(wildAtcMatch[0].haplotype == "CART")

        val wildMatch = victim.unmatchedHaplotype(0, listOf(wildCandidate), aminoAcidCount)
        Assert.assertEquals(0, wildMatch.size)
    }

    private fun create(sequences:HlaSequence):HlaSequenceLoci

    {
        return HlaSequenceLoci.create(sequences.allele, sequences.rawSequence, sequences.rawSequence)
    }

}
