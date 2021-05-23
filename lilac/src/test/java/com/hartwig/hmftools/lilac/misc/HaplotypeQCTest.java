package com.hartwig.hmftools.lilac.misc;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.qc.Haplotype;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.Test;

public class HaplotypeQCTest
{
    @Test
    public void testUnmatchedHaplotype()
    {
        Map<String,Integer> index0Map = Maps.newHashMap();
        index0Map.put("C", 4);
        index0Map.put("A", 4);

        Map<String,Integer> index1Map = Maps.newHashMap();
        index1Map.put("A", 4);
        index1Map.put("T", 4);

        Map<String,Integer> index2Map = Maps.newHashMap();
        index2Map.put("R", 8);

        Map<String,Integer> index3Map = Maps.newHashMap();
        index3Map.put("T", 4);
        index3Map.put("C", 4);

        Map<String,Integer>[] array = new Map[] { index0Map, index1Map, index2Map, index3Map };

        SequenceCount aminoAcidCount = new SequenceCount(1, array);

        List<Integer> aminoAcidIndexList = Lists.newArrayList(0, 1, 3);
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("CAT", 4);
        evidence.put("ATC", 5);
        PhasedEvidence victim = new PhasedEvidence(aminoAcidIndexList, evidence);

        HlaSequenceLoci catCandidate = createSequenceLoci("A*01:01", "CART");
        HlaSequenceLoci atcCandidate = createSequenceLoci("A*01:02", "ATRC");
        HlaSequenceLoci wildAtcCandidate = createSequenceLoci("A*01:03", "*TRC");
        HlaSequenceLoci wildCandidate = createSequenceLoci("A*01:04", "****");

        List<Haplotype> noMissing = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(catCandidate, atcCandidate), aminoAcidCount);
        Assert.assertEquals(0, noMissing.size());

        List<Haplotype> catMissing = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(atcCandidate), aminoAcidCount);
        Assert.assertEquals(1, catMissing.size());
        Assert.assertTrue(catMissing.get((0)).Haplotype.equals("CART"));

        List<Haplotype> catNotMissingBecauseOfMinEvidence = HaplotypeQC.unmatchedHaplotype(
                victim, 5, Lists.newArrayList(atcCandidate), aminoAcidCount);
        Assert.assertEquals(0, catNotMissingBecauseOfMinEvidence.size());

        List<Haplotype> wildAtcMatch = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(wildAtcCandidate), aminoAcidCount);
        Assert.assertEquals(1, wildAtcMatch.size());
        Assert.assertTrue(wildAtcMatch.get((0)).Haplotype.equals("CART"));

        List<Haplotype> wildMatch = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(wildCandidate), aminoAcidCount);
        Assert.assertEquals(0, wildMatch.size());
    }

    private HlaSequenceLoci createSequenceLoci(final String allele, final String sequenceStr)
    {
        HlaSequence sequence = new HlaSequence(HlaAllele.fromString(allele), sequenceStr);
        return HlaSequenceLoci.create(sequence.Allele, sequence.getRawSequence(), sequence.getRawSequence());
    }
}
