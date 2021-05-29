package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.lilac.LilacConstants.WARN_UNMATCHED_HAPLOTYPE_SUPPORT;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.SequenceCount;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.qc.AminoAcidQC;
import com.hartwig.hmftools.lilac.qc.Haplotype;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.commons.math3.util.Pair;
import org.junit.Assert;
import org.junit.Test;

public class QualityControlTest
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
                victim, 0, Lists.newArrayList(catCandidate, atcCandidate), aminoAcidCount, Lists.newArrayList());
        Assert.assertEquals(0, noMissing.size());

        List<Haplotype> catMissing = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(atcCandidate), aminoAcidCount, Lists.newArrayList());
        Assert.assertEquals(1, catMissing.size());
        Assert.assertTrue(catMissing.get((0)).Haplotype.equals("CART"));

        List<Haplotype> catNotMissingBecauseOfMinEvidence = HaplotypeQC.unmatchedHaplotype(
                victim, 5, Lists.newArrayList(atcCandidate), aminoAcidCount, Lists.newArrayList());
        Assert.assertEquals(0, catNotMissingBecauseOfMinEvidence.size());

        List<Haplotype> wildAtcMatch = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(wildAtcCandidate), aminoAcidCount, Lists.newArrayList());
        Assert.assertEquals(1, wildAtcMatch.size());
        Assert.assertTrue(wildAtcMatch.get((0)).Haplotype.equals("CART"));

        List<Haplotype> wildMatch = HaplotypeQC.unmatchedHaplotype(
                victim, 0, Lists.newArrayList(wildCandidate), aminoAcidCount, Lists.newArrayList());
        Assert.assertEquals(0, wildMatch.size());
    }

    @Test
    public void testUnmatchedAminoAcid()
    {
        // look for AAs not in unmatched haplotypes and not in wildcard sequences

        HlaSequenceLoci winningSeq1 = createSequenceLoci("A*01:01", "ABCDEFGHIJKLMNO");
        HlaSequenceLoci winningSeq2 = createSequenceLoci("A*01:02", "ABCDEFGHIJKLM**");
        HlaSequenceLoci winningSeq3 = createSequenceLoci("A*01:02", "A**DEFGHIJKLMNN");

        Map<String,Integer>[] sequenceCountsMap = new Map[winningSeq1.length()];

        for(int i = 0; i < winningSeq1.length(); ++i)
        {
            Map<String,Integer> locusMap = Maps.newHashMap();
            locusMap.put(winningSeq1.sequence(i), 10);
            sequenceCountsMap[i] = locusMap;
        }

        // add the unmatched AAs

        // in wildcard sequences, so not reported
        int warnLevel = (int)(WARN_UNMATCHED_HAPLOTYPE_SUPPORT * 1000);
        sequenceCountsMap[1].put("C", warnLevel); // expect B
        sequenceCountsMap[14].put("C", warnLevel); // expect anything

        // in an unmatched haplotype, so not reported
        sequenceCountsMap[7].put("C", warnLevel); // expect G

        // too low so not reported
        sequenceCountsMap[3].put("C", warnLevel - 1); // expect D

        // reportable
        sequenceCountsMap[0].put("B", warnLevel); // expect A
        sequenceCountsMap[12].put("B", warnLevel); // expect M

        SequenceCount aminoAcidCount = new SequenceCount(1, sequenceCountsMap);

        List<Haplotype> unmatchedHaplotypes = Lists.newArrayList(
                new Haplotype(5, 10, 10, "ABCDEF"));

        List<HlaSequenceLoci> winningSequences = Lists.newArrayList(winningSeq1, winningSeq2, winningSeq3);

        AminoAcidQC qc = AminoAcidQC.create(winningSequences, Lists.newArrayList(), aminoAcidCount, unmatchedHaplotypes, 1000);
        assertEquals(2, qc.UnusedAminoAcids);
    }

    private HlaSequenceLoci createSequenceLoci(final String allele, final String sequenceStr)
    {
        HlaSequence sequence = new HlaSequence(HlaAllele.fromString(allele), sequenceStr);
        return HlaSequenceLoci.create(sequence.Allele, sequence.getRawSequence(), sequence.getRawSequence());
    }
}
