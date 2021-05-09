package com.hartwig.hmftools.lilac.evidence;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class CombineEvidenceTest
{
    @Test
    public void testCheckSingles()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("T", 13);
        evidence.put("S", 71);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(344), evidence);

        evidence = Maps.newHashMap();
        evidence.put("S", 20);
        evidence.put("C", 20);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(348), evidence);

        evidence = Maps.newHashMap();
        evidence.put("SS", 9);
        evidence.put("SC", 4);
        PhasedEvidence combined = new PhasedEvidence(Lists.newArrayList(344, 348), evidence);

        assertFalse(CombineEvidence.canCombine(left, combined, right));
    }

    @Test
    public void testRealLife()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("AM", 36);
        evidence.put("LM", 16);
        evidence.put("RT", 11);
        evidence.put("RM", 35);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(1, 3), evidence);

        evidence = Maps.newHashMap();
        evidence.put("AMT", 31);
        evidence.put("AMA", 1);
        evidence.put("LMT", 15);
        evidence.put("RTT", 9);
        evidence.put("RTA", 1);
        evidence.put("RMA", 16);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(1,3,7), evidence);

        assertTrue(CombineEvidence.canCombine(left, right));
    }

    @Test
    public void canCombineExtraBaseToRight()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("AR", 2);
        evidence.put("TT", 3);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(3, 7), evidence);

        evidence = Maps.newHashMap();
        evidence.put("ARS", 2);
        evidence.put("TTH", 3);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(3, 7, 10), evidence);

        assertTrue(CombineEvidence.canCombine(left, right));
    }

    @Test
    public void cannotCombineExtraBaseToRight()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("AR", 2);
        evidence.put("TT", 3);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(3, 7), evidence);

        evidence = Maps.newHashMap();
        evidence.put("ARS", 2);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(3, 7, 10), evidence);
        assertFalse(CombineEvidence.canCombine(left, right));
    }

    @Test
    public void canCombineExtraBaseToLeft()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("HAR", 2);
        evidence.put("ATT", 3);

        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(2, 3, 7), evidence);
        evidence = Maps.newHashMap();
        evidence.put("AR", 2);
        evidence.put("TT", 3);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(3, 7), evidence);
        assertTrue(CombineEvidence.canCombine(left, right));
    }

    @Test
    public void canCombineOverlappingAndMatchingButNotWhenReversed()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("AR", 2);
        evidence.put("TT", 3);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(3, 7), evidence);

        evidence = Maps.newHashMap();
        evidence.put("RS", 2);
        evidence.put("TS", 3);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(7, 10), evidence);
        assertTrue(CombineEvidence.canCombine(left, right));
        assertFalse(CombineEvidence.canCombine(right, left));
    }

    @Test
    public void cannotCombineWhenSequencesAreNotAccountedFor()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("AR", 2);
        evidence.put("TT", 3);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(3, 7), evidence);

        evidence = Maps.newHashMap();
        evidence.put("SS", 2);
        evidence.put("TS", 3);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(7, 1), evidence);
        assertFalse(CombineEvidence.canCombine(left, right));
    }

    @Test
    public void cannotCombineNonOverlapping()
    {
        Map<String,Integer> evidence = Maps.newHashMap();
        evidence.put("AR", 2);
        evidence.put("TR", 3);
        PhasedEvidence left = new PhasedEvidence(Lists.newArrayList(3, 7), evidence);

        evidence = Maps.newHashMap();
        evidence.put("TS", 2);
        evidence.put("QR", 3);
        PhasedEvidence right = new PhasedEvidence(Lists.newArrayList(9, 10), evidence);
        assertFalse(CombineEvidence.canCombine(left, right));
    }

}
