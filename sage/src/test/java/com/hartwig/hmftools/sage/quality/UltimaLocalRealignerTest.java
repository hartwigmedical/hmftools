package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.pairHomopolymers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.Homopolymer;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.HomopolymerPair;

import org.junit.Test;

public class UltimaLocalRealignerTest
{
    @Test
    public void testPairHomopolymersEmpty()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList();
        List<Homopolymer> readHomopolymers = Lists.newArrayList();

        List<HomopolymerPair> pairs = pairHomopolymers(refHomopolymers, readHomopolymers);

        assertTrue(pairs.isEmpty());
    }

    @Test
    public void testPairHomopolymersSingleMatch()
    {
        Homopolymer a1 = new Homopolymer('A', 1);
        Homopolymer a3 = new Homopolymer('A', 3);

        List<Homopolymer> refHomopolymers = Lists.newArrayList(a1);
        List<Homopolymer> readHomopolymers = Lists.newArrayList(a3);

        List<HomopolymerPair> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
        List<HomopolymerPair> expectedPairs = Lists.newArrayList(new HomopolymerPair(a1, a3));

        assertEquals(expectedPairs, actualPairs);
    }

    // TODO: Re-enable.
//    @Test
//    public void testPairHomopolymersSingleMatchmatch()
//    {
//        Homopolymer a1 = new Homopolymer('A', 1);
//        Homopolymer t1 = new Homopolymer('T', 1);
//
//        List<Homopolymer> refHomopolymers = Lists.newArrayList(a1);
//        List<Homopolymer> readHomopolymers = Lists.newArrayList(t1);
//
//        List<HomopolymerPair> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//
//        assertTrue(
//                actualPairs.equals(Lists.newArrayList(new HomopolymerPair(a1, null), new HomopolymerPair(null, t1)))
//                || actualPairs.equals(Lists.newArrayList(new HomopolymerPair(null, t1), new HomopolymerPair(a1, null))));
//
//        // TODO: Remove me
//        assertTrue(false);
//    }

    @Test
    public void testPairHomopolymersMultipleMatch()
    {
        Homopolymer a1 = new Homopolymer('A', 1);
        Homopolymer a3 = new Homopolymer('A', 3);
        Homopolymer t1 = new Homopolymer('T', 1);
        Homopolymer t2 = new Homopolymer('T', 2);

        List<Homopolymer> refHomopolymers = Lists.newArrayList(a1, t2);
        List<Homopolymer> readHomopolymers = Lists.newArrayList(a3, t1);

        List<HomopolymerPair> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
        List<HomopolymerPair> expectedPairs = Lists.newArrayList(new HomopolymerPair(a1, a3), new HomopolymerPair(t2, t1));

        assertEquals(expectedPairs, actualPairs);
    }

    @Test
    public void testPairHomopolymersContractionAndSnv()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer('C', 2),
                new Homopolymer('T', 1),
                new Homopolymer('C', 3),
                new Homopolymer('G', 2),
                new Homopolymer('A', 1));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer('C', 2),
                new Homopolymer('T', 1),
                new Homopolymer('C', 2),
                new Homopolymer('T', 1),
                new Homopolymer('G', 1),
                new Homopolymer('A', 1));

        List<HomopolymerPair> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
        List<HomopolymerPair> expectedPairs = Lists.newArrayList(
                new HomopolymerPair(new Homopolymer('C', 2), new Homopolymer('C', 2)),
                new HomopolymerPair(new Homopolymer('T', 1), new Homopolymer('T', 1)),
                new HomopolymerPair(new Homopolymer('C', 3), new Homopolymer('C', 2)),
                new HomopolymerPair(null, new Homopolymer('T', 1)),
                new HomopolymerPair(new Homopolymer('G', 2), new Homopolymer('G', 1)),
                new HomopolymerPair(new Homopolymer('A', 1), new Homopolymer('A', 1)));

        assertEquals(expectedPairs, actualPairs);
    }

    @Test
    public void testPairHomopolymersDeletionAndSnv()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer('C', 1),
                new Homopolymer('T', 1),
                new Homopolymer('C', 3),
                new Homopolymer('G', 2));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer('C', 1),
                new Homopolymer('T', 2),
                new Homopolymer('G', 1));

        List<HomopolymerPair> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
        List<HomopolymerPair> expectedPairs = Lists.newArrayList(
                new HomopolymerPair(new Homopolymer('C', 1), new Homopolymer('C', 1)),
                new HomopolymerPair(new Homopolymer('T', 1), new Homopolymer('T', 2)),
                new HomopolymerPair(new Homopolymer('C', 3), null),
                new HomopolymerPair(new Homopolymer('G', 2), new Homopolymer('G', 1)));

        assertEquals(expectedPairs, actualPairs);
    }
}
