package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.pairHomopolymers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.Homopolymer;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.HomopolymerIndel;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.HomopolymerMatch;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner.HomopolymerPair;

import org.junit.Test;

public class UltimaLocalRealignerTest
{
    // TODO: LATER re-enable tests.

//    @Test
//    public void testPairHomopolymersEmpty()
//    {
//        List<Homopolymer> refHomopolymers = Lists.newArrayList();
//        List<Homopolymer> readHomopolymers = Lists.newArrayList();
//
//        Set<List<HomopolymerPair>> pairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//
//        assertTrue(pairs.isEmpty());
//    }
//
//    @Test
//    public void testPairHomopolymersSingleMatch()
//    {
//        Homopolymer a1 = new Homopolymer('A', 1);
//        Homopolymer a3 = new Homopolymer('A', 3);
//
//        List<Homopolymer> refHomopolymers = Lists.newArrayList(a1);
//        List<Homopolymer> readHomopolymers = Lists.newArrayList(a3);
//
//        Set<List<HomopolymerPair>> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//        Set<List<HomopolymerPair>> expectedPairs = Set.of(Lists.newArrayList(new HomopolymerMatch(a1, a3)));
//
//        assertEquals(expectedPairs, actualPairs);
//    }
//
//    @Test
//    public void testPairHomopolymersSingleMatchmatch()
//    {
//        Homopolymer a1 = new Homopolymer('A', 1);
//        Homopolymer t1 = new Homopolymer('T', 1);
//
//        List<Homopolymer> refHomopolymers = Lists.newArrayList(a1);
//        List<Homopolymer> readHomopolymers = Lists.newArrayList(t1);
//
//        Set<List<HomopolymerPair>> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//        Set<List<HomopolymerPair>> expectedPairs = Set.of(Lists.newArrayList(new HomopolymerIndel(List.of(a1), List.of(t1))));
//
//        assertEquals(expectedPairs, actualPairs);
//    }
//
//    @Test
//    public void testPairHomopolymersMultipleMatch()
//    {
//        Homopolymer a1 = new Homopolymer('A', 1);
//        Homopolymer a3 = new Homopolymer('A', 3);
//        Homopolymer t1 = new Homopolymer('T', 1);
//        Homopolymer t2 = new Homopolymer('T', 2);
//
//        List<Homopolymer> refHomopolymers = Lists.newArrayList(a1, t2);
//        List<Homopolymer> readHomopolymers = Lists.newArrayList(a3, t1);
//
//        Set<List<HomopolymerPair>> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//        Set<List<HomopolymerPair>> expectedPairs = Set.of(Lists.newArrayList(new HomopolymerMatch(a1, a3), new HomopolymerMatch(t2, t1)));
//
//        assertEquals(expectedPairs, actualPairs);
//    }
//
//    @Test
//    public void testPairHomopolymersContractionAndSnv()
//    {
//        List<Homopolymer> refHomopolymers = Lists.newArrayList(
//                new Homopolymer('C', 2),
//                new Homopolymer('T', 1),
//                new Homopolymer('C', 3),
//                new Homopolymer('G', 2),
//                new Homopolymer('A', 1));
//
//        List<Homopolymer> readHomopolymers = Lists.newArrayList(
//                new Homopolymer('C', 2),
//                new Homopolymer('T', 1),
//                new Homopolymer('C', 2),
//                new Homopolymer('T', 1),
//                new Homopolymer('G', 1),
//                new Homopolymer('A', 1));
//
//        Set<List<HomopolymerPair>> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//        Set<List<HomopolymerPair>> expectedPairs = Set.of(Lists.newArrayList(
//                new HomopolymerMatch(new Homopolymer('C', 2), new Homopolymer('C', 2)),
//                new HomopolymerMatch(new Homopolymer('T', 1), new Homopolymer('T', 1)),
//                new HomopolymerMatch(new Homopolymer('C', 3), new Homopolymer('C', 2)),
//                new HomopolymerIndel(Lists.newArrayList(), Lists.newArrayList(new Homopolymer('T', 1))),
//                new HomopolymerMatch(new Homopolymer('G', 2), new Homopolymer('G', 1)),
//                new HomopolymerMatch(new Homopolymer('A', 1), new Homopolymer('A', 1))));
//
//        assertEquals(expectedPairs, actualPairs);
//    }
//
//    @Test
//    public void testPairHomopolymersDeletionAndSnv()
//    {
//        List<Homopolymer> refHomopolymers = Lists.newArrayList(
//                new Homopolymer('C', 1),
//                new Homopolymer('T', 1),
//                new Homopolymer('C', 3),
//                new Homopolymer('G', 2));
//
//        List<Homopolymer> readHomopolymers = Lists.newArrayList(
//                new Homopolymer('C', 1),
//                new Homopolymer('T', 2),
//                new Homopolymer('G', 1));
//
//        Set<List<HomopolymerPair>> actualPairs = pairHomopolymers(refHomopolymers, readHomopolymers);
//        Set<List<HomopolymerPair>> expectedPairs = Set.of(Lists.newArrayList(
//                new HomopolymerMatch(new Homopolymer('C', 1), new Homopolymer('C', 1)),
//                new HomopolymerMatch(new Homopolymer('T', 1), new Homopolymer('T', 2)),
//                new HomopolymerIndel(Lists.newArrayList(new Homopolymer('C', 3)), Lists.newArrayList()),
//                new HomopolymerMatch(new Homopolymer('G', 2), new Homopolymer('G', 1))));
//
//        assertEquals(expectedPairs, actualPairs);
//    }
}
