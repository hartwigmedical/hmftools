package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.maskSandwichedSnvMnv;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Arrays;

import org.junit.Test;

import htsjdk.samtools.CigarOperator;

// TODO: re-align class name.
public class UltimaRealignedQualModelBuilderTest
{
    @Test
    public void testMaskSandwichedSnvMnvNoSandwichesSnv()
    {
        byte[] coreRefBases = new byte[] {5, 5, 5, 4, 4};
        byte[] coreReadBases = new byte[] {5, 5, 2, 4, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        assertArrayEquals(coreRefBases, coreRefBasesOut);
        assertArrayEquals(coreReadBases, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvNoSandwichesInsert()
    {
        byte[] coreRefBases = new byte[] {5, 5, 4, 4};
        byte[] coreReadBases = new byte[] {5, 5, 2, 4, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, I, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        assertArrayEquals(coreRefBases, coreRefBasesOut);
        assertArrayEquals(coreReadBases, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvNoSandwichesDelete()
    {
        byte[] coreRefBases = new byte[] {5, 5, 5, 4, 4};
        byte[] coreReadBases = new byte[] {5, 5, 4, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, D, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        assertArrayEquals(coreRefBases, coreRefBasesOut);
        assertArrayEquals(coreReadBases, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvSandwichedSnv()
    {
        byte[] coreRefBases = new byte[] {5, 5, 5, 5, 4};
        byte[] coreReadBases = new byte[] {5, 5, 2, 5, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        assertArrayEquals(coreRefBases, coreRefBasesOut);
        assertArrayEquals(new byte[] {5, 5, 5, 5, 4}, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvSandwichedMnv()
    {
        byte[] coreRefBases = new byte[] {5, 5, 2, 3, 2, 5, 4};
        byte[] coreReadBases = new byte[] {5, 5, 5, 5, 5, 5, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        assertArrayEquals(new byte[] {5, 5, 5, 5, 5, 5, 4}, coreRefBasesOut);
        assertArrayEquals(coreReadBases, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvSandwichedMultiple()
    {
        byte[] coreRefBases = new byte[] {5, 5, 5, 5, 5, 1, 2, 3, 5, 5};
        byte[] coreReadBases = new byte[] {5, 2, 5, 5, 5, 5, 5, 5, 5, 5};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        byte[] expectedOut = new byte[] {5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
        assertArrayEquals(expectedOut, coreRefBasesOut);
        assertArrayEquals(expectedOut, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvBreakOfRepeat()
    {
        byte[] coreRefBases = new byte[] {5, 5, 2, 3, 2, 5, 4, 4, 4};
        byte[] coreReadBases = new byte[] {5, 5, 5, 4, 5, 5, 4, 2, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        // TODO: created expectedFoo variables.
        assertArrayEquals(new byte[] {5, 5, 2, 3, 2, 5, 4, 4, 4}, coreRefBasesOut);
        assertArrayEquals(new byte[] {5, 5, 5, 4, 5, 5, 4, 4, 4}, coreReadBasesOut);

        // TODO: reverse case
        coreRefBasesOut = Arrays.copyArray(coreRefBases);
        coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreReadBasesOut, coreRefBasesOut);

        // TODO: created expectedFoo variables.
        assertArrayEquals(new byte[] {5, 5, 2, 3, 2, 5, 4, 4, 4}, coreRefBasesOut);
        assertArrayEquals(new byte[] {5, 5, 5, 4, 5, 5, 4, 4, 4}, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvSandwichedMissingBaseInRef()
    {
        byte[] coreRefBases = new byte[] {5, 5, 2, 2, 5, 4, 4, 4};
        byte[] coreReadBases = new byte[] {5, 5, 5, 5, 5, 5, 4, 2, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, I, M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        // TODO: created expectedFoo variables.
        assertArrayEquals(new byte[] {5, 5, 2, 2, 5, 4, 4, 4}, coreRefBasesOut);
        assertArrayEquals(new byte[] {5, 5, 5, 5, 5, 5, 4, 4, 4}, coreReadBasesOut);
    }

    @Test
    public void testMaskSandwichedSnvMnvSandwichedMissingBaseInRead()
    {
        byte[] coreRefBases = new byte[] {5, 5, 5, 5, 5, 5, 4, 2, 4};
        byte[] coreReadBases = new byte[] {5, 5, 2, 2, 5, 4, 4, 4};
        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, D, M, M, M, M, M);

        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
        maskSandwichedSnvMnv(coreCigarOps, coreRefBasesOut, coreReadBasesOut);

        // TODO: created expectedFoo variables.
        assertArrayEquals(new byte[] {5, 5, 5, 5, 5, 5, 4, 4, 4}, coreRefBasesOut);
        assertArrayEquals(new byte[] {5, 5, 2, 2, 5, 4, 4, 4}, coreReadBasesOut);
    }

    // TODO: re-enable tests.

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

    // TODO: Test inserts
    // TODO: test mnv variants.
    // TODO: test sandwiched snv/mnv variants.
}
