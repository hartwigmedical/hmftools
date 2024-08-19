package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.getHomopolymers;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.getQualVariants;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.getRealignedVariants;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.mergeSandwichedHomopolymers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.Homopolymer;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.MergedHomopolymers;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelBuilder.UltimaRealignedQualModel;

import org.junit.Test;

public class UltimaRealignedQualModelBuilderTest
{
    @Test
    public void testMergeSandwichedHomopolymersMultiHomopolymerSkip()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'C', 3),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'T', 10),
                new Homopolymer((byte) 'A', 1));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'T', 1),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'T', 10),
                new Homopolymer((byte) 'A', 1));

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(null, refHomopolymers, readHomopolymers);
        List<Homopolymer> expectedHomopolymers = refHomopolymers;

        assertEquals(expectedHomopolymers, mergedHomopolymers.RefHomopolymers);
        assertEquals(expectedHomopolymers, mergedHomopolymers.ReadHomopolymers);
    }

    @Test
    public void testMergeSandwichedHomopolymersValidationCases()
    {
        String refBases = "TCAAACAAACAAACAAACAAACAAAAAAAAAAAACT";
        String readBases = "TCAAACAAACAAACAAACAAACAAACAAAAAAACT";
        List<Homopolymer> refHomopolymers = getHomopolymers(refBases.getBytes(), 0, refBases.length() - 1);
        List<Homopolymer> readHomopolymers = getHomopolymers(readBases.getBytes(), 0, readBases.length() - 1);

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(null, refHomopolymers, readHomopolymers);
        List<Homopolymer> expectedMergedRefHomopolymers = refHomopolymers;
        List<Homopolymer> expectedMergedReadHomopolymers = Lists.newArrayList(refHomopolymers);
        expectedMergedReadHomopolymers.set(expectedMergedReadHomopolymers.size() - 3, new Homopolymer((byte) 'A', 11));

        assertEquals(expectedMergedRefHomopolymers, mergedHomopolymers.RefHomopolymers);
        assertEquals(expectedMergedReadHomopolymers, mergedHomopolymers.ReadHomopolymers);
    }

    private static class SandwichedHomopolymerTestCase
    {
        private final List<Homopolymer> RefHomopolymers;
        private final List<Homopolymer> ReadHomopolymers;
        private final List<Homopolymer> ExpectedHomopolymers;

        private SandwichedHomopolymerTestCase(final String refBases, final String readBases)
        {
            RefHomopolymers = getHomopolymers(refBases.getBytes(), 0, refBases.length() - 1);
            ReadHomopolymers = getHomopolymers(readBases.getBytes(), 0, readBases.length() - 1);
            ExpectedHomopolymers = Lists.newArrayList(new Homopolymer(refBases.getBytes()[0], refBases.length()));
        }

        public static SandwichedHomopolymerTestCase fromRefBases(final String refBases)
        {
            return new SandwichedHomopolymerTestCase(refBases, String.valueOf(refBases.charAt(0)).repeat(refBases.length()));
        }

        public static SandwichedHomopolymerTestCase fromReadBases(final String readBases)
        {
            return new SandwichedHomopolymerTestCase(String.valueOf(readBases.charAt(0)).repeat(readBases.length()), readBases);
        }

        public void check()
        {
            MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(null, RefHomopolymers, ReadHomopolymers);

            assertEquals(ExpectedHomopolymers, mergedHomopolymers.RefHomopolymers);
            assertEquals(ExpectedHomopolymers, mergedHomopolymers.ReadHomopolymers);
        }
    }

    @Test
    public void testSandwichedHomopolymer()
    {
        List<SandwichedHomopolymerTestCase> testCases = Lists.newArrayList(
                SandwichedHomopolymerTestCase.fromRefBases("ATA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATA"),
                SandwichedHomopolymerTestCase.fromRefBases("ATGA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATGA"),
                SandwichedHomopolymerTestCase.fromRefBases("ATTA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATTA"),
                SandwichedHomopolymerTestCase.fromRefBases("ATCGA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATCGA"),
                SandwichedHomopolymerTestCase.fromRefBases("ATTTA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATTTA"),
                SandwichedHomopolymerTestCase.fromRefBases("ATAGA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATAGA"),
                SandwichedHomopolymerTestCase.fromRefBases("ATATA"),
                SandwichedHomopolymerTestCase.fromReadBases("ATATA")
        );

        for(SandwichedHomopolymerTestCase testCase : testCases)
        {
            testCase.check();
        }
    }

    @Test
    public void testGetQualVariantsNonHomopolymerInsert()
    {
        int variantPos = 1000;
        SimpleVariant variant = new SimpleVariant(CHR_1, variantPos, "C", "CAT");
        SimpleVariant realignedVariant1 = new SimpleVariant(CHR_1, variantPos - 200, "CC", "C");
        SimpleVariant realignedVariant2 = new SimpleVariant(CHR_1, variantPos - 100, "C", "CCC");
        SimpleVariant realignedVariant3 = new SimpleVariant(CHR_1, variantPos, "C", "CA");
        SimpleVariant realignedVariant4 = new SimpleVariant(CHR_1, variantPos, "C", "CT");
        List<UltimaRealignedQualModel> realignedVariants = Lists.newArrayList(
                new UltimaRealignedQualModel(realignedVariant1),
                new UltimaRealignedQualModel(realignedVariant2),
                new UltimaRealignedQualModel(realignedVariant3),
                new UltimaRealignedQualModel(realignedVariant4));
        List<SimpleVariant> actualQualVariants = getQualVariants(false, variant, realignedVariants)
                .stream()
                .map(UltimaRealignedQualModel::variant)
                .collect(Collectors.toList());

        List<SimpleVariant> expectedQualVariants = Lists.newArrayList(realignedVariant2, realignedVariant3, realignedVariant4);

        assertEquals(expectedQualVariants.size(), actualQualVariants.size());
        for(int i = 0; i < expectedQualVariants.size(); i++)
        {
            assertTrue(expectedQualVariants.get(i).matches(actualQualVariants.get(i)));
        }
    }

    @Test
    public void testGetRealignedVariantsVariantStartForHomopolymerContraction()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'A', 5));
        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 2),
                new Homopolymer((byte) 'A', 5));
        int coreIndexStart = 0;
        int varIndex = 4;
        int corePositionStart = 100;
        int position = 104;
        SimpleVariant variant = new SimpleVariant(CHR_1, position,"ATTT", "A");
        VariantReadContext readContext = new VariantReadContext(variant, -1, -1, null, null, Lists.newArrayList(), coreIndexStart, varIndex, -1, null, null, null, corePositionStart, -1);
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(readContext, null, refHomopolymers, readHomopolymers);

        assertEquals(1, realignedVariants.size());
        assertTrue(variant.matches(realignedVariants.get(0).variant()));
    }

    @Test
    public void testGetRealignedVariantsVariantStartForHomopolymerExpansion()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 2),
                new Homopolymer((byte) 'A', 5));
        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'A', 5));
        int coreIndexStart = 0;
        int varIndex = 4;
        int corePositionStart = 100;
        int position = 104;
        SimpleVariant variant = new SimpleVariant(CHR_1, position,"A", "ATTT");
        VariantReadContext readContext = new VariantReadContext(variant, -1, -1, null, null, Lists.newArrayList(), coreIndexStart, varIndex, -1, null, null, null, corePositionStart, -1);
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(readContext, null, refHomopolymers, readHomopolymers);

        assertEquals(1, realignedVariants.size());
        assertTrue(variant.matches(realignedVariants.get(0).variant()));
    }

    // TODO: re-enable tests.

//    @Test
//    public void testMaskSandwichedSnvMnvNoSandwichesSnv()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 5, 4, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 2, 4, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(coreRefBases, coreRefBasesOut);
//        assertArrayEquals(coreReadBases, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvNoSandwichesInsert()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 4, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 2, 4, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, I, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(coreRefBases, coreRefBasesOut);
//        assertArrayEquals(coreReadBases, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvNoSandwichesDelete()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 5, 4, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 4, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, D, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(coreRefBases, coreRefBasesOut);
//        assertArrayEquals(coreReadBases, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvSandwichedSnv()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 5, 5, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 2, 5, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(coreRefBases, coreRefBasesOut);
//        assertArrayEquals(new byte[] {5, 5, 5, 5, 4}, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvSandwichedMnv()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 2, 3, 2, 5, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 5, 5, 5, 5, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(new byte[] {5, 5, 5, 5, 5, 5, 4}, coreRefBasesOut);
//        assertArrayEquals(coreReadBases, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvSandwichedMultiple()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 5, 5, 5, 1, 2, 3, 5, 5};
//        byte[] coreReadBases = new byte[] {5, 2, 5, 5, 5, 5, 5, 5, 5, 5};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        byte[] expectedOut = new byte[] {5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
//        assertArrayEquals(expectedOut, coreRefBasesOut);
//        assertArrayEquals(expectedOut, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvBreakOfRepeat()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 2, 3, 2, 5, 4, 4, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 5, 4, 5, 5, 4, 2, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        // TODO: created expectedFoo variables.
//        assertArrayEquals(new byte[] {5, 5, 2, 3, 2, 5, 4, 4, 4}, coreRefBasesOut);
//        assertArrayEquals(new byte[] {5, 5, 5, 4, 5, 5, 4, 4, 4}, coreReadBasesOut);
//
//        // TODO: reverse case
//        coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreReadBasesOut, coreRefBasesOut));
//
//        // TODO: created expectedFoo variables.
//        assertArrayEquals(new byte[] {5, 5, 2, 3, 2, 5, 4, 4, 4}, coreRefBasesOut);
//        assertArrayEquals(new byte[] {5, 5, 5, 4, 5, 5, 4, 4, 4}, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvSandwichedMissingBaseInRef()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 2, 2, 5, 4, 4, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 5, 5, 5, 5, 4, 2, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, I, M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        // TODO: created expectedFoo variables.
//        assertArrayEquals(new byte[] {5, 5, 2, 2, 5, 4, 4, 4}, coreRefBasesOut);
//        assertArrayEquals(new byte[] {5, 5, 5, 5, 5, 5, 4, 4, 4}, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvSandwichedMissingBaseInRead()
//    {
//        SimpleVariant variant = new SimpleVariant(CHR_1, 1, "A", "AA");
//        VariantReadContext readContext = new VariantReadContext(variant, 0, 0, null, (byte) 0, null, Lists.newArrayList(), 0, 0, 0, null, null, null, 0, 0);
//
//        byte[] coreRefBases = new byte[] {5, 5, 5, 5, 5, 5, 4, 2, 4};
//        byte[] coreReadBases = new byte[] {5, 5, 2, 2, 5, 4, 4, 4};
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, D, M, M, M, M, M);
//
//        byte[] coreRefBasesOut = Arrays.copyArray(coreRefBases);
//        byte[] coreReadBasesOut = Arrays.copyArray(coreReadBases);
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        // TODO: created expectedFoo variables.
//        assertArrayEquals(new byte[] {5, 5, 5, 5, 5, 5, 4, 4, 4}, coreRefBasesOut);
//        assertArrayEquals(new byte[] {5, 5, 2, 2, 5, 4, 4, 4}, coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvVariantIsNotSandwichedSnv()
//    {
//        String coreRef = "AAATT";
//        String coreRead = "AAGTT";
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M);
//
//        SimpleVariant variant = new SimpleVariant(CHR_1, 102, "A", "G");
//        VariantReadContext readContext = new VariantReadContext(variant, 1, 200, null, (byte) 0, null, Lists.newArrayList(), 0, 2, 4, null, null, null, 100, 104);
//
//        byte[] coreRefBasesOut = coreRef.getBytes();
//        byte[] coreReadBasesOut = coreRead.getBytes();
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(coreRef.getBytes(), coreRefBasesOut);
//        assertArrayEquals(coreRead.getBytes(), coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvVariantIsSandwichedSnv()
//    {
//        String coreRef = "AAAAA";
//        String coreRead = "AAGAA";
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M);
//
//        SimpleVariant variant = new SimpleVariant(CHR_1, 102, "A", "G");
//        VariantReadContext readContext = new VariantReadContext(variant, 1, 200, null, (byte) 0, null, Lists.newArrayList(), 0, 2, 4, null, null, null, 100, 104);
//
//        byte[] coreRefBasesOut = coreRef.getBytes();
//        byte[] coreReadBasesOut = coreRead.getBytes();
//        assertTrue(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        // TODO: expected output variables.
//        assertArrayEquals(coreRef.getBytes(), coreRefBasesOut);
//        assertArrayEquals(coreRef.getBytes(), coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvVariantIsNotSandwichedMnv()
//    {
//        String coreRef = "AAACATT";
//        String coreRead = "AAGTGTT";
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M);
//
//        SimpleVariant variant = new SimpleVariant(CHR_1, 102, "ACA", "GTG");
//        VariantReadContext readContext = new VariantReadContext(variant, 1, 200, null, (byte) 0, null, Lists.newArrayList(), 0, 2, 6, null, null, null, 100, 106);
//
//        byte[] coreRefBasesOut = coreRef.getBytes();
//        byte[] coreReadBasesOut = coreRead.getBytes();
//        assertFalse(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        assertArrayEquals(coreRef.getBytes(), coreRefBasesOut);
//        assertArrayEquals(coreRead.getBytes(), coreReadBasesOut);
//    }
//
//    @Test
//    public void testMaskSandwichedSnvMnvVariantIsSandwichedMnv()
//    {
//        String coreRef = "AAGTGAA";
//        String coreRead = "AAAAAAA";
//        List<CigarOperator> coreCigarOps = Lists.newArrayList(M, M, M, M, M, M, M);
//
//        SimpleVariant variant = new SimpleVariant(CHR_1, 102, "GTG","ACA");
//        VariantReadContext readContext = new VariantReadContext(variant, 1, 200, null, (byte) 0, null, Lists.newArrayList(), 0, 2, 6, null, null, null, 100, 106);
//
//        byte[] coreRefBasesOut = coreRef.getBytes();
//        byte[] coreReadBasesOut = coreRead.getBytes();
//        assertTrue(maskSandwichedSnvMnv(readContext, coreCigarOps, coreRefBasesOut, coreReadBasesOut));
//
//        // TODO: expected output variables.
//        assertArrayEquals(coreRead.getBytes(), coreRefBasesOut);
//        assertArrayEquals(coreRead.getBytes(), coreReadBasesOut);
//    }

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
