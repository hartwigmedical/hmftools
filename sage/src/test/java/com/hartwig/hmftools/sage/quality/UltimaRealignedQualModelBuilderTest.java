package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.HOMOPOLYMER_ADJUSTMENT;
import static com.hartwig.hmftools.sage.quality.UltimaModelType.HOMOPOLYMER_DELETION;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.getHomopolymers;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.getQualVariants;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.getRealignedVariants;
import static com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.mergeSandwichedHomopolymers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.UltimaQualCalculator.HomopolymerAdjustment;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.Homopolymer;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.MergedHomopolymers;
import com.hartwig.hmftools.sage.quality.UltimaRealignedQualModelsBuilder.RefMask;

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

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(0);

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(mockReadContext, refHomopolymers, readHomopolymers, true);
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

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(0);

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(mockReadContext, refHomopolymers, readHomopolymers, true);
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
            SimpleVariant mockVariant = mock(SimpleVariant.class);
            when(mockVariant.isIndel()).thenReturn(true);

            VariantReadContext mockReadContext = mock(VariantReadContext.class);
            when(mockReadContext.variant()).thenReturn(mockVariant);
            when(mockReadContext.corePositionStart()).thenReturn(0);

            MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(mockReadContext, RefHomopolymers, ReadHomopolymers, true);

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
            assertTrue(expectedQualVariants.get(i).matches(actualQualVariants.get(i)));
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
        VariantReadContext readContext = new VariantReadContext(variant, -1, -1, null, "AAAAATTAAAAA".getBytes(), Lists.newArrayList(), coreIndexStart, varIndex, -1, null, null, null, corePositionStart, -1);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(readContext, null, refHomopolymers, readHomopolymers, refMasks);

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
        VariantReadContext readContext = new VariantReadContext(variant, -1, -1, null, "AAAAATTTTTAAAAA".getBytes(), Lists.newArrayList(), coreIndexStart, varIndex, -1, null, null, null, corePositionStart, -1);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(readContext, null, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());
        assertTrue(variant.matches(realignedVariants.get(0).variant()));
    }

    @Test
    public void testSandwichRefGenome()
    {
        String refBases = "GTAAAAAAAAAAAAGAGAGAAAAAAC";
        String readCoreBases = "GTGAAAAAAAAAAGAGAGAAAAAAC";
        List<Homopolymer> refHomopolymers = getHomopolymers(refBases.getBytes(), 0, refBases.length() - 1);
        List<Homopolymer> readHomopolymers = getHomopolymers(readCoreBases.getBytes(), 0, readCoreBases.length() - 1);

        int corePositionStart = 13;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(mockReadContext, refHomopolymers, readHomopolymers, false);

        int refMaskIndex = 18;
        assertEquals(1, mergedHomopolymers.refMasks().size());

        RefMask refMask = mergedHomopolymers.refMasks().get(0);
        assertEquals((byte) 'A', refMask.BaseMask);
        assertEquals(1, refMask.PosEnd - refMask.PosStart + 1);
        assertEquals(corePositionStart + refMaskIndex, refMask.PosStart);
    }

    @Test
    public void testGetRealignedVariantsSingleHomopolymerDeletion()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 2),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));
    }

    @Test
    public void testGetRealignedVariantsSingleHomopolymerContraction()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        HomopolymerAdjustment qualModel = (HomopolymerAdjustment) realignedVariant.baseQualModel();

        assertEquals(2, qualModel.refAdjustCount());
        assertEquals(1, qualModel.hpStartIndex() + realignedVariant.varReadIndexOffset());
        assertEquals(3, qualModel.hpEndIndex() + realignedVariant.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsTwoHomopolymerDeletions()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 2),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(0, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 6, "TC", "A"));
    }

    @Test
    public void testGetRealignedVariantsTwoHomopolymerContractions()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        HomopolymerAdjustment qualModel1 = (HomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(3, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(3, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 9, "TC", "T"));

        HomopolymerAdjustment qualModel2 = (HomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(1, qualModel2.refAdjustCount());
        assertEquals(4, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(4, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsHomopolyerDeletionThenHomopolymerContraction()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATTTTT", "A"));

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 9, "TC", "A"));

        HomopolymerAdjustment qualModel2 = (HomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(1, qualModel2.refAdjustCount());
        assertEquals(1, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(1, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsHomopolyerContractionThenHomopolymerDeletion()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        HomopolymerAdjustment qualModel1 = (HomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(3, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(3, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 9, "TCC", "T"));
    }

    @Test
    public void testGetRealignedVariantsSingleHomopolymerCreation()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 2),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        HomopolymerAdjustment qualModel = (HomopolymerAdjustment) realignedVariant.baseQualModel();

        assertEquals(-2, qualModel.refAdjustCount());
        assertEquals(1, qualModel.hpStartIndex() + realignedVariant.varReadIndexOffset());
        assertEquals(2, qualModel.hpEndIndex() + realignedVariant.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsSingleHomopolymerExpansion()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        HomopolymerAdjustment qualModel = (HomopolymerAdjustment) realignedVariant.baseQualModel();

        assertEquals(-2, qualModel.refAdjustCount());
        assertEquals(1, qualModel.hpStartIndex() + realignedVariant.varReadIndexOffset());
        assertEquals(5, qualModel.hpEndIndex() + realignedVariant.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsTwoHomopolymerCreations()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 2),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        HomopolymerAdjustment qualModel1 = (HomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(2, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(2, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 4, "A", "TC"));

        HomopolymerAdjustment qualModel2 = (HomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(-1, qualModel2.refAdjustCount());
        assertEquals(3, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(3, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsTwoHomopolymerExpansions()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        HomopolymerAdjustment qualModel1 = (HomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(5, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(5, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 7, "T", "TC"));

        HomopolymerAdjustment qualModel2 = (HomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(-1, qualModel2.refAdjustCount());
        assertEquals(6, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(7, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsHomopolyerCreationThenHomopolymerExpansion()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATTTTT"));

        HomopolymerAdjustment qualModel1 = (HomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-5, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(5, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(5, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 4, "A", "TC"));

        HomopolymerAdjustment qualModel2 = (HomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(-1, qualModel2.refAdjustCount());
        assertEquals(6, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(7, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsHomopolyerExpansionThenHomopolymerCreation()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 5),
                new Homopolymer((byte) 'C', 2),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        HomopolymerAdjustment qualModel1 = (HomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(5, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(5, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 7, "T", "TCC"));

        HomopolymerAdjustment qualModel2 = (HomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(-2, qualModel2.refAdjustCount());
        assertEquals(6, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(7, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsSNV()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 1),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'C', 1),
                new Homopolymer((byte) 'G', 5));

        int refGenomePaddingLength = 99;
        String refGenomePadding = "C".repeat(refGenomePaddingLength);
        String refGenomeBases = refGenomePadding + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + refGenomePadding;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, refGenomeBases);
        refGenome.ChromosomeLengths.put(CHR_1, refGenomeBases.length());

        int flankLength = 10;
        String flankBases = "C".repeat(flankLength);
        int corePositionStart = refGenomePaddingLength + 1;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 4;
        String readBases = flankBases + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + flankBases;

        SimpleVariant mockVariant = mock(SimpleVariant.class);
        when(mockVariant.chromosome()).thenReturn(CHR_1);
        when(mockVariant.isIndel()).thenReturn(true);

        VariantReadContext mockReadContext = mock(VariantReadContext.class);
        when(mockReadContext.variant()).thenReturn(mockVariant);
        when(mockReadContext.corePositionStart()).thenReturn(corePositionStart);
        when(mockReadContext.coreIndexStart()).thenReturn(coreIndexStart);
        when(mockReadContext.varIndex()).thenReturn(varIndex);
        when(mockReadContext.readBasesBytes()).thenReturn(readBases.getBytes());

        UltimaQualCalculator ultimaQualCalculator = new UltimaQualCalculator(refGenome);
        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(mockReadContext, ultimaQualCalculator, refHomopolymers, readHomopolymers, refMasks);
        List<UltimaRealignedQualModel> delRealignedVariants = realignedVariants.stream().filter(x -> x.type() == HOMOPOLYMER_DELETION).collect(Collectors.toList());
        List<UltimaRealignedQualModel> adjRealignedVariants = realignedVariants.stream().filter(x -> x.type() == HOMOPOLYMER_ADJUSTMENT).collect(Collectors.toList());

        assertEquals(1, delRealignedVariants.size());
        assertEquals(1, adjRealignedVariants.size());

        UltimaRealignedQualModel delRealignedVariant = delRealignedVariants.get(0);
        UltimaRealignedQualModel adjRealignedVariant = adjRealignedVariants.get(0);

        assertTrue(delRealignedVariant.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(0, delRealignedVariant.varReadIndexOffset());
        assertTrue(delRealignedVariant.variant().matches(CHR_1, corePositionStart + 4, "AT", "A"));

        assertTrue(adjRealignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, adjRealignedVariant.varReadIndexOffset());
        assertTrue(adjRealignedVariant.variant().matches(CHR_1, corePositionStart + 4, "A", "AC"));

        HomopolymerAdjustment adjQualModel = (HomopolymerAdjustment) adjRealignedVariant.baseQualModel();

        assertEquals(-1, adjQualModel.refAdjustCount());
        assertEquals(1, adjQualModel.hpStartIndex() + adjRealignedVariant.varReadIndexOffset());
        assertEquals(1, adjQualModel.hpEndIndex() + adjRealignedVariant.varReadIndexOffset());
    }
}
