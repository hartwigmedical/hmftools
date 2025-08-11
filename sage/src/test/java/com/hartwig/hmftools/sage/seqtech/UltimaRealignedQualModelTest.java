package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.sage.common.TestUtils.setUltimaSequencing;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.seqtech.Homopolymer.getHomopolymers;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.HOMOPOLYMER_ADJUSTMENT;
import static com.hartwig.hmftools.sage.seqtech.UltimaModelType.HOMOPOLYMER_DELETION;
import static com.hartwig.hmftools.sage.seqtech.UltimaRealignedQualModelsBuilder.getRealignedVariants;
import static com.hartwig.hmftools.sage.seqtech.UltimaRealignedQualModelsBuilder.mergeSandwichedHomopolymers;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.isCleanSnv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.CigarElement;

public class UltimaRealignedQualModelTest
{
    public UltimaRealignedQualModelTest()
    {
        setUltimaSequencing();
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    @Test
    public void testGetHomopolymersInvalidIndices()
    {
        byte[] bases = new byte[] { (byte) 'A', (byte) 'A' };

        List<Homopolymer> homopolymers = getHomopolymers(bases, -1, 1);
        assertEquals(0, homopolymers.size());

        homopolymers = getHomopolymers(bases, 0, 2);
        assertEquals(0, homopolymers.size());

        homopolymers = getHomopolymers(bases, 1, 0);
        assertEquals(0, homopolymers.size());
    }

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

        SimpleVariant variant = new SimpleVariant(CHR_1, 100, "A", "AC");
        VariantReadContext readContext = createReadContext(variant);

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(readContext, refHomopolymers, readHomopolymers);
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

        SimpleVariant variant = new SimpleVariant(CHR_1, 100, "A", "AC");
        VariantReadContext readContext = createReadContext(variant);

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(readContext, refHomopolymers, readHomopolymers);
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
            SimpleVariant variant = new SimpleVariant(CHR_1, 50, "A", "AA");
            String refBases = REF_BASES_200.substring(30, 70);
            String readBases = REF_BASES_200.substring(30, 51) + "A" + REF_BASES_200.substring(51, 70);
            VariantReadContext readContext = createVariantReadContext(variant, readBases.getBytes(), refBases.getBytes());

            MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(readContext, RefHomopolymers, ReadHomopolymers);

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

        SimpleVariant variant = new SimpleVariant(CHR_1, position, "ATTT", "A");

        VariantReadContext readContext = new VariantReadContext(
                variant, -1, -1, null, "AAAAATTAAAAA".getBytes(), Lists.newArrayList(),
                coreIndexStart, varIndex, -1, null, null, null, corePositionStart, -1);

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, null, refHomopolymers, readHomopolymers, refMasks);

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

        SimpleVariant variant = new SimpleVariant(CHR_1, position, "A", "ATTT");
        VariantReadContext readContext = new VariantReadContext(
                variant, -1, -1, null, "AAAAATTTTTAAAAA".getBytes(), Lists.newArrayList(),
                coreIndexStart, varIndex, -1, null, null, null, corePositionStart, -1);

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, null, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());
        assertTrue(variant.matches(realignedVariants.get(0).variant()));
    }

    @Test
    public void testSandwichRefGenome()
    {
        String refBases =      "GTAAAAAAAAAAAAGAGAGAAAAAAC";
        String readCoreBases = "GTGAAAAAAAAAAGAGAGAAAAAAC";
        List<Homopolymer> refHomopolymers = getHomopolymers(refBases.getBytes(), 0, refBases.length() - 1);
        List<Homopolymer> readHomopolymers = getHomopolymers(readCoreBases.getBytes(), 0, readCoreBases.length() - 1);

        int corePositionStart = 13;

        SimpleVariant variant = new SimpleVariant(CHR_1, 15, "AA", "G");

        VariantReadContext readContext = createVariantReadContext(variant, readCoreBases.getBytes(), refBases.getBytes());

        MergedHomopolymers mergedHomopolymers = mergeSandwichedHomopolymers(readContext, refHomopolymers, readHomopolymers);

        int refMaskIndex = 18;
        assertEquals(1, mergedHomopolymers.refMasks().size());

        RefMask refMask = mergedHomopolymers.refMasks().get(0);
        assertEquals((byte) 'A', refMask.BaseMask);
        assertEquals(1, refMask.PosEnd - refMask.PosStart + 1);
        assertEquals(corePositionStart + refMaskIndex, refMask.PosStart);
    }

    // common state for following realigned qual model tests
    private static final String FLANK_BASES = "C".repeat(DEFAULT_FLANK_LENGTH);
    private static final String REF_BASE_PADDING = "C".repeat(99);
    private static final int CORE_INDEX_START = 10;
    private static final int CORE_POS_START = 100;
    private static final int VAR_CORE_START_OFFSET = 4;

    private static MockRefGenome REF_GENOME = new MockRefGenome(true);

    private static void populateRefGenomeBases(final List<Homopolymer> refHomopolymers)
    {
        String refBases = REF_BASE_PADDING + refHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + REF_BASE_PADDING;

        REF_GENOME.RefGenomeMap.put(CHR_1, refBases);
        REF_GENOME.ChromosomeLengths.put(CHR_1, refBases.length());
    }

    private static final UltimaQualCalculator ULTIMA_QUAL_CALCULATOR = new UltimaQualCalculator(REF_GENOME);

    private static VariantReadContext createTestReadContext(final List<Homopolymer> readHomopolymers)
    {
        int coreIndexStart = CORE_INDEX_START;
        int corePositionStart = CORE_POS_START;
        int varPosition = corePositionStart + VAR_CORE_START_OFFSET;
        int varIndex = CORE_INDEX_START + VAR_CORE_START_OFFSET;
        String readBases = FLANK_BASES + readHomopolymers.stream().map(Homopolymer::expand).collect(Collectors.joining()) + FLANK_BASES;

        SimpleVariant variant = new SimpleVariant(CHR_1, varPosition, "AA", "G");

        VariantReadContext readContext = createVariantReadContext(
                variant, readBases.getBytes(), readBases.getBytes(), coreIndexStart, varIndex);

        return readContext;
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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        UltimaHomopolymerAdjustment qualModel = (UltimaHomopolymerAdjustment) realignedVariant.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(3, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(3, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 9, "TC", "T"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_DELETION);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATTTTT", "A"));

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 9, "TC", "A"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "ATT", "A"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        UltimaHomopolymerAdjustment qualModel = (UltimaHomopolymerAdjustment) realignedVariant.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(1, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant = realignedVariants.get(0);

        assertTrue(realignedVariant.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant.varReadIndexOffset());
        assertTrue(realignedVariant.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        UltimaHomopolymerAdjustment qualModel = (UltimaHomopolymerAdjustment) realignedVariant.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(2, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(2, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 4, "A", "TC"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(5, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(5, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 7, "T", "TC"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATTTTT"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-5, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(5, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(5, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 4, "A", "TC"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);

        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "ATT"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-2, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(5, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(5, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 7, "T", "TCC"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(-2, qualModel2.refAdjustCount());
        assertEquals(6, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(7, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
    }

    @Test
    public void testGetRealignedVariantsHomopolyerCreationThenHomopolymerCreation()
    {
        List<Homopolymer> refHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'G', 5));

        List<Homopolymer> readHomopolymers = Lists.newArrayList(
                new Homopolymer((byte) 'A', 5),
                new Homopolymer((byte) 'G', 1),
                new Homopolymer((byte) 'A', 1),
                new Homopolymer((byte) 'T', 3),
                new Homopolymer((byte) 'G', 5));

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        assertEquals(2, realignedVariants.size());

        UltimaRealignedQualModel realignedVariant1 = realignedVariants.get(0);
        UltimaRealignedQualModel realignedVariant2 = realignedVariants.get(1);
        assertTrue(realignedVariant1.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(0, realignedVariant1.varReadIndexOffset());
        assertTrue(realignedVariant1.variant().matches(CHR_1, corePositionStart + 4, "A", "AG"));

        UltimaHomopolymerAdjustment qualModel1 = (UltimaHomopolymerAdjustment) realignedVariant1.baseQualModel();

        assertEquals(-1, qualModel1.refAdjustCount());
        assertEquals(1, qualModel1.hpStartIndex() + realignedVariant1.varReadIndexOffset());
        assertEquals(1, qualModel1.hpEndIndex() + realignedVariant1.varReadIndexOffset());

        assertTrue(realignedVariant2.baseQualModel().type() == HOMOPOLYMER_ADJUSTMENT);
        assertEquals(1, realignedVariant2.varReadIndexOffset());
        assertTrue(realignedVariant2.variant().matches(CHR_1, corePositionStart + 4, "A", "GA"));

        UltimaHomopolymerAdjustment qualModel2 = (UltimaHomopolymerAdjustment) realignedVariant2.baseQualModel();

        assertEquals(-1, qualModel2.refAdjustCount());
        assertEquals(2, qualModel2.hpStartIndex() + realignedVariant2.varReadIndexOffset());
        assertEquals(2, qualModel2.hpEndIndex() + realignedVariant2.varReadIndexOffset());
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

        populateRefGenomeBases(refHomopolymers);

        VariantReadContext readContext = createTestReadContext(readHomopolymers);
        int corePositionStart = CORE_POS_START;

        List<RefMask> refMasks = Lists.newArrayList();
        List<UltimaRealignedQualModel> realignedVariants = getRealignedVariants(
                readContext, ULTIMA_QUAL_CALCULATOR, refHomopolymers, readHomopolymers, refMasks);

        List<UltimaRealignedQualModel> delRealignedVariants = realignedVariants.stream()
                .filter(x -> x.type() == HOMOPOLYMER_DELETION)
                .collect(Collectors.toList());

        List<UltimaRealignedQualModel> adjRealignedVariants = realignedVariants.stream()
                .filter(x -> x.type() == HOMOPOLYMER_ADJUSTMENT)
                .collect(Collectors.toList());

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

        UltimaHomopolymerAdjustment adjQualModel = (UltimaHomopolymerAdjustment) adjRealignedVariant.baseQualModel();

        assertEquals(-1, adjQualModel.refAdjustCount());
        assertEquals(1, adjQualModel.hpStartIndex() + adjRealignedVariant.varReadIndexOffset());
        assertEquals(1, adjQualModel.hpEndIndex() + adjRealignedVariant.varReadIndexOffset());
    }

    @Test
    public void testIsCleanSnvIndel()
    {
        SimpleVariant variant = new SimpleVariant(CHR_1, 100, "A", "AC");

        String refBases = "AGACT";
        String readBases = "AGACCT";

        VariantReadContext readContext = createVariantReadContext(variant, readBases.getBytes(), refBases.getBytes());
        assertFalse(isCleanSnv(readContext));
    }

    @Test
    public void testIsCleanSnvRefCoreLengthMismatch()
    {
        String readBases = "AATAAA";
        String refBases = "AAAAAA";

        VariantReadContext readContext = createVariantReadContext(createSimpleSnv(), readBases.getBytes(), refBases.getBytes());

        assertFalse(isCleanSnv(readContext));

        assertEquals(5, readContext.coreLength());
        assertFalse(isCleanSnv(readContext));
    }

    @Test
    public void testIsCleanSnvClean()
    {
        int coreLength = 5;
        int flankLength = 5;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 1;
        int coreIndexEnd = coreIndexStart + coreLength - 1;
        String readBases = "A".repeat(coreLength + 2 * flankLength);

        StringBuilder refBasesBuilder = new StringBuilder(readBases);
        refBasesBuilder.setCharAt(varIndex, 'T');
        refBasesBuilder.delete(coreIndexEnd + 1, refBasesBuilder.length());
        refBasesBuilder.delete(0, coreIndexStart);
        String refBases = refBasesBuilder.toString();

        VariantReadContext readContext = createVariantReadContext(createSimpleSnv(), readBases.getBytes(), refBases.getBytes());

        assertFalse(isCleanSnv(readContext));
    }

    @Test
    public void testIsCleanSnvNotClean()
    {
        int coreLength = 5;
        int flankLength = 5;
        int coreIndexStart = flankLength;
        int varIndex = flankLength + 1;
        int coreIndexEnd = coreIndexStart + coreLength - 1;
        String readBases = "A".repeat(coreLength + 2 * flankLength);

        StringBuilder refBasesBuilder = new StringBuilder(readBases);
        refBasesBuilder.setCharAt(varIndex, 'T');
        refBasesBuilder.setCharAt(varIndex + 1, 'T');
        refBasesBuilder.delete(coreIndexEnd + 1, refBasesBuilder.length());
        refBasesBuilder.delete(0, coreIndexStart);
        String refBases = refBasesBuilder.toString();

        VariantReadContext readContext = createVariantReadContext(createSimpleSnv(), readBases.getBytes(), refBases.getBytes());

        assertFalse(isCleanSnv(readContext));
    }

    private static SimpleVariant createSimpleSnv()
    {
        return new SimpleVariant(CHR_1, 100, "A", "C");
    }

    private static VariantReadContext createVariantReadContext(final SimpleVariant variant, final byte[] readBases, final byte[] refBases)
    {
        int flankLength = DEFAULT_FLANK_LENGTH;
        int coreIndexStart = flankLength;
        int varIndex = coreIndexStart + 2;

        return createVariantReadContext(variant, readBases, refBases, coreIndexStart, varIndex);
    }
    private static VariantReadContext createVariantReadContext(
            final SimpleVariant variant, final byte[] readBases, final byte[] refBases, int coreIndexStart, int varIndex)
    {
        List<CigarElement> readCigar = CigarUtils.cigarElementsFromStr("25M");

        int corePosStart = variant.Position - (varIndex - coreIndexStart); // assumes no other indels

        int corePosEnd = variant.Position + variant.indelLength() + 2;

        int coreIndexEnd = varIndex + variant.indelLength() + 2;
        int alignStart = variant.Position - varIndex;
        int alignEnd = variant.Position + (coreIndexEnd - varIndex) + DEFAULT_FLANK_LENGTH;

        return new VariantReadContext(
                variant, alignStart, alignEnd, refBases, readBases, readCigar, coreIndexStart, varIndex, coreIndexEnd,
                null, null, Collections.emptyList(), corePosStart, corePosEnd);
    }

}
