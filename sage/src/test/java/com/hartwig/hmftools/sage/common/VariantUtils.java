package com.hartwig.hmftools.sage.common;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.QUALITY_CALCULATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class VariantUtils
{
    // Simple variant creation
    public static SimpleVariant createSimpleVariant(int position)
    {
        return new SimpleVariant(CHR_1, position, "A", "C");
    }

    public static SimpleVariant createSimpleVariant(int position, final String ref, final String alt)
    {
        return new SimpleVariant(CHR_1, position, ref, alt);
    }

    public static SimpleVariant createSimpleVariant(final String chromosome, int position, final String ref, final String alt)
    {
        return createSimpleVariant(chromosome, position, ref, alt);
    }

    // Variant read context creation               0123456789
    public static final String TEST_LEFT_FLANK  = "ACCGCTGACT"; // DEFAULT_READ_CONTEXT_FLANK_SIZE
    public static final String TEST_RIGHT_FLANK = "CTGAGACTCA";
    public static final String TEST_LEFT_CORE = "AC"; // based on MIN_CORE_DISTANCE
    public static final String TEST_RIGHT_CORE = "GT";

    // Read context counter creation

    public static VariantReadContext createReadContext(int position, final String ref, final String alt)
    {
        return createReadContext(createSimpleVariant(position, ref, alt));
    }

    public static VariantReadContext createReadContext(final SimpleVariant variant)
    {
        return createReadContext(variant, TEST_LEFT_CORE, TEST_RIGHT_CORE);
    }

    public static VariantReadContext createReadContext(
            final SimpleVariant variant, final String leftCore, final String rightCore)
    {
        return createReadContext(variant, leftCore, rightCore, TEST_LEFT_FLANK, TEST_RIGHT_FLANK);
    }

    public static VariantReadContext createReadContext(
            final SimpleVariant variant, final String leftCore, final String rightCore, final String leftFlank, final String rightFlank)
    {
        // create a context with no repeat or homology
        int coreIndexStart = leftFlank.length();
        int varReadIndex = coreIndexStart + leftCore.length();
        int coreIndexEnd = varReadIndex + variant.alt().length() - 1 + rightCore.length();
        String refBases = leftFlank + leftCore + variant.ref() + rightCore + rightFlank;
        String readBases = leftFlank + leftCore + variant.alt() + rightCore + rightFlank;

        int alignmentStart = variant.Position - varReadIndex;

        int coreRightPosStart = variant.Position + min(variant.ref().length(), variant.alt().length());
        int alignmentEnd = coreRightPosStart + rightCore.length() + rightFlank.length() - 1;

        List<CigarElement> readCigar = List.of(new CigarElement(readBases.length(), CigarOperator.M));

        return new VariantReadContext(
                variant, alignmentStart, alignmentEnd, refBases.getBytes(), readBases.getBytes(), readCigar,
                coreIndexStart, varReadIndex, coreIndexEnd, null, null);
    }

    // Read context counter
    public static ReadContextCounter createReadCounter(final int id, final VariantReadContext readContext)
    {
        /* CLEAN-UP
        SimpleVariant variant = createSimpleVariant(position);

        IndexedBases indexBases = new IndexedBases(position, 10, "ACGTACGTACGT".getBytes());
        ReadContext readContext = new ReadContext(position, "", 0, "", indexBases, false);
        */

        return new ReadContextCounter(
                id, readContext, VariantTier.LOW_CONFIDENCE,
                100, 1, TEST_CONFIG, QUALITY_CALCULATOR, null);
    }

    // Sage variant creation
    public static SageVariant createSageVariant(int position, final String ref, final String alt)
    {
        SimpleVariant variant = createSimpleVariant(position, ref, alt);

        /* CLEAN-UP
        String readBases = buildReadContextBases(alt);

        // LF          L   I   RC          RF
        // 0123456789  01  2  34  0123456789
        int leftCoreIndex = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        int index = leftCoreIndex + MIN_CORE_DISTANCE;
        int rightCoreIndex = index + alt.length() - 1 + MIN_CORE_DISTANCE;

        IndexedBases indexBases = new IndexedBases(
                position, index, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases.getBytes());
        */

        VariantReadContext readContext = createReadContext(variant);
        return createSageVariant(readContext);
    }

    public static String buildReadContextBases(final String alt)
    {
        String flank = generateRandomBases(DEFAULT_FLANK_LENGTH);
        String core = generateRandomBases(MIN_CORE_DISTANCE);
        return flank + core + alt + core + flank;
    }

    public static SageVariant createSageVariant(final VariantReadContext readContext)
    {
        ReadContextCounter readCounter = new ReadContextCounter(
                0, readContext, VariantTier.LOW_CONFIDENCE, 100, 1,
                TestUtils.TEST_CONFIG, QUALITY_CALCULATOR, null);

        List<ReadContextCounter> tumorCounters = Lists.newArrayList(readCounter);

        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, tumorCounters.get(0).readContext(), 1, 1);

        List<ReadContextCounter> normalCounters = Lists.newArrayList();

        return new SageVariant(candidate, normalCounters, tumorCounters);
    }

    public static SageVariant sageVariantFromReadContextCounter(final ReadContextCounter readContextCounter)
    {
        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, readContextCounter.readContext(), 1, 1);

        return new SageVariant(candidate, Collections.emptyList(), Lists.newArrayList(readContextCounter));
    }


    // CLEAN-UP: old context and variant creation methods
    /*
    public static ReadContext createReadContext(
            int refPosition, int readIndex, int leftCentreIndex, int rightCentreIndex, String readBases, String microhomology)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length() - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBasesIndexed = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, 0, readBases.getBytes());

        return new ReadContext(refPosition, "", 0, microhomology, readBasesIndexed, incompleteCore);
    }

    @Deprecated
    public static SageVariant createVariantOld(
            final String chromosome, int position, final String ref, final String alt, final IndexedBases indexBases)
    {
        SimpleVariant variant = new SimpleVariant(chromosome, position, ref, alt);

        ReadContext readContext = new ReadContext(position, "", 0, "", indexBases, false);

        ReadContextCounter readCounter = new ReadContextCounter(
                0, null, VariantTier.LOW_CONFIDENCE, 100, 1,
                TestUtils.TEST_CONFIG, QUALITY_CALCULATOR, null);

        List<ReadContextCounter> tumorCounters = Lists.newArrayList(readCounter);

        Candidate candidate = new Candidate(
                VariantTier.HIGH_CONFIDENCE, tumorCounters.get(0).readContext(), 1, 1);

        List<ReadContextCounter> normalCounters = Lists.newArrayList();

        return new SageVariant(candidate, normalCounters, tumorCounters);
    }
    */
}
