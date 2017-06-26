package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.purple.copynumber.BaseCopyNumberBuilder.allowedCopyNumberDeviation;
import static com.hartwig.hmftools.common.purple.copynumber.HighConfidenceCopyNumberBuilder.isEven;
import static com.hartwig.hmftools.common.purple.copynumber.HighConfidenceCopyNumberBuilder.purityAdjustedBAF;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HighConfidencePurpleCopyNumberBuilderTest {
    private static final double EPSILON = 1e-10;

    @Test
    public void copyNumberAllowances() {
        assertEquals(1.3, allowedCopyNumberDeviation(0), EPSILON);
        assertEquals(0.8, allowedCopyNumberDeviation(5), EPSILON);
        assertEquals(0.3, allowedCopyNumberDeviation(10), EPSILON);
        assertEquals(0.3, allowedCopyNumberDeviation(50), EPSILON);
    }


    @Test
    public void testIsEvenCopyNumber() {
        assertFalse(isEven(0.5));
        assertFalse(isEven(0.75));
        assertFalse(isEven(1));
        assertFalse(isEven(1.74));

        assertTrue(isEven(1.75));
        assertTrue(isEven(2));
        assertTrue(isEven(2.25));

        assertFalse(isEven(2.5));
        assertFalse(isEven(2.75));
        assertFalse(isEven(3));
        assertFalse(isEven(3.74));

        assertTrue(isEven(3.75));
        assertTrue(isEven(4));
        assertTrue(isEven(4.25));
    }

    @Test
    public void testPurityAdjustedBaf() {
        testPurityAdjustedBaf(0.1);
        testPurityAdjustedBaf(0.2);
        testPurityAdjustedBaf(0.3);
        testPurityAdjustedBaf(0.4);
        testPurityAdjustedBaf(0.5);
        testPurityAdjustedBaf(0.6);
        testPurityAdjustedBaf(0.7);
        testPurityAdjustedBaf(0.8);
        testPurityAdjustedBaf(0.9);
    }

    private void testPurityAdjustedBaf(double purity) {
        testPurityAdjustedBaf(purity, 2, 1);
        testPurityAdjustedBaf(purity, 2, 2);
        testPurityAdjustedBaf(purity, 3, 2);
        testPurityAdjustedBaf(purity, 3, 3);
        testPurityAdjustedBaf(purity, 4, 2);
        testPurityAdjustedBaf(purity, 4, 3);
        testPurityAdjustedBaf(purity, 4, 4);
        testPurityAdjustedBaf(purity, 5, 3);
        testPurityAdjustedBaf(purity, 5, 4);
    }

    private static void testPurityAdjustedBaf(final double purity, final int ploidy, final int alleleCount) {
        double expectedPurityAdjustedBAF = 1d * alleleCount / ploidy;
        double observedBAF = modelBAF(purity, ploidy, alleleCount);
        assertEquals(expectedPurityAdjustedBAF, purityAdjustedBAF(purity, ploidy, observedBAF), EPSILON);
    }

    @Test
    public void averageOnLengthUntilNonZeroBafCount() {
        final HighConfidenceCopyNumberBuilder builder = new HighConfidenceCopyNumberBuilder(1, create(1, 100_000_000, 3));
        assertAverages(builder, 0, 3);

        builder.extendRegion(create(100_000_001, 200_000_000, 4));
        assertAverages(builder, 0, 3.5);

        builder.extendRegion(create(200_000_001, 200_000_010, 1, 0.5, 3));
        assertAverages(builder, 0.5, 3);

        builder.extendRegion(create(200_000_011, 300_000_000, 3, 1, 4d));
        assertAverages(builder, 0.875, 3.75);
    }

    @Test
    public void averageOnLengthForNonZeroRatio() {
        HighConfidenceCopyNumberBuilder builder = new HighConfidenceCopyNumberBuilder(1, create(1, 100, 3));
        assertAverages(builder, 0, 3);

        builder.extendRegion(create(101, 200, 0));
        assertAverages(builder, 0, 3);
    }

    @Test
    public void doNotIncludeZeroRatio() {
        final FittedRegion startRegion = create(1, 100, 200, 0.5, 0);
        HighConfidenceCopyNumberBuilder builder = new HighConfidenceCopyNumberBuilder(1, startRegion);
        assertAverages(builder, 0.5, 0);

        builder.extendRegion(create(201, 300, 200, 1, 2));
        assertAverages(builder, 0.75, 2);
    }

    private static void assertAverages(@NotNull HighConfidenceCopyNumberBuilder victim, double expectedBAF,
            double expectedRatio) {
        assertAverages(victim.build(), expectedBAF, expectedRatio);
    }

    private static void assertAverages(@NotNull PurpleCopyNumber victim, double expectedBAF, double expectedRatio) {
        assertEquals(expectedBAF, victim.averageObservedBAF(), EPSILON);
        assertEquals(expectedRatio, victim.averageTumorCopyNumber(), EPSILON);
    }

    private static FittedRegion create(long start, long end, double copyNumber) {
        return create("1", start, end, 0, 0, copyNumber);
    }

    private static FittedRegion create(long start, long end, int bafCount, double baf, double copyNumber) {
        return create("1", start, end, bafCount, baf, copyNumber);
    }

    @NotNull
    public static FittedRegion create(@NotNull String chromosome, long start, long end, int bafCount, double baf,
            double tumorCopyNumber) {
        return ImmutableFittedRegion.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .bafCount(bafCount)
                .observedBAF(baf)
                .tumorCopyNumber(tumorCopyNumber)
                .broadBAF(0)
                .broadTumorCopyNumber(0)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .observedNormalRatio(1.0)
                .observedNormalRatio(1.0)
                .cnvDeviation(0)
                .deviation(0)
                .fittedPloidy(0)
                .modelBAF(0)
                .observedTumorRatio(0)
                .modelTumorRatio(0)
                .status(FreecStatus.UNKNOWN)
                .refNormalisedCopyNumber(tumorCopyNumber)
                .ratioSupport(true)
                .structuralVariantSupport(StructuralVariantSupport.NONE)
                .bafDeviation(0)
                .build();
    }

    private static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        assert (alleleCount >= ploidy / 2d);
        return (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2));
    }
}
