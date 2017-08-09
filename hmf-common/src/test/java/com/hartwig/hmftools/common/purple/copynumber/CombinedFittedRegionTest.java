package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CombinedFittedRegionTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void averageOnLengthUntilNonZeroBafCount() {
        final CombinedFittedRegion region = createCombinedFittedRegion(1, 100_000_000, 3);
        assertAverages(region, 0, 3);

        region.combine(create(100_000_001, 200_000_000, 4));
        assertAverages(region, 0, 3.5);

        region.combine(create(200_000_001, 200_000_010, 1, 0.5, 3));
        assertAverages(region, 0.5, 3);

        region.combine(create(200_000_011, 300_000_000, 3, 1, 4d));
        assertAverages(region, 0.875, 3.75);
    }

    @Test
    public void averageOnLengthForNonZeroRatio() {
        CombinedFittedRegion builder = createCombinedFittedRegion(1, 100, 3);
        assertAverages(builder, 0, 3);

        builder.combine(create(101, 200, 0));
        assertAverages(builder, 0, 3);
    }

    @Test
    public void doNotIncludeZeroCopyNumber() {
        final FittedRegion startRegion = create(1, 100, 200, 0.5, 0);
        CombinedFittedRegion builder = new CombinedFittedRegion(true, startRegion);
        assertAverages(builder, 0.5, 0);

        builder.combine(create(201, 300, 200, 1, 2));
        assertAverages(builder, 0.75, 2);
    }

    private static void assertAverages(@NotNull CombinedFittedRegion victim, double expectedBAF, double expectedCopyNumber) {
        assertAverages(victim.region(), expectedBAF, expectedCopyNumber);
    }

    private static void assertAverages(@NotNull FittedRegion victim, double expectedBAF, double expectedCopyNumber) {
        assertEquals(expectedBAF, victim.observedBAF(), EPSILON);
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    private CombinedFittedRegion createCombinedFittedRegion(long start, long end, double copyNumber) {
        return new CombinedFittedRegion(true, create(start, end, copyNumber));
    }

    private static FittedRegion create(long start, long end, double copyNumber) {
        return create("1", start, end, 0, 0, copyNumber);
    }

    private static FittedRegion create(long start, long end, int bafCount, double baf, double copyNumber) {
        return create("1", start, end, bafCount, baf, copyNumber);
    }

    @NotNull
    private static FittedRegion create(@NotNull String chromosome, long start, long end, int bafCount, double baf, double tumorCopyNumber) {
        return PurpleDatamodelTest.createDefaultFittedRegion(chromosome, start, end)
                .bafCount(bafCount)
                .observedBAF(baf)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(tumorCopyNumber)
                .build();
    }
}
