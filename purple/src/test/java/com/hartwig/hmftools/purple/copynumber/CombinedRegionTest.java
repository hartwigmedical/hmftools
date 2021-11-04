package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createDefaultFittedRegion;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.purple.copynumber.CombinedRegion;
import com.hartwig.hmftools.purple.copynumber.CombinedRegionImpl;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class CombinedRegionTest
{

    private static final double EPSILON = 1e-10;

    @Ignore
    @Test
    public void averageOnLengthUntilNonZeroBafCount()
    {
        final CombinedRegion region = createCombinedFittedRegion(1, 100_000_000, 3);
        assertAverages(region, 0, 3);

        region.extendWithWeightedAverage(create(100_000_001, 200_000_000, 4));
        assertAverages(region, 0, 3.5);

        region.extendWithWeightedAverage(create(200_000_001, 200_000_010, 1, 0.5, 3));
        assertAverages(region, 0.5, 3);

        region.extendWithWeightedAverage(create(200_000_011, 300_000_000, 3, 1, 4d));
        assertAverages(region, 0.875, 3.75);
    }

    @Test
    public void testDepthWindowCountSummationOnlyAppliesToSomatic()
    {
        final FittedRegion somaticRegion = createDefaultFittedRegion("1", 2001, 3000)
                .depthWindowCount(2)
                .germlineStatus(GermlineStatus.DIPLOID)
                .build();
        final CombinedRegion region = new CombinedRegionImpl(somaticRegion);
        assertEquals(2, region.region().depthWindowCount());

        final FittedRegion amplificationRegion = createDefaultFittedRegion("1", 1, 1000)
                .depthWindowCount(2)
                .germlineStatus(GermlineStatus.AMPLIFICATION)
                .build();

        region.extend(amplificationRegion);
        assertEquals(2, region.region().depthWindowCount());

        final FittedRegion germlineRegion = createDefaultFittedRegion("1", 1001, 2000)
                .depthWindowCount(2)
                .germlineStatus(GermlineStatus.AMPLIFICATION)
                .build();
        region.extend(germlineRegion);
        assertEquals(2, region.region().depthWindowCount());
    }

    @Test
    public void averageOnLengthForNonZeroRatio()
    {
        CombinedRegion builder = createCombinedFittedRegion(1, 100, 3);
        assertAverages(builder, 0, 3);

        builder.extendWithWeightedAverage(create(101, 200, 0));
        assertAverages(builder, 0, 3);
    }

    @Test
    public void doNotIncludeZeroCopyNumber()
    {
        final FittedRegion startRegion = create(1, 100, 200, 0.5, 0);
        CombinedRegion builder = new CombinedRegionImpl(startRegion);
        assertAverages(builder, 0.5, 0);

        builder.extendWithWeightedAverage(create(201, 300, 200, 1, 2));
        assertAverages(builder, 0.75, 2);
    }

    private static void assertAverages(@NotNull CombinedRegion victim, double expectedBAF, double expectedCopyNumber)
    {
        assertAverages(victim.region(), expectedBAF, expectedCopyNumber);
    }

    private static void assertAverages(@NotNull FittedRegion victim, double expectedBAF, double expectedCopyNumber)
    {
        assertEquals(expectedBAF, victim.observedBAF(), EPSILON);
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    @NotNull
    private static CombinedRegion createCombinedFittedRegion(long start, long end, double copyNumber)
    {
        return new CombinedRegionImpl(create(start, end, copyNumber));
    }

    @NotNull
    private static FittedRegion create(long start, long end, double copyNumber)
    {
        return create("1", start, end, 0, 0, copyNumber);
    }

    @NotNull
    private static FittedRegion create(long start, long end, int bafCount, double baf, double copyNumber)
    {
        return create("1", start, end, bafCount, baf, copyNumber);
    }

    @NotNull
    private static FittedRegion create(@NotNull String chromosome, long start, long end, int bafCount, double baf, double tumorCopyNumber)
    {
        return createDefaultFittedRegion(chromosome, start, end)
                .bafCount(bafCount)
                .observedBAF(baf)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(tumorCopyNumber)
                .build();
    }
}
