package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.MiscTestUtils.createDefaultFittedRegion;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;

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
        final ObservedRegion somaticRegion = createDefaultFittedRegion("1", 2001, 3000);
        somaticRegion.setDepthWindowCount(2);
        somaticRegion.setGermlineStatus(GermlineStatus.DIPLOID);

        final CombinedRegion region = new CombinedRegion(somaticRegion);
        assertEquals(2, region.region().depthWindowCount());

        final ObservedRegion amplificationRegion = createDefaultFittedRegion("1", 1, 1000);
        amplificationRegion.setDepthWindowCount(2);
        amplificationRegion.setGermlineStatus(GermlineStatus.AMPLIFICATION);

        region.extend(amplificationRegion);
        assertEquals(2, region.region().depthWindowCount());

        final ObservedRegion germlineRegion = createDefaultFittedRegion("1", 1001, 2000);
        germlineRegion.setDepthWindowCount(2);
        germlineRegion.setGermlineStatus(GermlineStatus.AMPLIFICATION);
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
        final ObservedRegion startRegion = create(1, 100, 200, 0.5, 0);
        CombinedRegion builder = new CombinedRegion(startRegion);
        assertAverages(builder, 0.5, 0);

        builder.extendWithWeightedAverage(create(201, 300, 200, 1, 2));
        assertAverages(builder, 0.75, 2);
    }

    private static void assertAverages(@NotNull CombinedRegion victim, double expectedBAF, double expectedCopyNumber)
    {
        assertAverages(victim.region(), expectedBAF, expectedCopyNumber);
    }

    private static void assertAverages(@NotNull ObservedRegion victim, double expectedBAF, double expectedCopyNumber)
    {
        assertEquals(expectedBAF, victim.observedBAF(), EPSILON);
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    @NotNull
    private static CombinedRegion createCombinedFittedRegion(int start, int end, double copyNumber)
    {
        return new CombinedRegion(create(start, end, copyNumber));
    }

    @NotNull
    private static ObservedRegion create(int start, int end, double copyNumber)
    {
        return create("1", start, end, 0, 0, copyNumber);
    }

    @NotNull
    private static ObservedRegion create(int start, int end, int bafCount, double baf, double copyNumber)
    {
        return create("1", start, end, bafCount, baf, copyNumber);
    }

    @NotNull
    private static ObservedRegion create(@NotNull String chromosome, int start, int end, int bafCount, double baf, double tumorCopyNumber)
    {
        ObservedRegion region = createDefaultFittedRegion(chromosome, start, end);
        region.setBafCount(bafCount);
        region.setObservedBAF(baf);
        region.setTumorCopyNumber(tumorCopyNumber);
        region.setRefNormalisedCopyNumber(tumorCopyNumber);
        return region;
    }
}
