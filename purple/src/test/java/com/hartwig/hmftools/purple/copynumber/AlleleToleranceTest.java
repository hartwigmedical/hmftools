package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.MiscTestUtils.createDefaultFittedRegion;
import static com.hartwig.hmftools.purple.MiscTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.junit.Before;
import org.junit.Test;

public class AlleleToleranceTest
{
    private AlleleTolerance victim;

    @Before
    public void setup()
    {
        victim = new AlleleTolerance(buildPurityAdjuster(Gender.FEMALE, 1, 1));
    }

    @Test
    public void testMinObservedBafDeviation()
    {
        ObservedRegion left = createBafRegion(1000, 1, 0);
        ObservedRegion right = createBafRegion(1000, 0.5, 0);
        assertTrue(victim.inTolerance(left, right));

        left = ObservedRegion.from(left);
        left.setObservedBAF(0.031);
        assertFalse(victim.inTolerance(left, right));
    }

    @Test
    public void testMinBafCount()
    {
        ObservedRegion left = createBafRegion(0, 1, 0.031);
        ObservedRegion right = createBafRegion(1000, 0.5, 0);
        assertTrue(victim.inTolerance(left, right));

        left = ObservedRegion.from(left);
        left.setBafCount(1000);
        assertFalse(victim.inTolerance(left, right));
    }

    @Test
    public void testRelativeChange()
    {
        ObservedRegion left = createCopyNumberRegion(1000, 10, 10);
        ObservedRegion right = createCopyNumberRegion(1000, 9.1, 9.1);
        assertTrue(victim.inTolerance(left, right));

        left = createCopyNumberRegion(1000, 2, 2);
        right = createCopyNumberRegion(1000, 1.1, 1.1);
        assertFalse(victim.inTolerance(left, right));
    }

    @Test
    public void testRelativeChangeToNegativeCopyNumber()
    {
        assertEquals(6, AlleleTolerance.relativeCopyNumberChange(2.5, -0.5), 0.01);
    }

    private static ObservedRegion createBafRegion(int bafCount, double tumorBaf, double observedBaf)
    {
        ObservedRegion region = createDefaultFittedRegion("1", 100, 200);
        region.setBafCount(bafCount);
        region.setTumorBAF(tumorBaf);
        region.setObservedBAF(observedBaf);
        region.setTumorCopyNumber(2);
        region.setRefNormalisedCopyNumber(2);
        region.setDepthWindowCount(100);
        return region;
    }

    private static ObservedRegion createCopyNumberRegion(int depthWindowCount, double copyNumber, double refCopyNumber)
    {
        ObservedRegion region = createDefaultFittedRegion("1", 100, 200);
        region.setBafCount(0);
        region.setTumorCopyNumber(copyNumber);
        region.setRefNormalisedCopyNumber(refCopyNumber);
        region.setDepthWindowCount(depthWindowCount);
        return region;
    }
}
