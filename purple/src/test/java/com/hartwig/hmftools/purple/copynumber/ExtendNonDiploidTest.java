package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.MiscTestUtils.createDefaultFittedRegion;
import static com.hartwig.hmftools.purple.copynumber.ExtendLongArmTest.assertCombinedRegion;
import static com.hartwig.hmftools.purple.copynumber.ExtendLongArmTest.createCombinedRegion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendNonDiploidTest
{
    @Test
    public void testIsEligible()
    {
        final ObservedRegion eligible = createFittedRegion(GermlineStatus.HET_DELETION, SegmentSupport.BND);
        final ObservedRegion centromere = ObservedRegion.from(eligible);
        centromere.setSupport(SegmentSupport.CENTROMERE);

        assertTrue(ExtendNonDiploid.isEligible(eligible, null));
        assertTrue(ExtendNonDiploid.isEligible(eligible, eligible));

        ObservedRegion newRegion = ObservedRegion.from(eligible);
        newRegion.setGermlineStatus(GermlineStatus.AMPLIFICATION);
        assertTrue(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion.setGermlineStatus(GermlineStatus.HOM_DELETION);
        assertTrue(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion.setGermlineStatus(GermlineStatus.UNKNOWN);
        assertTrue(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion.setGermlineStatus(GermlineStatus.HOM_DELETION);
        newRegion.setSupport(SegmentSupport.NONE);
        assertTrue(ExtendNonDiploid.isEligible(newRegion, eligible));

        assertFalse(ExtendNonDiploid.isEligible(eligible, centromere));
        assertFalse(ExtendNonDiploid.isEligible(centromere, eligible));

        newRegion = ObservedRegion.from(eligible);
        newRegion.setGermlineStatus(GermlineStatus.NOISE);
        assertFalse(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion = ObservedRegion.from(eligible);
        newRegion.setSupport(SegmentSupport.NONE);
        assertFalse(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion = ObservedRegion.from(eligible);
        newRegion.setGermlineStatus(GermlineStatus.DIPLOID);
        assertFalse(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion = ObservedRegion.from(eligible);
        newRegion.setDepthWindowCount(0);
        assertFalse(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion = ObservedRegion.from(eligible);
        newRegion.setObservedTumorRatio(0.9);
        assertFalse(ExtendNonDiploid.isEligible(newRegion, null));

        newRegion = ObservedRegion.from(eligible);
        newRegion.setObservedNormalRatio(1.2);
        assertFalse(ExtendNonDiploid.isEligible(newRegion, null));
    }

    @Test
    public void testSingleExtendLeft()
    {
        final CombinedRegion leftMost = createNonDiploidImplied(2001, 3000, 3, 0.5, SegmentSupport.BND);
        final CombinedRegion middle1 = createNonDiploidImplied(3001, 4000, 4, 0.6, SegmentSupport.NONE);
        final CombinedRegion middle2 = createNonDiploidImplied(4001, 5000, 4, 0.6, SegmentSupport.NONE);
        final CombinedRegion middle3 = createNonDiploidImplied(5001, 6000, 6, 0.9, SegmentSupport.NONE);
        final CombinedRegion rightMost = createCombinedRegion(6001, 7000, 7, 1, SegmentSupport.BND);

        final List<CombinedRegion> regions = Lists.newArrayList(leftMost, middle1, middle2, middle3, rightMost);
        final List<CombinedRegion> result = ExtendNonDiploid.nonDiploid(regions);

        assertEquals(2, result.size());
        assertCombinedRegion(2001, 6000, 4.25, CopyNumberMethod.NON_DIPLOID, result.get(0));
        assertCombinedRegion(6001, 7000, 7, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testSingleExtendRight()
    {
        final CombinedRegion leftMost = createNonDiploidImplied(2001, 3000, 3, 0.5, SegmentSupport.BND);
        final CombinedRegion middle1 = createNonDiploidImplied(3001, 4000, 4, 0.6, SegmentSupport.NONE);
        final CombinedRegion middle2 = createNonDiploidImplied(4001, 5000, 4, 0.6, SegmentSupport.NONE);
        final CombinedRegion middle3 = createNonDiploidImplied(5001, 6000, 6, 0.9, SegmentSupport.NONE);
        final CombinedRegion rightMost = createCombinedRegion(6001, 7000, 7, 1, SegmentSupport.BND);

        final List<CombinedRegion> regions = Lists.newArrayList(leftMost, middle1, middle2, middle3, rightMost);
        final List<CombinedRegion> result = ExtendNonDiploid.nonDiploid(regions);

        assertEquals(2, result.size());
        assertCombinedRegion(2001, 6000, 4.25, CopyNumberMethod.NON_DIPLOID, result.get(0));
        assertCombinedRegion(6001, 7000, 7, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testTwoSVsExtendBetween()
    {
        final CombinedRegion leftMost = createNonDiploidImplied(2001, 3000, 3, 0.5, SegmentSupport.BND);
        final CombinedRegion middle1 = createNonDiploidImplied(3001, 4000, 4, 0.6, SegmentSupport.NONE);
        final CombinedRegion middle2 = createNonDiploidImplied(4001, 5000, 5, 0.6, SegmentSupport.NONE);
        final CombinedRegion middle3 = createCombinedRegion(5001, 6000, 6, 0.9, SegmentSupport.NONE);
        final CombinedRegion rightMost = createCombinedRegion(6001, 7000, 7, 1, SegmentSupport.BND);

        final List<CombinedRegion> regions = Lists.newArrayList(leftMost, middle1, middle2, middle3, rightMost);
        final List<CombinedRegion> result = ExtendNonDiploid.nonDiploid(regions);

        assertEquals(2, result.size());
        assertCombinedRegion(2001, 6000, 4, CopyNumberMethod.NON_DIPLOID, result.get(0));
    }

    @NotNull
    private static CombinedRegion createNonDiploidImplied(int start, int end, double copyNumber, double baf, SegmentSupport support)
    {
        final CombinedRegion region = createCombinedRegion(start, end, copyNumber, baf, support);
        region.setCopyNumberMethod(CopyNumberMethod.NON_DIPLOID);
        return region;
    }

    @NotNull
    private static ObservedRegion createFittedRegion(GermlineStatus status, SegmentSupport support)
    {
        ObservedRegion region = createDefaultFittedRegion("1", 1, 1000);
        region.setDepthWindowCount(1);
        region.setObservedNormalRatio(1);
        region.setObservedTumorRatio(1.1);
        region.setGermlineStatus(status);
        region.setSupport(support);
        return region;
    }
}
