package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.TestUtils.createDefaultFittedRegion;
import static com.hartwig.hmftools.purple.copynumber.ExtendLongArmTest.assertCombinedRegion;
import static com.hartwig.hmftools.purple.copynumber.ExtendLongArmTest.createCombinedRegion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.purple.region.ImmutableFittedRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendNonDiploidTest
{

    @Test
    public void testIsEligible()
    {
        final FittedRegion eligible = createFittedRegion(GermlineStatus.HET_DELETION, SegmentSupport.BND);
        final FittedRegion centromere = ImmutableFittedRegion.builder().from(eligible).support(SegmentSupport.CENTROMERE).build();

        assertTrue(ExtendNonDiploid.isEligible(eligible, null));
        assertTrue(ExtendNonDiploid.isEligible(eligible, eligible));
        assertTrue(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder()
                .from(eligible)
                .germlineStatus(GermlineStatus.AMPLIFICATION)
                .build(), null));
        assertTrue(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder()
                .from(eligible)
                .germlineStatus(GermlineStatus.HOM_DELETION)
                .build(), null));
        assertTrue(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder()
                .from(eligible)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .build(), null));

        assertTrue(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder()
                .from(eligible)
                .support(SegmentSupport.NONE)
                .build(), eligible));

        assertFalse(ExtendNonDiploid.isEligible(eligible, centromere));
        assertFalse(ExtendNonDiploid.isEligible(centromere, eligible));

        assertFalse(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder().from(eligible).germlineStatus(GermlineStatus.NOISE).build(), null));
        assertFalse(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder().from(eligible).support(SegmentSupport.NONE).build(), null));
        assertFalse(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder()
                .from(eligible)
                .germlineStatus(GermlineStatus.DIPLOID)
                .build(), null));
        assertFalse(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder().from(eligible).depthWindowCount(0).build(), null));
        assertFalse(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder().from(eligible).observedTumorRatio(0.9).build(), null));
        assertFalse(ExtendNonDiploid.isEligible(ImmutableFittedRegion.builder().from(eligible).observedNormalRatio(1.2).build(), null));
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
    private static FittedRegion createFittedRegion(GermlineStatus status, SegmentSupport support)
    {
        return createDefaultFittedRegion("1", 1, 1000)
                .depthWindowCount(1)
                .observedNormalRatio(1)
                .observedTumorRatio(1.1)
                .germlineStatus(status)
                .support(support)
                .build();
    }
}
