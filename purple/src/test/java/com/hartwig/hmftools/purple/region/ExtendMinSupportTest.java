package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendMinSupportTest
{
    @Test
    public void testDefault()
    {
        final ObservedRegion prior = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ObservedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior, target));
        assertEquals(1000, target.minStart());
    }

    @Test
    public void testExtendOnlyDiploid()
    {
        final ObservedRegion prior = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ObservedRegion target = create(2000, GermlineStatus.NOISE, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior, target));
        assertEquals(2000, target.minStart());
    }

    @Test
    public void testExtendOnlyNoSVSupport()
    {
        final ObservedRegion prior = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ObservedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.BND);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior, target));
        assertEquals(2000, target.minStart());
    }

    @Test
    public void testKeepExtending()
    {
        final ObservedRegion prior1 = create(500, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ObservedRegion prior2 = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ObservedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior1, prior2, target));
        assertEquals(500, target.minStart());
    }

    @Test
    public void testStopExtendingAtSV()
    {
        final ObservedRegion prior1 = create(500, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ObservedRegion prior2 = create(1000, GermlineStatus.NOISE, SegmentSupport.INS);
        final ObservedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior1, prior2, target));
        assertEquals(1000, target.minStart());
    }

    @NotNull
    private ObservedRegion create(int start, @NotNull final GermlineStatus status, @NotNull final SegmentSupport support)
    {
        return new ObservedRegion("", start, 0, false, support, 0, 0,
                0, 0, 0, 0, status, false,
                0, start, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0);
    }
}
