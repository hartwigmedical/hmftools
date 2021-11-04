package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ModifiableEnrichedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendMinSupportTest
{

    @Test
    public void testDefault()
    {
        final ModifiableEnrichedRegion prior = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ModifiableEnrichedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior, target));
        assertEquals(1000, target.minStart());
    }

    @Test
    public void testExtendOnlyDiploid()
    {
        final ModifiableEnrichedRegion prior = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ModifiableEnrichedRegion target = create(2000, GermlineStatus.NOISE, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior, target));
        assertEquals(2000, target.minStart());
    }

    @Test
    public void testExtendOnlyNoSVSupport()
    {
        final ModifiableEnrichedRegion prior = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ModifiableEnrichedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.BND);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior, target));
        assertEquals(2000, target.minStart());
    }

    @Test
    public void testKeepExtending()
    {
        final ModifiableEnrichedRegion prior1 = create(500, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ModifiableEnrichedRegion prior2 = create(1000, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ModifiableEnrichedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior1, prior2, target));
        assertEquals(500, target.minStart());
    }

    @Test
    public void testStopExtendingAtSV()
    {
        final ModifiableEnrichedRegion prior1 = create(500, GermlineStatus.NOISE, SegmentSupport.NONE);
        final ModifiableEnrichedRegion prior2 = create(1000, GermlineStatus.NOISE, SegmentSupport.INS);
        final ModifiableEnrichedRegion target = create(2000, GermlineStatus.DIPLOID, SegmentSupport.NONE);

        assertEquals(2000, target.minStart());
        ObservedRegionFactory.extendMinSupport(Lists.newArrayList(prior1, prior2, target));
        assertEquals(1000, target.minStart());
    }

    @NotNull
    private ModifiableEnrichedRegion create(long start, @NotNull final GermlineStatus status, @NotNull final SegmentSupport support)
    {
        return ModifiableEnrichedRegion.create().setStart(start).setMinStart(start).setGermlineStatus(status).setSupport(support);
    }
}
