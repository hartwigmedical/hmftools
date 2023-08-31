package com.hartwig.hmftools.purple.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltTestUtils;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;

public class PopulateUnknownTest
{
    @Test
    public void testJoinUnknownRegions()
    {
        final CobaltChromosomes female = CobaltTestUtils.female();

        final CombinedRegion region1 = ExtendLongArmTest.createCombinedRegion(1, 1000, 0, 0, SegmentSupport.NONE);
        final CombinedRegion region2 = ExtendLongArmTest.createCombinedRegion(1001, 2000, 0, 0, SegmentSupport.NONE);
        final CombinedRegion region3 = ExtendLongArmTest.createCombinedRegion(2001, 3000, 0, 0, SegmentSupport.NONE);
        final CombinedRegion region4 = ExtendLongArmTest.createCombinedRegion(3001, 4000, 0, 0, SegmentSupport.NONE);
        region4.setTumorCopyNumber(CopyNumberMethod.BAF_WEIGHTED, 2);

        List<CombinedRegion> result = ExtendUtils.populateUnknown(Lists.newArrayList(region1, region2, region3, region4), female);
        assertEquals(2, result.size());
    }
}
