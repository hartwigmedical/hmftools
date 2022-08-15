package com.hartwig.hmftools.purple.copynumber;

import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;

final class Extend
{
    static boolean doNotExtend(final CombinedRegion target, final ObservedRegion neighbour)
    {
        return breakForCentromereStart(target, neighbour) || breakForStructuralVariant(target, neighbour);
    }

    private static boolean breakForCentromereStart(final CombinedRegion target, final ObservedRegion neighbour)
    {
        if(target.start() < neighbour.start())
        {
            return neighbour.support() == SegmentSupport.CENTROMERE;
        }

        return target.region().support() == SegmentSupport.CENTROMERE;
    }

    private static boolean breakForStructuralVariant(final CombinedRegion target, final ObservedRegion neighbour)
    {
        if(target.start() < neighbour.start())
        {
            return neighbour.support().isSV();
        }

        return target.region().support().isSV();
    }
}
