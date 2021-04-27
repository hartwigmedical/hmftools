package com.hartwig.hmftools.purple.copynumber;

import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

final class Extend
{
    static boolean doNotExtend(@NotNull final CombinedRegion target, @NotNull final FittedRegion neighbour)
    {
        return breakForCentromereStart(target, neighbour) || breakForStructuralVariant(target, neighbour);
    }

    private static boolean breakForCentromereStart(@NotNull final CombinedRegion target, @NotNull final FittedRegion neighbour)
    {
        if(target.start() < neighbour.start())
        {
            return neighbour.support() == SegmentSupport.CENTROMERE;
        }

        return target.region().support() == SegmentSupport.CENTROMERE;
    }

    private static boolean breakForStructuralVariant(@NotNull final CombinedRegion target, @NotNull final FittedRegion neighbour)
    {
        if(target.start() < neighbour.start())
        {
            return neighbour.support().isSV();
        }

        return target.region().support().isSV();
    }
}
