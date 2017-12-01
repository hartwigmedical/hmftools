package com.hartwig.hmftools.common.purple.copynumber;

import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

class Extend {

    static boolean doNotExtend(@NotNull final CombinedRegion target, @NotNull final FittedRegion neighbour) {
        return breakForCentromereStart(target, neighbour) || breakForStructuralVariant(target, neighbour);
    }

    static boolean breakForCentromereStart(@NotNull final CombinedRegion target, @NotNull final FittedRegion neighbour) {
        if (target.start() < neighbour.start()) {
            return neighbour.support() == SegmentSupport.CENTROMERE;
        }

        return target.region().support() == SegmentSupport.CENTROMERE;
    }

    static boolean breakForStructuralVariant(@NotNull final CombinedRegion target, @NotNull final FittedRegion neighbour) {
        if (target.start() < neighbour.start()) {
            return neighbour.support() != SegmentSupport.NONE;
        }

        return target.region().support() != SegmentSupport.NONE;
    }

}
