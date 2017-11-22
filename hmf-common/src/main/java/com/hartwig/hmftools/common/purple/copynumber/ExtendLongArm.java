package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

class ExtendLongArm {

    @NotNull
    static List<CombinedRegion> extendLongArm(@NotNull final List<CombinedRegion> regions) {

        final int centromereIndex = findCentromere(regions);
        if (centromereIndex > 1) {
            extendLeft(centromereIndex, regions);
        }

        return regions;
    }

    private static int findCentromere(@NotNull final List<CombinedRegion> regions) {
        for (int i = 0; i < regions.size(); i++) {
            if (regions.get(i).region().support() == SegmentSupport.CENTROMERE) {
                return i;
            }
        }

        return -1;
    }

    private static void extendLeft(int startIndex, @NotNull final List<CombinedRegion> regions) {
        assert (startIndex < regions.size());
        if (startIndex > 1) {
            final CombinedRegion source = regions.get(startIndex);
            final CombinedRegion target = regions.get(startIndex - 1);
            if (!target.isProcessed()) {
                target.setTumorCopyNumber(CombinedRegionMethod.LONG_ARM, source.tumorCopyNumber());
            } else {
                return;
            }

            int targetIndex = startIndex - 2;
            while (targetIndex >= 0) {
                final CombinedRegion neighbour = regions.get(targetIndex);
                if (neighbour.isProcessed()) {
                    return;
                } else {
                    target.extend(regions.get(targetIndex).region());
                }

                regions.remove(targetIndex);
                targetIndex--;
            }
        }

    }

}
