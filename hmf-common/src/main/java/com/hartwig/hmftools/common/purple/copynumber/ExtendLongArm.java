package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

final class ExtendLongArm {

    @NotNull
    static List<CombinedRegion> extendLongArm(@NotNull final List<CombinedRegion> regions) {

        final int centromereIndex = findCentromere(regions);
        if (centromereIndex > 0) {
            final double copyNumber = regions.get(centromereIndex).tumorCopyNumber();
            extendLeft(copyNumber, centromereIndex - 1, regions);
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

    private static void extendLeft(double copyNumber, int targetIndex, @NotNull final List<CombinedRegion> regions) {
        if (targetIndex < 0 || regions.get(targetIndex).isProcessed()) {
            return;
        }

        final CombinedRegion target = regions.get(targetIndex);
        target.setTumorCopyNumber(CopyNumberMethod.LONG_ARM, copyNumber);
        if (target.region().support() != SegmentSupport.NONE) {
            extendLeft(copyNumber, targetIndex - 1, regions);
            return;
        }

        targetIndex--;
        while (targetIndex >= 0) {
            final CombinedRegion neighbour = regions.get(targetIndex);
            if (neighbour.isProcessed()) {
                return;
            }

            target.extend(neighbour.region());
            regions.remove(targetIndex);
            targetIndex--;

            if (target.region().support() != SegmentSupport.NONE) {
                extendLeft(copyNumber, targetIndex, regions);
                return;
            }
        }

    }
}
