package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

class ExtendDiploidBAF {

    @NotNull
    static List<CombinedRegion> extendBAF(@NotNull final List<CombinedRegion> regions) {

        final Set<CombinedRegion> processed = Sets.newHashSet();

        int nextIndex = nextIndex(processed, regions);
        while (nextIndex != -1) {

            extendRight(nextIndex, regions);
            extendLeft(nextIndex, regions);

            processed.add(regions.get(nextIndex));
            nextIndex = nextIndex(processed, regions);
        }

        return regions;
    }

    private static void extendRight(int startIndex, @NotNull final List<CombinedRegion> regions) {
        assert (startIndex < regions.size());
        final CombinedRegion source = regions.get(startIndex);

        for (int i = startIndex + 1; i < regions.size(); i++) {
            final CombinedRegion target = regions.get(i);
            if (target.bafCount() > 0 || target.isInferredBAF()) {
                return;
            }
            target.setInferredTumorBAF(estimateBAF(target.tumorCopyNumber(), source.tumorBAF(), source.tumorCopyNumber()));
        }
    }

    private static void extendLeft(int startIndex, @NotNull final List<CombinedRegion> regions) {
        assert (startIndex < regions.size());
        final CombinedRegion source = regions.get(startIndex);

        for (int i = startIndex - 1; i >= 0; i--) {
            final CombinedRegion target = regions.get(i);
            if (target.bafCount() > 0 || target.isInferredBAF()) {
                return;
            }
            target.setInferredTumorBAF(estimateBAF(target.tumorCopyNumber(), source.tumorBAF(), source.tumorCopyNumber()));
        }
    }

    private static int nextIndex(@NotNull final Set<CombinedRegion> processed, @NotNull final List<CombinedRegion> regions) {

        int indexOfLargestBaf = -1;
        int biggestBaddestBAF = 0;

        for (int i = 0; i < regions.size(); i++) {
            final CombinedRegion combined = regions.get(i);
            if (combined.bafCount() > biggestBaddestBAF && !processed.contains(combined)) {
                biggestBaddestBAF = combined.bafCount();
                indexOfLargestBaf = i;
            }
        }

        return indexOfLargestBaf;
    }

    static double estimateBAF(double copyNumber, double neighbourBAF, double neighbourCopyNumber) {

        if (Doubles.isZero(copyNumber)) {
            return 0;
        }

        if (Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }

        double neighbourMajorAllele = Math.min(1d, neighbourBAF) * neighbourCopyNumber;

        double regionMajorAllele = Math.max(0, neighbourMajorAllele - neighbourCopyNumber + copyNumber);
        double regionBAF = regionMajorAllele / copyNumber;

        return 0.5 + Math.abs(regionBAF - 0.5);
    }
}
