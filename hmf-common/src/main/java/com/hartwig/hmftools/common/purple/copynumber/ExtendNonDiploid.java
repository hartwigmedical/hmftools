package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ExtendNonDiploid {

    @NotNull
    static List<CombinedRegion> nonDiploid(@NotNull final List<CombinedRegion> regions) {

        for (int i = 0; i < regions.size(); i++) {
            CombinedRegion region = regions.get(i);
            FittedRegion next = i < regions.size() - 1 ? regions.get(i + 1).region() : null;
            if (region.copyNumberMethod().equals(CopyNumberMethod.UNKNOWN) && isEligible(region.region(), next)) {
                region.setTumorCopyNumber(CopyNumberMethod.NON_DIPLOID, region.region().refNormalisedCopyNumber());
            }
        }
        return extendNonDiploid(regions);
    }

    static boolean isEligible(@NotNull final FittedRegion region, @Nullable final FittedRegion neighbour) {
        return Doubles.greaterThan(region.observedTumorRatio(), region.observedNormalRatio()) && !region.status()
                .equals(GermlineStatus.DIPLOID) && !region.status().equals(GermlineStatus.NOISE) && region.depthWindowCount() > 0
                && !isBoundByCentromere(region, neighbour) && isBoundBySV(region, neighbour);
    }

    private static boolean isBoundByCentromere(@NotNull final FittedRegion region, @Nullable final FittedRegion neighbour) {
        return region.support().equals(SegmentSupport.CENTROMERE) || (neighbour != null && neighbour.support()
                .equals(SegmentSupport.CENTROMERE));
    }

    private static boolean isBoundBySV(@NotNull final FittedRegion region, @Nullable final FittedRegion neighbour) {
        return region.support().isSV() || (neighbour != null && neighbour.support().isSV());
    }

    @NotNull
    static List<CombinedRegion> extendNonDiploid(@NotNull final List<CombinedRegion> regions) {

        for (int i = 0; i < regions.size(); i++) {

            CombinedRegion region = regions.get(i);
            if (region.copyNumberMethod() == CopyNumberMethod.NON_DIPLOID) {
                extendRight(i, regions);
                i -= extendLeft(i, regions);
            }
        }

        return regions;
    }

    private static void extendRight(int startIndex, @NotNull final List<CombinedRegion> regions) {
        assert (startIndex < regions.size());
        final CombinedRegion target = regions.get(startIndex);
        int targetIndex = startIndex + 1;

        while (targetIndex < regions.size()) {
            final CombinedRegion neighbour = regions.get(targetIndex);

            if (Extend.doNotExtend(target, neighbour.region())) {
                break;
            }

            if (neighbour.copyNumberMethod() == CopyNumberMethod.NON_DIPLOID) {
                target.extendWithWeightedAverage(neighbour.region());
            } else {
                break;
            }

            regions.remove(targetIndex);
        }
    }

    private static int extendLeft(int startIndex, @NotNull final List<CombinedRegion> regions) {
        assert (startIndex < regions.size());
        final CombinedRegion target = regions.get(startIndex);

        int targetIndex = startIndex - 1;
        while (targetIndex >= 0) {
            final CombinedRegion neighbour = regions.get(targetIndex);
            if (Extend.doNotExtend(target, neighbour.region())) {
                break;
            }

            if (neighbour.copyNumberMethod() == CopyNumberMethod.NON_DIPLOID) {
                target.extendWithWeightedAverage(neighbour.region());
            } else {
                break;
            }

            regions.remove(targetIndex);
            targetIndex--;
        }

        return startIndex - targetIndex - 1;
    }

}
