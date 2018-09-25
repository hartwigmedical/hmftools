package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.combine.CombinedRegion;

import org.jetbrains.annotations.NotNull;

class ExtendStructuralVariant {

    @NotNull
    static List<CombinedRegion> extendStructuralVariants(@NotNull final List<CombinedRegion> regions) {

        for (int i = 0; i < regions.size(); i++) {

            CombinedRegion region = regions.get(i);
            if (region.copyNumberMethod() == CopyNumberMethod.STRUCTURAL_VARIANT) {
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

            if (!neighbour.isProcessed()) {
                target.extend(neighbour.region());
            } else if (neighbour.copyNumberMethod() == CopyNumberMethod.STRUCTURAL_VARIANT) {
                target.extendWithUnweightedAverage(neighbour.region());
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

            if (!neighbour.isProcessed()) {
                target.extend(neighbour.region());
            } else if (neighbour.copyNumberMethod() == CopyNumberMethod.STRUCTURAL_VARIANT) {
                target.extendWithUnweightedAverage(neighbour.region());
            } else {
                break;
            }

            regions.remove(targetIndex);
            targetIndex--;
        }

        return startIndex - targetIndex - 1;
    }

}
