package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import org.jetbrains.annotations.NotNull;

abstract class ExtendRegion {

    private final CopyNumberMethod method;

    ExtendRegion(final CopyNumberMethod method) {
        this.method = method;
    }

    @NotNull
    public List<CombinedRegion> extend(@NotNull final List<CombinedRegion> regions) {

        for (int i = 0; i < regions.size(); i++) {

            CombinedRegion region = regions.get(i);
            if (region.copyNumberMethod().equals(method)) {
                extendRight(i, regions);
                i -= extendLeft(i, regions);
            }
        }

        return regions;
    }

    protected abstract void extend(final CombinedRegion target, final CombinedRegion neighbour);

    private void extendRight(int startIndex, @NotNull final List<CombinedRegion> regions) {
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
            } else if (neighbour.copyNumberMethod().equals(method)) {
                extend(target, neighbour);
            } else {
                break;
            }

            regions.remove(targetIndex);
        }
    }

    private int extendLeft(int startIndex, @NotNull final List<CombinedRegion> regions) {
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
            } else if (neighbour.copyNumberMethod().equals(method)) {
                extend(target, neighbour);
            } else {
                break;
            }

            regions.remove(targetIndex);
            targetIndex--;
        }

        return startIndex - targetIndex - 1;
    }
}
