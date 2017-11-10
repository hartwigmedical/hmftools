package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.function.BiPredicate;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

@Deprecated
class CombinedFittedRegions {

    @NotNull
    public static List<CombinedFittedRegion> mergeLeft(@NotNull final List<CombinedFittedRegion> regions,
            @NotNull final BiPredicate<CombinedFittedRegion, CombinedFittedRegion> merge) {
        final List<CombinedFittedRegion> result = Lists.newArrayList();

        if (!regions.isEmpty()) {
            CombinedFittedRegion left = regions.get(0);
            for (int i = 1; i < regions.size(); i++) {
                CombinedFittedRegion right = regions.get(i);
                if (merge.test(left, right)) {
                    left.combine(right.region());
                } else {
                    result.add(left);
                    left = right;
                }
            }
            result.add(left);
        }

        return result;
    }

}
