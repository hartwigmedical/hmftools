package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.function.BiPredicate;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

@Deprecated
class CombinedFittedRegions {

    @NotNull
    public static List<CombinedRegion> mergeLeft(@NotNull final List<CombinedRegion> regions,
            @NotNull final BiPredicate<CombinedRegion, CombinedRegion> merge) {
        final List<CombinedRegion> result = Lists.newArrayList();

        if (!regions.isEmpty()) {
            CombinedRegion left = regions.get(0);
            for (int i = 1; i < regions.size(); i++) {
                CombinedRegion right = regions.get(i);
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
