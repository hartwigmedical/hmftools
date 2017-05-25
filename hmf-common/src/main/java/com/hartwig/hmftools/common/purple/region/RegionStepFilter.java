package com.hartwig.hmftools.common.purple.region;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.signum;

import static com.hartwig.hmftools.common.numeric.Doubles.equal;
import static com.hartwig.hmftools.common.numeric.Doubles.lessThan;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

enum MergeRegionBreaks {
    ;

    private static final long MIN_SIZE = 1000;

    @NotNull
    static List<ConsolidatedRegion> merge(@NotNull final List<ConsolidatedRegion> region) {
        List<ConsolidatedRegion> results = Lists.newArrayList();

        int i = 0;
        while (i < region.size()) {
            ConsolidatedRegion first = region.get(i);
            ConsolidatedRegion second = i + 1 < region.size() ? region.get(i + 1) : null;
            ConsolidatedRegion third = i + 2 < region.size() ? region.get(i + 2) : null;

            if (second != null && third != null && second.bases() == MIN_SIZE && isSameChromosome(first, third)) {
                double myFirstSecondDifference = first.averageTumorCopyNumber() - second.averageTumorCopyNumber();
                double mySecondThirdDifference = second.averageTumorCopyNumber() - third.averageTumorCopyNumber();

                if (equal(signum(myFirstSecondDifference), signum(mySecondThirdDifference))) {
                    if (lessThan(abs(myFirstSecondDifference), abs(mySecondThirdDifference))) {
                        results.add(merge(first, second));
                        i += 2;

                    } else {
                        results.add(first);
                        results.add(merge(third, second));
                        i += 3;
                    }
                    continue;
                }
            }

            results.add(first);
            i++;
        }

        return results;
    }

    @NotNull
    private static ConsolidatedRegion merge(@NotNull final ConsolidatedRegion primary,
            @NotNull final ConsolidatedRegion secondary) {
        return ImmutableConsolidatedRegion.builder().from(primary).start(min(primary.start(), secondary.start())).end(
                max(primary.end(), secondary.end())).build();
    }

    private static boolean isSameChromosome(@NotNull final GenomeRegion region1, @NotNull final GenomeRegion region2) {
        return region1.chromosome().equals(region2.chromosome());
    }
}
