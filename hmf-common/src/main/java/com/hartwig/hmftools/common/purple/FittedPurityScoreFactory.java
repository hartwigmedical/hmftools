package com.hartwig.hmftools.common.purple;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.ConsolidatedRegion;

import org.jetbrains.annotations.NotNull;

public class FittedPurityScoreFactory {

    private static final double PERCENT_RANGE = 0.1;
    private static final double POLYCLONAL_DISTANCE = 0.25;

    @NotNull
    public static FittedPurityScore score(@NotNull final List<FittedPurity> purities, @NotNull final List<ConsolidatedRegion> regions) {
        ImmutableFittedPurityScore.Builder builder = ImmutableFittedPurityScore.builder()
                .minPloidy(0)
                .minPurity(0)
                .maxPloidy(0)
                .maxPurity(0)
                .polyclonalProportion(polyclonalProproption(regions));

        if (!purities.isEmpty()) {
            Collections.sort(purities);

            FittedPurity best = purities.get(0);

            List<FittedPurity> withinPercent = purities.stream()
                    .filter(inPercentRange(PERCENT_RANGE, best.score()))
                    .collect(Collectors.toList());

            withinPercent.stream()
                    .max(FittedPurityScoreFactory::comparePloidy)
                    .ifPresent(x -> builder.maxPloidy(x.ploidy()));
            withinPercent.stream()
                    .min(FittedPurityScoreFactory::comparePloidy)
                    .ifPresent(x -> builder.minPloidy(x.ploidy()));
            withinPercent.stream()
                    .max(FittedPurityScoreFactory::comparePurity)
                    .ifPresent(x -> builder.maxPurity(x.purity()));
            withinPercent.stream()
                    .min(FittedPurityScoreFactory::comparePurity)
                    .ifPresent(x -> builder.minPurity(x.purity()));
        }

        return builder.build();
    }

    private static int comparePurity(@NotNull FittedPurity o1, @NotNull FittedPurity o2) {
        return Double.compare(o1.purity(), o2.purity());
    }

    private static int comparePloidy(@NotNull FittedPurity o1, @NotNull FittedPurity o2) {
        return Double.compare(o1.ploidy(), o2.ploidy());
    }

    @NotNull
    private static Predicate<FittedPurity> inPercentRange(double percent, double score) {
        return fittedPurity -> Doubles.lessOrEqual(Math.abs((fittedPurity.score() - score) / score), percent);
    }

    private static double polyclonalProproption(@NotNull Collection<ConsolidatedRegion> regions) {
        int polyclonalCount = 0;
        int totalCount = 0;

        for (ConsolidatedRegion region : regions) {
            totalCount += region.bafCount();
            if (isPolyclonal(region.averageTumorCopyNumber())) {
                polyclonalCount += region.bafCount();
            }
        }

        return totalCount == 0 ? 0 : 1d * polyclonalCount / totalCount;
    }

    @VisibleForTesting
    static boolean isPolyclonal(double copyNumber) {
        double remainder = Math.abs(copyNumber - Math.round(copyNumber));
        return Doubles.greaterThan(remainder, POLYCLONAL_DISTANCE);
    }
}
