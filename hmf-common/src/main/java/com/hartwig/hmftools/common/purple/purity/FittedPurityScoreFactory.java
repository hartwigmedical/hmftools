package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;
import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;

public class FittedPurityScoreFactory {

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0025;
    private static final double POLYCLONAL_DISTANCE = 0.25;

    @NotNull
    public static FittedPurityScore score(@NotNull final List<FittedPurity> purities,
            @NotNull final List<PurpleCopyNumber> regions) {
        ImmutableFittedPurityScore.Builder builder = ImmutableFittedPurityScore.builder()
                .minPloidy(0)
                .minPurity(0)
                .maxPloidy(0)
                .maxPurity(0)
                .polyclonalProportion(polyclonalProproption(regions));

        if (!purities.isEmpty()) {
            Collections.sort(purities);

            final FittedPurity best = purities.get(0);

            final List<FittedPurity> withinRange = purities.stream()
                    .filter(inRange(best.score()))
                    .collect(Collectors.toList());

            withinRange.stream()
                    .max(FittedPurityScoreFactory::comparePloidy)
                    .ifPresent(x -> builder.maxPloidy(x.ploidy()));
            withinRange.stream()
                    .min(FittedPurityScoreFactory::comparePloidy)
                    .ifPresent(x -> builder.minPloidy(x.ploidy()));
            withinRange.stream()
                    .max(FittedPurityScoreFactory::comparePurity)
                    .ifPresent(x -> builder.maxPurity(x.purity()));
            withinRange.stream()
                    .min(FittedPurityScoreFactory::comparePurity)
                    .ifPresent(x -> builder.minPurity(x.purity()));
        }

        return builder.build();
    }

    private static int comparePurity(@NotNull final FittedPurity o1, @NotNull final FittedPurity o2) {
        return Double.compare(o1.purity(), o2.purity());
    }

    private static int comparePloidy(@NotNull final FittedPurity o1, @NotNull final FittedPurity o2) {
        return Double.compare(o1.ploidy(), o2.ploidy());
    }

    @NotNull
    private static Predicate<FittedPurity> inRange(final double score) {
        return fittedPurity -> {
            double absDifference = Math.abs(fittedPurity.score() - score);
            double relDifference = Math.abs(absDifference / score);
            return lessOrEqual(absDifference, ABS_RANGE) || lessOrEqual(relDifference, PERCENT_RANGE);
        };
    }

    private static double polyclonalProproption(@NotNull final Collection<PurpleCopyNumber> regions) {
        int polyclonalCount = 0;
        int totalCount = 0;

        for (final PurpleCopyNumber region : regions) {
            totalCount += region.bafCount();
            if (isPolyclonal(region.averageTumorCopyNumber())) {
                polyclonalCount += region.bafCount();
            }
        }

        return totalCount == 0 ? 0 : 1d * polyclonalCount / totalCount;
    }

    @VisibleForTesting
    static boolean isPolyclonal(final double copyNumber) {
        double remainder = Math.abs(copyNumber - Math.round(copyNumber));
        return greaterThan(remainder, POLYCLONAL_DISTANCE);
    }
}
