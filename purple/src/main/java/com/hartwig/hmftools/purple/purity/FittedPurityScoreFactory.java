package com.hartwig.hmftools.purple.purity;

import static com.hartwig.hmftools.common.utils.Doubles.greaterThan;

import java.util.Collection;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurityScore;

import org.jetbrains.annotations.NotNull;

public final class FittedPurityScoreFactory
{
    private static final double POLYCLONAL_DISTANCE = 0.25;

    @NotNull
    public static FittedPurityScore score(@NotNull final List<FittedPurity> purities)
    {
        ImmutableFittedPurityScore.Builder builder = ImmutableFittedPurityScore.builder()
                .minPloidy(0)
                .minPurity(0)
                .minDiploidProportion(0)
                .maxPloidy(0)
                .maxPurity(0)
                .maxDiploidProportion(0);

        purities.stream().max(FittedPurityScoreFactory::comparePloidy).ifPresent(x -> builder.maxPloidy(x.ploidy()));
        purities.stream().min(FittedPurityScoreFactory::comparePloidy).ifPresent(x -> builder.minPloidy(x.ploidy()));
        purities.stream().max(FittedPurityScoreFactory::comparePurity).ifPresent(x -> builder.maxPurity(x.purity()));
        purities.stream().min(FittedPurityScoreFactory::comparePurity).ifPresent(x -> builder.minPurity(x.purity()));
        purities.stream()
                .max(FittedPurityScoreFactory::compareDiploidProportion)
                .ifPresent(x -> builder.maxDiploidProportion(x.diploidProportion()));
        purities.stream()
                .min(FittedPurityScoreFactory::compareDiploidProportion)
                .ifPresent(x -> builder.minDiploidProportion(x.diploidProportion()));

        return builder.build();
    }

    public static double polyclonalProportion(@NotNull final Collection<PurpleCopyNumber> regions)
    {
        int polyclonalCount = 0;
        int totalCount = 0;

        for(final PurpleCopyNumber region : regions)
        {
            totalCount += region.bafCount();
            if(isPolyclonal(region.averageTumorCopyNumber()))
            {
                polyclonalCount += region.bafCount();
            }
        }

        return totalCount == 0 ? 0 : 1d * polyclonalCount / totalCount;
    }

    private static int comparePurity(@NotNull final FittedPurity o1, @NotNull final FittedPurity o2)
    {
        return Double.compare(o1.purity(), o2.purity());
    }

    private static int comparePloidy(@NotNull final FittedPurity o1, @NotNull final FittedPurity o2)
    {
        return Double.compare(o1.ploidy(), o2.ploidy());
    }

    private static int compareDiploidProportion(@NotNull final FittedPurity o1, @NotNull final FittedPurity o2)
    {
        return Double.compare(o1.diploidProportion(), o2.diploidProportion());
    }

    @VisibleForTesting
    static boolean isPolyclonal(final double copyNumber)
    {
        double remainder = Math.abs(copyNumber - Math.round(copyNumber));
        return greaterThan(remainder, POLYCLONAL_DISTANCE);
    }
}