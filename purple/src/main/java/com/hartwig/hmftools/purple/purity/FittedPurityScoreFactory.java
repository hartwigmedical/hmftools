package com.hartwig.hmftools.purple.purity;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Doubles.greaterThan;

import java.util.Collection;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityScore;

import org.jetbrains.annotations.NotNull;

public final class FittedPurityScoreFactory
{
    private static final double POLYCLONAL_DISTANCE = 0.25;

    public static FittedPurityScore score(final List<FittedPurity> purities)
    {
        ImmutableFittedPurityScore.Builder builder = ImmutableFittedPurityScore.builder()
                .minPloidy(0)
                .minPurity(0)
                .minDiploidProportion(0)
                .maxPloidy(0)
                .maxPurity(0)
                .maxDiploidProportion(0);

        double maxPloidy = 0;
        double minPloidy = Double.MAX_VALUE;
        double maxPurity = 0;
        double minPurity = Double.MAX_VALUE;
        double maxDp = 0;
        double minDp = Double.MAX_VALUE;

        for(FittedPurity fittedPurity : purities)
        {
            maxPloidy = max(maxPloidy, fittedPurity.ploidy());
            minPloidy = min(minPloidy, fittedPurity.ploidy());

            maxPurity = max(maxPurity, fittedPurity.purity());
            minPurity = min(minPurity, fittedPurity.purity());

            maxDp = max(maxDp, fittedPurity.diploidProportion());
            minDp = min(minDp, fittedPurity.diploidProportion());
        }

        builder.maxPurity(maxPurity)
                .minPurity(minPurity)
                .maxPloidy(maxPloidy)
                .minPloidy(minPloidy)
                .maxDiploidProportion(maxDp)
                .minDiploidProportion(minDp);

        FittedPurityScore fittedPurityScore = builder.build();
        return fittedPurityScore;

        /*
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

        FittedPurityScore fittedPurityScore = builder.build();

        if(minPurity != fittedPurityScore.minPurity() || maxPurity != fittedPurityScore.maxPurity()
        || minPloidy != fittedPurityScore.minPloidy() || maxPloidy != fittedPurityScore.maxPloidy()
        || minDp != fittedPurityScore.minDiploidProportion() || maxDp != fittedPurityScore.maxDiploidProportion())
        {
            // assert(false);
        }
        */
    }

    public static double polyclonalProportion(final Collection<PurpleCopyNumber> regions)
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

    private static int comparePurity(final FittedPurity o1, final FittedPurity o2)
    {
        return Double.compare(o1.purity(), o2.purity());
    }

    private static int comparePloidy(final FittedPurity o1, final FittedPurity o2)
    {
        return Double.compare(o1.ploidy(), o2.ploidy());
    }

    private static int compareDiploidProportion(final FittedPurity o1, final FittedPurity o2)
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