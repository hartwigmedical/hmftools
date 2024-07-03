package com.hartwig.hmftools.purple.fitting;

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

    @VisibleForTesting
    public static boolean isPolyclonal(final double copyNumber)
    {
        double remainder = Math.abs(copyNumber - Math.round(copyNumber));
        return greaterThan(remainder, POLYCLONAL_DISTANCE);
    }
}