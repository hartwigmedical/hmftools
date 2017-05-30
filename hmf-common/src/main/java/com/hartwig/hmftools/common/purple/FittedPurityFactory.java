package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.purity.PurityAdjustment;

import org.jetbrains.annotations.NotNull;

public class FittedPurityFactory {

    private final int maxPloidy;
    private final double minPurity;
    private final double maxPurity;
    private final double purityIncrements;
    private final double minNormFactor;
    private final double maxNormFactor;
    private final double normFactorIncrements;
    @NotNull
    private final FittedRegionFactory fittedRegionFactory;

    public FittedPurityFactory(final int maxPloidy, final double minPurity, final double maxPurity,
            final double purityIncrements, final double minNormFactor, final double maxNormFactor,
            final double normFactorIncrements, @NotNull final FittedRegionFactory fittedRegionFactory) {
        this.maxPloidy = maxPloidy;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
        this.purityIncrements = purityIncrements;
        this.minNormFactor = minNormFactor;
        this.maxNormFactor = maxNormFactor;
        this.normFactorIncrements = normFactorIncrements;
        this.fittedRegionFactory = fittedRegionFactory;
    }

    @NotNull
    public List<FittedPurity> fitPurity(@NotNull final Collection<EnrichedRegion> copyNumbers) {
        final List<FittedPurity> result = Lists.newArrayList();

        int totalBAFCount = 0;
        final List<EnrichedRegion> filteredCopyNumbers = Lists.newArrayList();
        for (final EnrichedRegion copyNumber : copyNumbers) {
            if (copyNumber.mBAFCount() > 0 && positiveOrZero(copyNumber.tumorRatio())
                    && Chromosomes.asInt(copyNumber.chromosome()) <= 22) {
                totalBAFCount += copyNumber.mBAFCount();
                filteredCopyNumbers.add(copyNumber);
            }
        }

        for (double purity = minPurity; lessOrEqual(purity, maxPurity); purity += purityIncrements) {
            for (double normFactor = minNormFactor; lessOrEqual(normFactor,
                    maxNormFactor); normFactor += normFactorIncrements) {

                double impliedPloidy = PurityAdjustment.purityAdjustedCopynumber(purity, normFactor, 1);

                if (greaterOrEqual(impliedPloidy, 1) && lessOrEqual(impliedPloidy, maxPloidy)) {
                    result.add(fitPurity(purity, normFactor, totalBAFCount, filteredCopyNumbers));
                }
            }
        }

        Collections.sort(result);
        return result;
    }

    @NotNull
    private FittedPurity fitPurity(double purity, double normFactor, double sumWeight,
            Collection<EnrichedRegion> enrichedRegions) {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double modelDeviation = 0;
        double diploidProportion = 0;
        double modelBAFDeviation = 0;

        for (EnrichedRegion enrichedRegion : enrichedRegions) {
            final FittedRegion fittedRegion = fittedRegionFactory.fittedCopyNumber(purity, normFactor,
                    enrichedRegion);
            modelDeviation += enrichedRegion.mBAFCount() / sumWeight * fittedRegion.deviation();
            modelBAFDeviation += enrichedRegion.mBAFCount() / sumWeight * fittedRegion.bafDeviation();
            if (fittedRegion.fittedPloidy() == 2) {
                diploidProportion += enrichedRegion.mBAFCount() / sumWeight;
            }
        }

        return builder.score(modelDeviation).modelBAFDeviation(modelBAFDeviation).diploidProportion(
                diploidProportion).build();
    }
}
