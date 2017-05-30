package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.purity.PurityAdjustment;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFactory;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;

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
    public List<FittedPurity> fitPurity(@NotNull final Collection<ObservedRegion> observedRegions) {
        final List<FittedPurity> result = Lists.newArrayList();

        int totalBAFCount = 0;
        final List<ObservedRegion> filteredCopyNumbers = Lists.newArrayList();
        for (final ObservedRegion copyNumber : observedRegions) {
            if (copyNumber.bafCount() > 0 && positiveOrZero(copyNumber.observedTumorRatio())
                    && Chromosomes.asInt(copyNumber.chromosome()) <= 22) {
                totalBAFCount += copyNumber.bafCount();
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
            Collection<ObservedRegion> observedRegions) {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double modelDeviation = 0;
        double diploidProportion = 0;
        double modelBAFDeviation = 0;

        for (ObservedRegion enrichedRegion : observedRegions) {
            final FittedRegion fittedRegion = fittedRegionFactory.fitRegion(purity, normFactor,
                    enrichedRegion);
            modelDeviation += enrichedRegion.bafCount() / sumWeight * fittedRegion.deviation();
            modelBAFDeviation += enrichedRegion.bafCount() / sumWeight * fittedRegion.bafDeviation();
            if (fittedRegion.fittedPloidy() == 2) {
                diploidProportion += enrichedRegion.bafCount() / sumWeight;
            }
        }

        return builder.score(modelDeviation).modelBAFDeviation(modelBAFDeviation).diploidProportion(
                diploidProportion).build();
    }
}
