package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;
import static com.hartwig.hmftools.common.numeric.Doubles.positiveOrZero;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;

public class FittedPurityFactory {

    private final int maxPloidy;
    private final double minPurity;
    private final double maxPurity;
    private final double purityIncrements;
    private final double minNormFactor;
    private final double maxNormFactor;
    private final double normFactorIncrements;
    private final FittedCopyNumberFactory fittedCopyNumberFactory;

    public FittedPurityFactory(int maxPloidy, double minPurity, double maxPurity, double purityIncrements, double minNormFactor,
                               double maxNormFactor, double normFactorIncrements, FittedCopyNumberFactory fittedCopyNumberFactory) {
        this.maxPloidy = maxPloidy;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
        this.purityIncrements = purityIncrements;
        this.minNormFactor = minNormFactor;
        this.maxNormFactor = maxNormFactor;
        this.normFactorIncrements = normFactorIncrements;
        this.fittedCopyNumberFactory = fittedCopyNumberFactory;
    }

    public List<FittedPurity> fitPurity(Collection<EnrichedCopyNumber> copyNumbers) {
        final List<FittedPurity> result = Lists.newArrayList();

        int totalBAFCount = 0;
        List<EnrichedCopyNumber> filteredCopyNumbers = Lists.newArrayList();
        for (EnrichedCopyNumber copyNumber : copyNumbers) {
            if (copyNumber.mBAFCount() > 0 && positiveOrZero(copyNumber.tumorRatio()) && Chromosomes.asInt(copyNumber.chromosome()) <= 22) {
                totalBAFCount += copyNumber.mBAFCount();
                filteredCopyNumbers.add(copyNumber);
            }
        }

        for (double purity = minPurity; lessOrEqual(purity, maxPurity); purity += purityIncrements) {
            for (double normFactor = minNormFactor; lessOrEqual(normFactor, maxNormFactor); normFactor += normFactorIncrements) {

                double impliedPloidy = impliedPloidy(normFactor, purity);

                if (greaterOrEqual(impliedPloidy, 1) && lessOrEqual(impliedPloidy, maxPloidy)) {
                    result.add(fitPurity(purity, normFactor, totalBAFCount, filteredCopyNumbers));
                }
            }
        }

        Collections.sort(result);
        return result;
    }

    private static double impliedPloidy(double normFactor, double purity) {
        return (1 - normFactor) / purity / normFactor * 2 + 2;
    }

    private FittedPurity fitPurity(double purity, double normFactor, double sumWeight, Collection<EnrichedCopyNumber> copyNumbers) {
        ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder().purity(purity).normFactor(normFactor);
        double modelDeviation = 0;
        double diploidProportion = 0;
        double modelBAFDeviation = 0;

        for (EnrichedCopyNumber copyNumber : copyNumbers) {

            final FittedCopyNumber fittedCopyNumber = fittedCopyNumberFactory.fittedCopyNumber(purity, normFactor, copyNumber);

            modelDeviation += copyNumber.mBAFCount() / sumWeight * fittedCopyNumber.deviation();
            modelBAFDeviation += copyNumber.mBAFCount() / sumWeight * fittedCopyNumber.bafDeviation();
            if (fittedCopyNumber.fittedPloidy() == 2) {
                diploidProportion += copyNumber.mBAFCount() / sumWeight;
            }
        }

        return builder
                .score(modelDeviation)
                .modelBAFDeviation(modelBAFDeviation)
                .diplodProportion(diploidProportion)
                .build();
    }

}
