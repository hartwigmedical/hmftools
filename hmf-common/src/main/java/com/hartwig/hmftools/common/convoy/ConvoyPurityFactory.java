package com.hartwig.hmftools.common.convoy;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import static com.hartwig.hmftools.common.numeric.Doubles.lessOrEqual;

public class ConvoyPurityFactory {

    public static List<ConvoyPurity> create(
            double minPurity, double maxPurity, double purityIncrements,
            double minNormFactor, double maxNormFactor, double normFactorIncrements,
            double cnvRatioWeighFactor, int maxPloidy, Collection<ConvoyCopyNumber> copyNumbers) {
        final List<ConvoyPurity> result = Lists.newArrayList();

        double sumWeight = copyNumbers.stream().mapToLong(ConvoyCopyNumber::mBAFCount).sum();

        for (double purity = minPurity; lessOrEqual(purity, maxPurity); purity += purityIncrements) {
            for (double normFactor = minNormFactor; lessOrEqual(normFactor, maxNormFactor); normFactor += normFactorIncrements) {
                result.add(create(purity, normFactor, cnvRatioWeighFactor, maxPloidy, sumWeight, copyNumbers));
            }
        }

        Collections.sort(result);
        return result;
    }

    private static ConvoyPurity create(double purity, double normFactor, double cnvRatioWeighFactor, int maxPloidy,
                double sumWeight, Collection<ConvoyCopyNumber> copyNumbers) {
        ImmutableConvoyPurity.Builder builder = ImmutableConvoyPurity.builder().purity(purity).normFactor(normFactor);
        double modelDeviation = 0;
        double diploidProportion = 0;
        double modelBAFDeviation = 0;

        for (ConvoyCopyNumber copyNumber : copyNumbers) {

            int fittedPloidy = 0;
            double minDeviation = 0;
            double minBafDeviation = 0;

            for (int ploidy = 1; ploidy <= maxPloidy; ploidy++) {
                double modelCNVRatio = modelCNVRatio(purity, normFactor, ploidy);
                double cnvDeviation = cnvDeviation(cnvRatioWeighFactor, modelCNVRatio, copyNumber.tumorRatio());
                double modelBAFToMinimizeDeviation = modelBAFToMinimizeDeviation(purity, ploidy, copyNumber.mBAF());

                double bafDeviation = Math.abs(copyNumber.mBAF() - modelBAFToMinimizeDeviation);
                double deviation = ploidy / 2.0 * (bafDeviation + cnvDeviation);

                if (ploidy == 1 || deviation < minDeviation) {
                    minDeviation = deviation;
                    minBafDeviation = bafDeviation;
                    fittedPloidy = ploidy;
                }
            }

            modelDeviation += copyNumber.mBAFCount() / sumWeight * minDeviation;
            modelBAFDeviation += copyNumber.mBAFCount() / sumWeight * minBafDeviation;
            if (fittedPloidy == 2) {
                diploidProportion += copyNumber.mBAFCount() / sumWeight;
            }
        }

        return builder
                .score(modelDeviation)
                .modelBAFDeviation(modelBAFDeviation)
                .diplodProportion(diploidProportion)
                .build();
    }

    @VisibleForTesting
    static double modelCNVRatio(double purity, double normFactor, int ploidy) {
        return normFactor + (ploidy - CopyNumber.NORMAL_HUMAN_COPY_NUMBER) * purity * normFactor / 2d;
    }

    @VisibleForTesting
    static double cnvDeviation(double cnvRatioWeighFactor, double modelCNVRatio, double actualRatio) {
        return cnvRatioWeighFactor * Math.abs(modelCNVRatio - actualRatio);
    }

    @VisibleForTesting
    static double modelBAFToMinimizeDeviation(double purity, int ploidy, double actualBAF) {
        double result = 0;
        double deviation = 0;

        int minBetaAllele = (int) Math.round(ploidy / 2d);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            double modelBAF = modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = Math.abs(modelBAF - actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < deviation) {
                result = modelBAF;
                deviation = modelDeviation;
            }
        }

        return result;
    }

    @VisibleForTesting
    static double modelBAF(double purity, int ploidy, int betaAllele) {
        if (ploidy / betaAllele == 2) {
            return 0.533; // TODO Check why this is with pete!
        }

        return (1 + purity * (betaAllele - 1)) / (2 + purity * (ploidy - 2));
    }
}
