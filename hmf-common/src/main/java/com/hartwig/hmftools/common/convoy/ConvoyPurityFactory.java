package com.hartwig.hmftools.common.convoy;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.copynumber.CopyNumber;

import java.util.Collection;

public class ConvoyPurityFactory {

    public static ConvoyPurity create(double purity, double normFactory, double cnvRatioWeighFactor, Collection<ConvoyCopyNumber> copyNumbers) {

        double sumWeight = copyNumbers.stream().mapToLong(ConvoyCopyNumber::mBAFCount).sum();
        return create(purity, normFactory, cnvRatioWeighFactor, 10, sumWeight, copyNumbers);
    }

    static ConvoyPurity create(double purity, double normFactor, double cnvRatioWeighFactor, int maxPloidy, double sumWeight, Collection<ConvoyCopyNumber> copyNumbers) {

        ImmutableConvoyPurity.Builder builder = ImmutableConvoyPurity.builder().purity(purity).normFactor(normFactor);
        double modelDeviation = 0, diploidProportion = 0, modelBAFDeviation = 0;

        for (ConvoyCopyNumber copyNumber : copyNumbers) {

            if (copyNumber.mBAFCount() > 0 && copyNumber.tumorRatio() >= 0) {

                int fittedPloidy = 0;
                double minDeviation = 0;
                double minBafDeviation = 0;

                for (int ploidy = 1; ploidy <= maxPloidy; ploidy++) {
                    double modelCNVRatio = modelCNVRatio(purity, normFactor, ploidy);
                    double cnvDeviation = cnvDeviation(cnvRatioWeighFactor, ploidy, modelCNVRatio, copyNumber.ratioOfRatio()); // TODO: switch to tumor ratio?
                    double modelBAFToMinimizeDeviation = modelBAFToMinimizeDeviation(purity, ploidy, copyNumber.mBAF());

                    double bafDeviation = Math.abs(copyNumber.mBAF() - modelBAFToMinimizeDeviation);
                    double deviation = bafDeviation + cnvDeviation;

                    if (ploidy == 1 || deviation < minDeviation) {
                        minDeviation = deviation;
                        minBafDeviation = bafDeviation;
                        fittedPloidy = ploidy;
                    }
                }

                modelDeviation += copyNumber.mBAFCount() / sumWeight * minDeviation;
                modelBAFDeviation += copyNumber.mBAFCount() / sumWeight * minBafDeviation;
                if (fittedPloidy == 2) {
                    diploidProportion += diploidProportion + copyNumber.mBAFCount() / sumWeight;
                }
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
    static double cnvDeviation(double cnvRatioWeighFactor, int ploidy, double modelCNVRatio, double actualRatio) {
        return ploidy / 2.0 * cnvRatioWeighFactor * Math.abs(modelCNVRatio - actualRatio);
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
