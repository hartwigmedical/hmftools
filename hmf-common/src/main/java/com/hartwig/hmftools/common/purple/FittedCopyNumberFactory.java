package com.hartwig.hmftools.common.purple;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;

public class FittedCopyNumberFactory {

    static final double NORMAL_BAF = 0.533;

    private final int maxPloidy;
    private final double cnvRatioWeightFactor;

    public FittedCopyNumberFactory(int maxPloidy, double cnvRatioWeightFactor) {
        this.maxPloidy = maxPloidy;
        this.cnvRatioWeightFactor = cnvRatioWeightFactor;
    }

    public List<FittedCopyNumber> fittedCopyNumber(double purity, double normFactor,
            Collection<EnrichedCopyNumber> copyNumbers) {
        return copyNumbers.stream().map(x -> fittedCopyNumber(purity, normFactor, x)).collect(Collectors.toList());
    }

    FittedCopyNumber fittedCopyNumber(double purity, double normFactor, EnrichedCopyNumber copyNumber) {

        double minDeviation = 0;
        double observedBAF = copyNumber.mBAF();
        double observedTumorRatio = copyNumber.tumorRatio();
        double tumorCopyNumber = copyNumber(purity, normFactor, observedTumorRatio);

        ImmutableFittedCopyNumber.Builder builder = ImmutableFittedCopyNumber.builder()
                .from(copyNumber)
                .status(FreecStatus.fromNormalRatio(copyNumber.normalRatio()))
                .bafCount(copyNumber.mBAFCount())
                .observedBAF(observedBAF)
                .observedTumorRatio(observedTumorRatio)
                .observedNormalRatio(copyNumber.normalRatio())
                .purityAdjustedBAF(purityAdjustedBAF(purity, observedBAF))
                .broadBAF(0)
                .broadTumorCopyNumber(0)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(Doubles.replaceNaNWithZero(observedTumorRatio / copyNumber.normalRatio() / normFactor * 2))
                .value((int) Math.round(tumorCopyNumber));

        for (int ploidy = 1; ploidy <= maxPloidy; ploidy++) {
            double modelRatio = modelRatio(purity, normFactor, ploidy);
            double cnvDeviation = cnvDeviation(cnvRatioWeightFactor, modelRatio, observedTumorRatio);

            double modelBAF =
                    copyNumber.mBAFCount() == 0 ? 0 : modelBAFToMinimizeDeviation(purity, ploidy, observedBAF);
            double bafDeviation = bafDeviation(modelBAF, observedBAF);

            double deviation =
                    Math.pow(Math.max(ploidy, 1.5) / 2.0, 0.85) * (bafDeviation + cnvDeviation) * observedBAF;

            if (ploidy == 1 || deviation < minDeviation) {
                builder.fittedPloidy(ploidy)
                        .modelBAF(modelBAF)
                        .modelTumorRatio(modelRatio)
                        .bafDeviation(bafDeviation)
                        .cnvDeviation(cnvDeviation)
                        .deviation(deviation);
                minDeviation = deviation;
            }
        }

        return builder.build();
    }

    @VisibleForTesting
    static double modelRatio(double purity, double normFactor, int ploidy) {
        return normFactor + (ploidy - 2) * purity * normFactor / 2d;
    }

    @VisibleForTesting
    static double copyNumber(double purity, double normFactor, double tumorRatio) {
        return 2 + 2 * (tumorRatio - normFactor) / purity / normFactor;
    }

    @VisibleForTesting
    static double cnvDeviation(double cnvRatioWeighFactor, double modelCNVRatio, double actualRatio) {
        return cnvRatioWeighFactor * Math.abs(modelCNVRatio - actualRatio);
    }

    @VisibleForTesting
    static double bafDeviation(double modelBAF, double actualBAF) {
        if (Doubles.equal(modelBAF, NORMAL_BAF) && Doubles.lessOrEqual(actualBAF, NORMAL_BAF)) {
            return 0;
        }

        return Math.abs(modelBAF - actualBAF);
    }

    @VisibleForTesting
    static double modelBAFToMinimizeDeviation(double purity, int ploidy, double actualBAF) {
        double result = 0;
        double deviation = 0;

        int minBetaAllele = (int) Math.round(ploidy / 2d);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            double modelBAF = modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = bafDeviation(modelBAF, actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < deviation) {
                result = modelBAF;
                deviation = modelDeviation;
            }
        }

        return result;
    }

    @VisibleForTesting
    static double modelBAF(double purity, int ploidy, int alleleCount) {
        if (ploidy / alleleCount == 2) {
            return NORMAL_BAF;
        }

        return (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2));
    }

    @VisibleForTesting
    static double purityAdjustedBAF(double purity, double observedBAF) {
        assert (Doubles.greaterThan(purity, 0));
        return (observedBAF - (1 - purity) * NORMAL_BAF) / purity;
    }

}
