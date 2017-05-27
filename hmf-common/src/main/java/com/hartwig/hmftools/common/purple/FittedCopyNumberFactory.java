package com.hartwig.hmftools.common.purple;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;

public class FittedCopyNumberFactory {

    static final double NORMAL_BAF = 0.533;

    private final int maxPloidy;
    private final double cnvRatioWeightFactor;

    public FittedCopyNumberFactory(int maxPloidy, double cnvRatioWeightFactor) {
        this.maxPloidy = maxPloidy;
        this.cnvRatioWeightFactor = cnvRatioWeightFactor;
    }

    @NotNull
    public List<FittedCopyNumber> fittedCopyNumber(final double purity, final double normFactor,
            @NotNull final Collection<EnrichedRegion> copyNumbers) {
        return copyNumbers.stream().map(x -> fittedCopyNumber(purity, normFactor, x)).collect(Collectors.toList());
    }

    FittedCopyNumber fittedCopyNumber(double purity, double normFactor, EnrichedRegion copyNumber) {
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
                .purityAdjustedBAF(purityAdjustedBAF(purity, 0, observedBAF))
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

            double[] modelBAFWithDeviation = copyNumber.mBAFCount() == 0
                    ? new double[] { 0, 0 }
                    : modelBAFToMinimizeDeviation(purity, ploidy, observedBAF);

            double modelBAF = modelBAFWithDeviation[0];
            double bafDeviation = modelBAFWithDeviation[1];

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

    private static double copyNumber(double purity, double normFactor, double tumorRatio) {
        return 2 + 2 * (tumorRatio - normFactor) / purity / normFactor;
    }

    @VisibleForTesting
    static double cnvDeviation(double cnvRatioWeighFactor, double modelCNVRatio, double actualRatio) {
        return cnvRatioWeighFactor * Math.abs(modelCNVRatio - actualRatio);
    }

    @VisibleForTesting
    static double bafDeviation(boolean heterozygous, double modelBAF, double actualBAF) {
        if (heterozygous && Doubles.lessOrEqual(actualBAF, NORMAL_BAF)) {
            return 0;
        }

        return Math.abs(modelBAF - actualBAF);
    }

    @VisibleForTesting
    static double[] modelBAFToMinimizeDeviation(final double purity, final int ploidy, final double actualBAF) {
        double result = 0;
        double deviation = 0;

        int minBetaAllele = (int) Math.round(ploidy / 2d);
        for (int betaAllele = minBetaAllele; betaAllele < ploidy + 1; betaAllele++) {

            boolean isHetrozygous = ploidy / betaAllele == 2;
            double modelBAF = isHetrozygous ? NORMAL_BAF : modelBAF(purity, ploidy, betaAllele);
            double modelDeviation = bafDeviation(isHetrozygous, modelBAF, actualBAF);

            if (betaAllele == minBetaAllele || modelDeviation < deviation) {
                result = modelBAF;
                deviation = modelDeviation;
            }
        }

        return new double[]{result, deviation};
    }

    @VisibleForTesting
    static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        if (ploidy / alleleCount == 2) {
            return NORMAL_BAF;
        }

        return (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2));
    }

    static double purityAdjustedBAF(final double purity, final double ploidy, final double observedBAF) {
        assert (Doubles.greaterThan(purity, 0));
        return (observedBAF - (1 - purity) * NORMAL_BAF) / purity;
    }
}
